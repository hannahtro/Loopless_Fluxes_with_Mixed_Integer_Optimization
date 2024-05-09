using SparseArrays, LinearAlgebra
using JuMP, HiGHS, SCIP
using MathOptInterface
using COBREXA
using DataFrames, JSON, DelimitedFiles

const MOI = MathOptInterface

include("optimization_model.jl")
include("conic_matrix.jl")
include("cobrexa.jl")
include("disjunctive_programming.jl")

# variable and consraint report
function variable_report(m, xi)
    return (
        name = name(xi),
        lower_bound = has_lower_bound(xi) ? lower_bound(xi) : -Inf,
        value = value(xi),
        upper_bound = has_upper_bound(xi) ? upper_bound(xi) : Inf,
        reduced_cost = reduced_cost(xi),
        obj_coefficient = coefficient(objective_function(m), xi),
        allowed_decrease = report[xi][1],
        allowed_increase = report[xi][2],
    )
end

function constraint_report(m, c::ConstraintRef)
    return (
        name = name(c),
        value = value(c),
        rhs = normalized_rhs(c),
        slack = normalized_rhs(c) - value(c),
        shadow_price = shadow_price(c),
        allowed_decrease = report[c][1],
        allowed_increase = report[c][2],
    )
end

function extract_basis(m)
    var_vec = all_variables(m)
    cs = all_constraints(m, include_variable_in_set_constraints=false)
    var_mapping = Dict([(val, idx) for (idx, val) in enumerate(var_vec)])

    # basic_cons = filter(cs) do c
    #     MOI.get(m, MOI.ConstraintBasisStatus(), c) == MOI.BASIC
    # end
    
    nonbasic_cons = filter(cs) do c
        MOI.get(m, MOI.ConstraintBasisStatus(), c) != MOI.BASIC
    end
    nonbasic_cons_idxs = [idx for (idx, val) in enumerate(cs) if val in nonbasic_cons]
    
    basic_vars = filter(var_vec) do v
        MOI.get(m, MOI.VariableBasisStatus(), v) == MOI.BASIC
    end
    
    # nonbasic_vars = filter(var_vec) do v
    #     MOI.get(m, MOI.VariableBasisStatus(), v) == MOI.NONBASIC
    # end
    # nonbasic_vars_idxs = [var_mapping[var_name] for var_name in nonbasic_vars]

    nonbasic_lower_vars = filter(var_vec) do v
        MOI.get(m, MOI.VariableBasisStatus(), v) == MOI.NONBASIC_AT_LOWER
    end
    
    nonbasic_upper_vars = filter(var_vec) do v
        MOI.get(m, MOI.VariableBasisStatus(), v) == MOI.NONBASIC_AT_UPPER
    end

    super_basic_vars = filter(var_vec) do v
        MOI.get(m, MOI.VariableBasisStatus(), v) == MOI.SUPER_BASIC
    end
    @show basic_vars
    # @show nonbasic_lower_vars, nonbasic_upper_vars, super_basic_vars

    # ensure that results match JuMP basis extraction
    conic_form = conic_matrix(m)
    @assert Matrix(conic_form.A) == Matrix(lp_matrix_data(m).A)
    @assert conic_form.b == lp_matrix_data(m).b_upper
    solution = value.(all_variables(m))

    # extract indices of active bounds
    # only works if bounds set with set_lower_bound()
    # vars_at_lower_bound_idxs = []
    # vars_at_upper_bound_idxs = []
    # for (idx,var) in enumerate(var_vec)
    #     # @show row_idx, length(cs)
    #     if has_lower_bound(var)
    #         if value(var) == lower_bound(var)
    #             push!(vars_at_lower_bound_idxs, (idx, lower_bound(var)))
    #         end
    #     end 
    #     if has_upper_bound(var)
    #         if value(var) == upper_bound(var)
    #             push!(vars_at_upper_bound_idxs, (idx, upper_bound(var)))
    #         end
    #     end
    # end 
    # @show [val for (idx, val) in enumerate(cs) if val in nonbasic_cons]

    vars_at_lower_bound_cs = []
    vars_at_upper_bound_cs = []
    for (idx,c) in enumerate(cs)
        if !(c in nonbasic_cons)
            if MOI.get(m, MOI.ConstraintName(), c) == "lbs"
                var = constraint_object(c).func.terms.keys[1]
                var_idx = var_mapping[var]
                coeff = constraint_object(c).func.terms[var]
                bound = constraint_object(c).set.upper
                if coeff * solution[var_idx] == bound
                    push!(vars_at_lower_bound_cs, (idx, c))
                end
            elseif MOI.get(m, MOI.ConstraintName(), c) == "ubs"
                var = constraint_object(c).func.terms.keys[1]
                var_idx = var_mapping[var]
                coeff = constraint_object(c).func.terms[var]
                bound = constraint_object(c).set.upper
                if coeff * solution[var_idx] == bound
                    push!(vars_at_upper_bound_cs, (idx, c))
                end
            end
        end
    end
    # @show vars_at_lower_bound_cs, vars_at_upper_bound_cs

    cons_idxs_to_select = vcat(nonbasic_cons_idxs)
    vars_idxs_to_select = [var_mapping[var_name] for var_name in basic_vars]

    # @show conic_form.A, conic_form.b
    A = Matrix(conic_form.A)[cons_idxs_to_select, :]
    A = A[:, vars_idxs_to_select]
    b = conic_form.b[cons_idxs_to_select]

    # ensure that basis is nxn
    @assert size(A)[1] == size(A)[2]
    @assert size(A)[1] == length(b)

    return (
        A,
        b, 
        basic_vars = basic_vars, 
        nonbasic_cons = nonbasic_cons, 
    )
end

function basis_original_space(S, lb, ub, internal_rxn_idxs)
    m = build_fba_model_inequalitiy_form(S, lb, ub, optimizer=HiGHS.Optimizer, set_objective=true)
    _, num_reactions = size(S)
    _, _, solution_x, _, _ = optimize_model(m)
    @show solution_x
    basis = extract_basis(m) 
    @assert rank(basis.A) == num_reactions
 
    # extract basis of auxiliary LP
    m_aux = build_auxiliary_lp(S, internal_rxn_idxs, lb, ub)    
    _, _, solution_G, _, _ = optimize_model(m_aux)
    basis_aux = extract_basis(m_aux) 
    @show solution_G
    @test rank(basis_aux.A) == length(internal_rxn_idxs)

    # merge bases 
    upper_right = zeros(size(basis.A)[1], size(basis_aux.A)[1])
    lower_left = zeros(size(basis_aux.A)[1], size(basis.A)[1])
    A = hcat(vcat(basis.A, lower_left), vcat(upper_right, basis_aux.A))
    b = vcat(basis.b, basis_aux.b)
    @test isapprox(A * vcat(solution_x, solution_G), b)
    return (A=A, b=b, solution_x, solution_G)
end 

"""
feasible region is set P
"""
function build_relaxed_problem(model, S, internal_rxn_idxs; set_silent=false)
    if set_silent
        set_attribute(model, MOI.Silent(), true)
    else 
        set_attribute(model, MOI.Silent(), false)
    end
    S_int = Array(S[:, internal_rxn_idxs])

    x = model[:x]
    G = @variable(model, G[1:length(internal_rxn_idxs)])
    μ = @variable(model, μ[1:size(S)[1]])

    @constraint(model, G' .== μ' * S_int)

    println("model with " * string(length(x)) * " flux variables and " * string(num_variables(model)) * " variables in total")
end

function build_relaxed_problem_inequality_form(model, S, internal_rxn_idxs; set_silent=false)
    if set_silent
        set_attribute(model, MOI.Silent(), true)
    else 
        set_attribute(model, MOI.Silent(), false)
    end
    S_int = Array(S[:, internal_rxn_idxs])

    x = model[:x]
    G = @variable(model, G[1:length(internal_rxn_idxs)])
    mu = @variable(model, mu[1:size(S)[1]])

    @constraint(model, G' .<= mu' * S_int)
    @constraint(model, -G' .<= -mu' * S_int)

    println("model with " * string(length(x)) * " flux variables and " * string(num_variables(model)) * " variables in total")
end

function build_relaxed_problem_nullspace_inequality_form(model, S, internal_rxn_idxs; set_silent=false)
    if set_silent
        set_attribute(model, MOI.Silent(), true)
    else 
        set_attribute(model, MOI.Silent(), false)
    end
    S_int = Array(S[:, internal_rxn_idxs])
    N_int = nullspace(S_int)

    x = model[:x]
    G = @variable(model, G[1:length(internal_rxn_idxs)])

    @constraint(model, N_int' * G .<= 0)
    @constraint(model, - N_int' * G .<= 0)

    println("model with " * string(length(x)) * " flux variables and " * string(num_variables(model)) * " variables in total")
end

function build_relaxed_problem_standard_form(model, S, internal_rxn_idxs; set_silent=false)
    if set_silent
        set_attribute(model, MOI.Silent(), true)
    else 
        set_attribute(model, MOI.Silent(), false)
    end
    S_int = Array(S[:, internal_rxn_idxs])

    @variable(model, G_pos[1:length(internal_rxn_idxs)])    
    @variable(model, G_neg[1:length(internal_rxn_idxs)])
    @variable(model, μ_pos[1:size(S)[1]])
    @variable(model, μ_neg[1:size(S)[1]])

    @constraint(model, (G_pos - G_neg)' .== (μ_pos - μ_neg)' * S_int)
end

"""
LP with Δμ variables, with bound and nullspace constraints 
"""
function build_auxiliary_lp(S, internal_rxn_idxs, lb, ub; set_silent=false, optimizer=HiGHS.Optimizer)
    model = Model(optimizer)
    if set_silent
        set_attribute(model, MOI.Silent(), true)
    else 
        set_attribute(model, MOI.Silent(), false)
    end
    S_int = Array(S[:, internal_rxn_idxs])
    N_int = nullspace(S_int)
    max_flux_bound = maximum(abs.(vcat(lb, ub)))

    G = @variable(model, G[1:length(internal_rxn_idxs)])

    # bounds on Δμ and nullspace constraint
    # @variable(model, G_pos[1:length(internal_rxn_idxs)])    
    # @variable(model, G_neg[1:length(internal_rxn_idxs)])

    @constraint(model, N_int' * G .<= 0)
    @constraint(model, - N_int' * G .<= 0)

    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        @constraint(model, -G[cidx] <= max_flux_bound)
        @constraint(model, G[cidx] <= max_flux_bound)
    end

    # @objective(model, Max, sum(G))  
    return model
end

"""
get indices for internal reactions with greatest and smallest flux
"""
function get_reaction_index(solution_x, solution_G, internal_rxn_idxs; selected_min_reactions=[], selected_max_reactions=[])
    solution_x_internal = solution_x[internal_rxn_idxs]

    # select reactions that violate thermodynamic feasibility
    infeasible_idxs = []
    for (idx, val) in enumerate(solution_G)
        # WARNING: if flux is zero, we do not select reaction
        # C2
        if solution_x_internal[idx] < -1e-6
            if solution_G[idx] < 1 - 1e-6  
                push!(infeasible_idxs, idx)
                if (idx in [i[1] for i in selected_min_reactions])
                    for c in [i for i in selected_min_reactions if i[1]== idx]
                        if isapprox(solution_x_internal[idx], c[2], atol=1e-3)
                            @assert solution_G[idx] >= c[3]
                        end
                        if isapprox(solution_G[idx], c[2], atol=1e-3)
                            @assert solution_G[idx] >= c[3]
                        end
                    end 
                end 
            end
        # C1
        elseif solution_x_internal[idx] > 1e-6
            if solution_G[idx] > -1 + 1e-6
                push!(infeasible_idxs, idx)
                if (idx in [i[1] for i in selected_max_reactions])
                    for c in [i for i in selected_max_reactions if i[1]== idx]
                        if isapprox(solution_x_internal[idx], c[2], atol=1e-3)
                            if solution_G[idx] > c[3]
                                @show solution_G[idx], c[3]
                                @show solution_x_internal[idx], c[2]
                                @infiltrate
                                @assert solution_G[idx] <= c[3]
                            end
                        end
                        if isapprox(solution_G[idx], c[2], atol=1e-3)
                            @assert solution_G[idx] <= c[3]
                        end
                    end 
                end 
            end
        end
    end

    # get index of reaction in solution_x plus value
    infeasible_values = []
    for (idx, val) in enumerate(internal_rxn_idxs)
        if idx in infeasible_idxs
            push!(infeasible_values, (val, solution_x[val]))
        end
    end

    if !isempty(infeasible_values)
        C1_x_idx = infeasible_values[findmax(i[2] for i in infeasible_values)[2]][1]
        C2_x_idx = infeasible_values[findmin(i[2] for i in infeasible_values)[2]][1]
        C1_G_idx = [idx for (idx, val) in enumerate(internal_rxn_idxs) if val==C1_x_idx][1]
        C2_G_idx = [idx for (idx, val) in enumerate(internal_rxn_idxs) if val==C2_x_idx][1]

        @show solution_G[C1_G_idx], solution_x[C1_x_idx]
        @show solution_G[C2_G_idx], solution_x[C2_x_idx]

        if solution_x[C2_x_idx] > 0
            C2_x_idx = NaN 
            C2_G_idx = NaN
        end
        if solution_x[C1_x_idx] < 0
            C1_x_idx = NaN 
            C1_G_idx = NaN
        end 
    else 
        C2_x_idx = NaN 
        C2_G_idx = NaN
        C1_x_idx = NaN 
        C1_G_idx = NaN
    end
    if !isnan(C1_G_idx)
        C1_G_idx = C1_G_idx + length(solution_x)
    end 
    if !isnan(C2_G_idx)
        C2_G_idx = C2_G_idx + length(solution_x)
    end
    return C1_x_idx, C1_G_idx, C2_x_idx, C2_G_idx
end

"""
add intersection cut if solutoin is in C1 or C2
"""
function intersection_cut(m, A, b, solution; C1_x_idx=NaN, C1_G_idx=NaN, C2_x_idx=NaN, C2_G_idx=NaN)
    @assert (isnan(C2_G_idx) || isnan(C2_x_idx)) || (isnan(C1_G_idx) || isnan(C1_x_idx)) # cut is derived on C1 or C2
    # @show C1_x_idx, C1_G_idx, C2_x_idx, C2_G_idx
    @assert isapprox(A\b, solution) # ensure that inv(A) * b = x

    # compute intersection
    println("compute intersection")
    n = size(A)[1]
    λ = ones(n) * Inf
    λ_1 = ones(n) * Inf
    λ_2 = ones(n) * Inf
    F = lu(A) # lu(big.(basis.A))
    intersecting_points = []

    # compute intersections for each ray
    for idx in 1:n
        e = I[1:n, idx]
        r = - (F\e)
        
        # WARNING: computation can take long
        if r != -inv(A)[:, idx]
            for (i,_) in enumerate(r)
                if r[i] != -inv(A)[:, idx][i]
                    @assert isapprox(r[i], -inv(A)[:, idx][i], atol=1e-8)
                end 
            end
        end

        # intersection with C1
        if !isnan(C1_x_idx)
            # @show r, r[C1_x_idx], r[C1_G_idx] 
            # y-axis -1
            if isapprox(solution[C1_x_idx], -1, atol=1e-8)
                λ_temp = 0
            else 
                λ_temp = -1/r[C1_G_idx] * (1 + solution[C1_G_idx])
            end
            if λ_temp >= 0
                λ_2[idx] = λ_temp
            end 
            # @show "lambda 2 ", λ_temp

            # x-axis
            if isapprox(solution[C1_x_idx], 0, atol=1e-8)
                λ_temp = 0
            else 
                λ_temp = - solution[C1_x_idx]/r[C1_x_idx]
            end
            if λ_temp >= 0
                λ_1[idx] = λ_temp
            end
            # @show "lambda 1 ", λ_temp

        # intersection with C2
        elseif !isnan(C2_x_idx)
            # @show r[min_idx_basic_G], r[min_idx_basic_v]

            # y-axis +1
            if isapprox(solution[C2_G_idx], 1, atol=1e-8)
                    λ_temp = 0
            else 
                λ_temp = 1/r[C2_G_idx] * (1 - solution[C2_G_idx])
            end
            if λ_temp >= 0
                λ_2[idx] = λ_temp
            end 
            # @show "lambda 2 ", λ_temp

            # x-axis
            if isapprox(solution[C2_x_idx], 1, atol=1e-8)
                λ_temp = 0
            else 
                λ_temp = - solution[C2_x_idx]/r[C2_x_idx]
            end
            if λ_temp >= 0
                λ_1[idx] = λ_temp
            end
            # @show "lambda 1 ", λ_temp

        else 
            @error "no reaction selected"
        end

        # compute intersection points
        if !isinf(λ_1[idx]) || !isinf(λ_2[idx])
            dist_1 = sqrt(sum((solution - (solution + λ_1[idx] * r)) .^ 2))
            dist_2 = sqrt(sum((solution - (solution + λ_2[idx] * r)) .^ 2))
            # intersection point with minimal distance is selected
            if dist_1 > dist_2 || isinf(λ_1[idx])
                λ[idx] = λ_2[idx]
            else 
                λ[idx] = λ_1[idx]
            end
            # @show λ[idx], λ_1[idx], λ_2[idx]

            p_1 = solution + λ_1[idx] * r
            p_2 = solution + λ_2[idx] * r
            # @show p_1, p_2

            push!(intersecting_points, round.(solution + λ[idx] * r, digits=7))
        end
    end

    @assert !isempty(intersecting_points)
    @show intersecting_points

    if !(isinf.(λ_1) == ones(n) || isinf.(λ_2) == ones(n))
        if !isnan(C1_x_idx)
            v_vals = [i[C1_x_idx] for i in intersecting_points]
            G_vals = [i[C1_G_idx] for i in intersecting_points]
            @assert [isapprox(v_vals[i], 0, atol=1e-8) || isapprox(G_vals[i], -1, atol=1e-8) for i in 1:length(v_vals)] == ones(length(v_vals))
        else 
            v_vals = [i[C2_x_idx] for i in intersecting_points]
            G_vals = [i[C2_G_idx] for i in intersecting_points]
            @assert [isapprox(v_vals[i], 0, atol=1e-8) || isapprox(G_vals[i], 1, atol=1e-8) for i in 1:length(v_vals)] == ones(length(v_vals))
        end 
        @show v_vals, G_vals
    
        # add intersection cut
        weighted_rows = []
        for (idx, row) in enumerate(eachrow(A))
            if !isinf(λ[idx])
                if isapprox(λ[idx], 0, atol=1e-3)
                    @warn "λ is almost zero"
                    @show λ[idx]
                end
                push!(weighted_rows, (row' * all_variables(m) - b[idx]) / λ[idx])
            end
        end
        # if !isnan(C1_x_idx)
        #     @assert [Bool(i) for i in [isapprox(v_vals[i], 0, atol=1e-8) || isapprox(G_vals[i], -1, atol=1e-8) for i in 1:length(v_vals)]] == ones(length(weighted_rows))
        # elseif !isnan(C1_G_idx)
        #     @assert [Bool(i) for i in [isapprox(v_vals[i], 0, atol=1e-8) || isapprox(G_vals[i], 11, atol=1e-8) for i in 1:length(v_vals)]] == ones(length(weighted_rows))
        # end

        # ensure that efficacy of cut is positive
        eff = efficacy(A, b, solution, λ)
        @show eff
        @assert eff >= 1e-3

        println("add cut")
        @constraint(m, sum(weighted_rows) <= -1)

        # ensure that previous solution is cut off
        weighted_rows = []
        for (idx, row) in enumerate(eachrow(A))
            push!(weighted_rows, (row' * solution - b[idx]) / λ[idx])
        end
        @show sum(weighted_rows)
        @assert !(sum(weighted_rows) <= -1)

    else
        if isinf.(λ_2) == ones(n)
            if !isnan(C2_x_idx)
                println("v < 0 cut off")
                @constraint(m, -m[:x][C2_x_idx] <= 0)
            elseif !isnan(C1_x_idx)
                println("v > 0 cut off")
                @constraint(m, m[:x][C1_x_idx] <= 0)
            end
        elseif isinf.(λ_1) == ones(n)
            if !isnan(C2_x_idx)
                println("G < 1 cut off")
                @constraint(m, -m[:G][C2_G_idx] <= -1)
            elseif !isnan(C1_x_idx)
                println("G > -1 cut off")
                @constraint(m, m[:G][C1_G_idx] <= -1)
            end
        end
    end
end 

# measure cut at optimal solution of current relaxed LP
# Euclidean distance of hyperplane to optimal LP solution
function efficacy(A, b, x, λ)
    n = size(A)[1]
    α_0 = sum([1/λ[i] * b[i] for i in 1:n]) - 1
    α = zeros(n) 
    for j in 1:n
        α[j] = sum([1/λ[i] * A[i,j] for i in 1:n])
    end
    return (dot(α, x) - α_0) / sqrt(sum([i^2 for i in α]))
end 

function intersection_cut_relaxed_problem(organism; verbose=true, C2=false, C1=true, max_cuts=10, set_silent=true, compare_to_ll_fba_solution=false, S=[], lb=[], ub=[], internal_rxn_idxs=[], solution_dict=[])
    # load organism if S, lb, ub not specified
    if !isempty(organism)
        molecular_model = load_model("../molecular_models/" * organism * ".json")
        S = stoichiometry(molecular_model)
        m, _ = size(S)
        internal_rxn_idxs = internal_reactions(molecular_model)
        m = make_optimization_model_ineq_form(molecular_model, HiGHS.Optimizer)
        if verbose 
            print_model(molecular_model, organism)
        end
    else 
        @assert !isempty(S) && !isempty(lb) && !isempty(ub) && !isempty(internal_rxn_idxs)
        m = build_fba_model_inequalitiy_form(S, lb, ub, set_objective=true, optimizer=HiGHS.Optimizer)
    end

    build_relaxed_problem_nullspace_inequality_form(m, S, internal_rxn_idxs, set_silent=set_silent)
    
    # optimal solution from ll FBA
    # TODO: use nullspace formulation
    # if compare_to_ll_fba_solution
    #     if isempty(solution_dict)
    #         file_name = "loopless_fba_1800"
    #         dict = JSON.parse(open("../experiments/json/" * organism * "_" * file_name * ".json"))
    #         @assert dict["termination"] == "OPTIMAL"
    #         ll_fba_solution_G = dict["G"]
    #         ll_fba_solution_μ = dict["μ"]
    #         ll_fba_solution_x = dict["x"]
    #     else 
    #         ll_fba_solution_G = solution_dict[:G]
    #         ll_fba_solution_μ = solution_dict[:μ]
    #         ll_fba_solution_x = solution_dict[:x]
    #     end 
    #     solution_dict = Dict(:G => ll_fba_solution_G, :μ => ll_fba_solution_μ, :x => ll_fba_solution_x)
    #     key_vector = vcat(m[:x], m[:G], m[:μ])
    #     value_vector = vcat(solution_dict[:x], solution_dict[:G], solution_dict[:μ])
    #     point = Dict(key_vector .=> value_vector)
    #     # @show vcat(ll_fba_solution_x, ll_fba_solution_G, ll_fba_solution_μ)
    # end   
    
    cuts = 0
    selected_C1_reactions = []
    selected_C2_reactions = []

    # soluton of relaxed problem
    objective_value, _, _, _, termination = optimize_model(m; silent=set_silent)
    @show termination
    @assert termination == MOI.OPTIMAL
    @show objective_value
    solution_x = [value(var) for var in m[:x]]
    solution_G = [value(var) for var in m[:G]]
    @show solution_x, solution_x[internal_rxn_idxs,] solution_G

    if !isempty(organism)
        @warn "not defined"
    else 
        # check if solution is in interior or in convex hull
        @assert is_basic_feasible_solution(solution_x, solution_G,internal_rxn_idxs; S=S, lb=lb, ub=ub)
        @infiltrate
    end

    C1_x_idx, C1_G_idx, C2_x_idx, C2_G_idx = get_reaction_index(solution_x, solution_G, internal_rxn_idxs)
    @show C1_x_idx, C1_G_idx, C2_x_idx, C2_G_idx

    if C1
        if !isnan(C1_x_idx)
            reaction_to_cut_on = true 
        end
    elseif C2
        if !(isnan(C2_x_idx))
            reaction_to_cut_on = true
        end
    end

    while cuts <= max_cuts && reaction_to_cut_on # stop if no constraint violated
        println("")
        @show cuts
        @show C1_x_idx, C1_G_idx, C2_x_idx, C2_G_idx
    
        # TODO: index of concatenated solution vector
        if compare_to_ll_fba_solution
            if !isnan(C1_x_idx)
                @show ll_fba_solution_x[C1_x_idx], ll_fba_solution_G[C1_G_idx]
                @show solution_x[C1_x_idx], solution_G[C1_G_idx]
            end
            if !isnan(C2_x_idx)
                @show ll_fba_solution_x[C2_x_idx], ll_fba_solution_G[C2_G_idx]
                @show solution_x[C2_x_idx], solution_G[C2_G_idx]
            end
            report = primal_feasibility_report(m, point, atol=0.001)
            if !isempty(report)
                writedlm("report.txt", report)
                @assert isempty(report)
            end
        end

        # extract basis 
        # TODO: extract basis of organism
        basis = basis_original_space(S, lb, ub, internal_rxn_idxs)
        @show basis.solution_x, basis.solution_G

        # add cut
        if C2 
            @assert !isnan(C2_x_idx)
            intersection_cut(m, basis.A, basis.b, vcat(basis.solution_x, basis.solution_G), C2_x_idx=C2_x_idx, C2_G_idx=C2_G_idx)
            push!(selected_C2_reactions, (C2_G_idx, basis.solution_x[C2_x_idx], vcat(basis.solution_x, basis.solution_G)[C2_G_idx])) 
        elseif C1 
            @assert !isnan(C1_x_idx)
            intersection_cut(m, basis.A, basis.b, vcat(basis.solution_x, basis.solution_G), C1_x_idx=C1_x_idx, C1_G_idx=C1_G_idx)
            push!(selected_C1_reactions, (C1_G_idx, basis.solution_x[C1_x_idx], vcat(basis.solution_x, basis.solution_G)[C1_G_idx])) 
        else 
            @error "no reaction selected for intersection cut"
        end

        # optimal solution of relaxed LP after adding intersection cut
        println("optimize after adding cut")
        optimize!(m)
        status = termination_status(m)
        @show status
        @assert status == MOI.OPTIMAL

        # ensure that previous solution violates cut
        key_vector = vcat(m[:x], m[:G])
        value_vector = vcat(basis.solution_x, basis.solution_G)
        point = Dict(key_vector .=> value_vector)
        report = primal_feasibility_report(m, point, atol=0.001)
        @assert !isempty(report)

        # check that solution changed after cut
        # solution_after = [value(var) for var in all_variables(m)]
        # @assert round.(solution, digits=7) != round.(solution_after, digits=7)
        objective_value_after = MOI.get(m, MOI.ObjectiveValue())
        # @show objective_value_after, objective_value
        @assert objective_value_after <= objective_value + 1e-8
        # if C1
        #     @assert solution_x[C1_x_idx] != [value(var) for var in m[:x]][C1_x_idx] || vcat(solution_x, solution_G)[C1_G_idx] != [value(var) for var in m[:G]][C1_G_idx]
        #     @show [value(var) for var in m[:x]][C1_x_idx], [value(var) for var in m[:G]][C1_G_idx]
        # elseif C2
        #     @assert solution_x[C2_x_idx] != [value(var) for var in m[:x]][C2_x_idx] || solution_G[C2_G_idx] != [value(var) for var in m[:G]][C2_G_idx]
        # end
        cuts = cuts +1

        solution_x = [value(var) for var in m[:x]]
        solution_G = [value(var) for var in m[:G]]
        objective_value = objective_value_after
        @show solution_x, solution_x[internal_rxn_idxs], solution_G

        # index of concatenated solution vector
        C1_x_idx, C1_G_idx, C2_x_idx, C2_G_idx = get_reaction_index(solution_x, solution_G, internal_rxn_idxs)
        @show C1_x_idx, C1_G_idx, C2_x_idx, C2_G_idx

        if C1
            if isnan(C1_x_idx)
                reaction_to_cut_on = false 
            end
        elseif C2 
            if (isnan(C2_x_idx))
                reaction_to_cut_on = false
            end
        end
        @show reaction_to_cut_on
    end
    return (
        x = value.(m[:x]),
        G = value.(m[:G]),
    )
end