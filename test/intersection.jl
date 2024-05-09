using Test 
using SCIP_PaPILO_jll
using Infiltrator 

include("../src/intersection.jl")
include("../src/cuts_decomposition.jl")
include("../src/generate_instances.jl")
include("../src/cobrexa.jl")
include("../src/PaPILO.jl")

function extract_basis_1(m)
    cs = all_constraints(m, include_variable_in_set_constraints=true)
    solution = value.(all_variables(m))
    conic_form = conic_matrix(m)
    @assert conic_form.A * solution <= conic_form.b .+ 1e-6
    @show size(conic_form.A)
    @show rank(conic_form.A), length(solution)
    @assert rank(conic_form.A) == length(solution)
    active_constraint_idxs = []
    for (idx,row) in enumerate(eachrow(conic_form.A))
        if isapprox.(dot(row, solution), conic_form.b[idx], atol=1e-6)
            push!(active_constraint_idxs, idx)
        end
    end
    A = conic_form.A[active_constraint_idxs,:]
    b = conic_form.b[active_constraint_idxs]
    @assert isapprox.(A * solution, b, atol=1e-6) == ones(length(b))

    # # filter redundant constraints 
    # A_full_rank = []
    # for row_A in eachrow(A)
    #     redundant = [1 for row in A_full_rank if row_A==-row]
    #     print(redundant)
    #     if isempty(redundant)
    #         push!(A_full_rank, row_A)
    #     end
    # end

    return (
        A=A,
        b=b
    )
end

function compute_radius(solution, nearest_feasible_point)
    difference = solution .- nearest_feasible_point
    difference = [i for i in difference if !iszero(i)]
    r = minimum(abs.(difference))
    return r
end

@testset "extract basis of model with loop" begin 
    S = [[1,0,0] [-1,1,0] [0,-1,1] [-1,0,1] [0,0,-1]]
    lb = [0,-10,-10,-10,0]
    ub = [10,10,10,10,10]
    internal_rxn_idxs = [2,3,4]
    m, num_reactions = size(S)
    @show m, num_reactions
    m = build_fba_model_inequalitiy_form(S, lb, ub, set_objective=true, optimizer=HiGHS.Optimizer)
    build_relaxed_problem_inequality_form(m, S, internal_rxn_idxs)
    # print(m)

    objective_value, _, solution, _, termination = optimize_model(m)
    @show length(solution), solution
    conic_form = conic_matrix(m)
    @test conic_form.A * solution <= conic_form.b
    @show size(conic_form.A)

    basis = extract_basis(m) 
    basic_var_idxs = [idx for (idx,var) in enumerate(all_variables(m)) if var in basis.basic_vars]
    @test basis.A * solution[basic_var_idxs] == basis.b[basic_var_idxs]
end 

## EXAMPLE 2 - figure
# setup relaxed LP
@testset "example from visualization" begin 
    println("")
    println("EXAMPLE 2 with conic matrix")
    println("--------------------------------------------------------")

    m = Model(HiGHS.Optimizer)
    @variable(m, x)
    @variable(m, y)

    @constraint(m, x <= 3)
    @constraint(m, -x <= 0)
    @constraint(m, y <= 4)
    @constraint(m, -y <= -0.5)
    @constraint(m, c_1, y <= 1.5x + 0.5) # c_1
    @constraint(m, c_2, y <= -0.5x + 4.5) # c_2
    n = length(all_variables(m))

    @objective(m, Max, y)
    # @objective(m, Min, x)

    var_vec = all_variables(m)
    cs = all_constraints(m, include_variable_in_set_constraints=false)

    # optimal solution of relaxed LP
    set_attribute(m, MOI.Silent(), true)
    optimize!(m)
    solution_before = [value(var) for var in all_variables(m)]
    @show solution_before
    objective_value_before = MOI.get(m, MOI.ObjectiveValue())

    basis = extract_basis_1(m)

    # extreme rays are not normalized
    radius = compute_radius(solution_before, [2.0, 3.0]) # TODO: compute nearest feasible point

    # TODO: infinity case
    rays = [- col for col in eachcol(inv(Matrix(basis.A)))]
    λ = radius ./ norm.(rays)

    # add intersection cut
    @show solution_before + λ[1] * rays[1]
    @show solution_before + λ[2] * rays[2]

    weighted_rows = []
    for (idx, row) in enumerate(eachrow(basis.A))
            @show row
        push!(weighted_rows, (dot(row, vcat(x,y)) - basis.b[idx]) / λ[idx])
    end
    @constraint(m, sum(weighted_rows) <= -1)

    # optimal solution of relaxed LP after adding intersection cut
    optimize!(m)
    solution_after_nonbasic_cons = [value(var) for var in all_variables(m)]
    @show solution_after_nonbasic_cons
    @test round.(solution_before, digits=4) != round.(solution_after_nonbasic_cons, digits=4)
    objective_value_after = MOI.get(m, MOI.ObjectiveValue())
    @show objective_value_after, objective_value_before
    @test objective_value_after <= objective_value_before
end

@testset "example in 3D" begin 
    println("")
    println("EXAMPLE in 3D")
    println("--------------------------------------------------------")

    m = Model(HiGHS.Optimizer)
    @variable(m, x)
    @variable(m, y)
    @variable(m, z)

    @constraint(m, -x <= -0.5)
    @constraint(m, -y <= -0.5)
    @constraint(m, -z <= -0.5)
    @constraint(m, x + 2y + 3z <= 10) # c_2
    n = length(all_variables(m))

    @objective(m, Max, y)
    # @objective(m, Min, x)

    var_vec = all_variables(m)
    cs = all_constraints(m, include_variable_in_set_constraints=false)

    # optimal solution of relaxed LP
    set_attribute(m, MOI.Silent(), true)
    optimize!(m)
    solution_before = [value(var) for var in all_variables(m)]
    @show solution_before
    objective_value_before = MOI.get(m, MOI.ObjectiveValue())

    basis = extract_basis_1(m)

    # extreme rays are not normalized
    radius = 0.5

    # TODO: infinity case
    rays = [- col for col in eachcol(inv(Matrix(basis.A)))]
    λ = radius ./ norm.(rays)

    # add intersection cut
    @show solution_before + λ[1] * rays[1]
    @show solution_before + λ[2] * rays[2]
    @show solution_before + λ[3] * rays[3]

    weighted_rows = []
    for (idx, row) in enumerate(eachrow(basis.A))
            @show row
        push!(weighted_rows, (dot(row, vcat(x,y,z)) - basis.b[idx]) / λ[idx])
    end
    @show weighted_rows
    @show sum(weighted_rows)
    @constraint(m, sum(weighted_rows) <= -1)

    # optimal solution of relaxed LP after adding intersection cut
    optimize!(m)
    solution_after_nonbasic_cons = [value(var) for var in all_variables(m)]
    @show solution_after_nonbasic_cons
    @test round.(solution_before, digits=4) != round.(solution_after_nonbasic_cons, digits=4)
    objective_value_after = MOI.get(m, MOI.ObjectiveValue())
    @show objective_value_after, objective_value_before
    @test objective_value_after <= objective_value_before
end

@testset "extract basis of model with loop" begin 
    S = [[1,0,0] [-1,1,0] [0,-1,1] [-1,0,1] [0,0,-1]]
    lb = [0,-10,-10,-10,0]
    ub = [10,10,10,10,10]
    internal_rxn_idxs = [2,3,4]
    m, num_reactions = size(S)
    @show m, num_reactions

    # extract basis inequality form
    m = build_fba_model_inequalitiy_form(S, lb, ub, set_objective=true, optimizer=HiGHS.Optimizer)
    # print(m)
    build_relaxed_problem_inequality_form(m, S, internal_rxn_idxs)
    # print(m)
    objective_value, _, solution, _, termination = optimize_model(m)
    # basis = extract_basis_1(m) fails because rank(A) != n

    # extract basis standard form
    m_standard_form = build_fba_model_standard_form(S, lb, ub, set_objective=true, optimizer=HiGHS.Optimizer)
    build_relaxed_problem_standard_form(m_standard_form, S, internal_rxn_idxs)
    objective_value, _, solution, _, termination = optimize_model(m_standard_form)
    # basis = extract_basis_1(m_standard_form) fails because rank(A) != n

    # extract basis standard form without presolve
    m_standard_form = build_fba_model_standard_form(S, lb, ub, set_objective=true, optimizer=optimizer_with_attributes(HiGHS.Optimizer, "presolve" => "off"))
    build_relaxed_problem_standard_form(m_standard_form, S, internal_rxn_idxs)
    objective_value, _, solution, _, termination = optimize_model(m_standard_form)
    # basis = extract_basis_1(m_standard_form) fails because rank(A) != n

    # write mps file, presolve with PaPILO
    m = build_fba_model_inequalitiy_form(S, lb, ub, set_objective=true, optimizer=HiGHS.Optimizer)
    build_relaxed_problem_inequality_form(m, S, internal_rxn_idxs)
    # print(m)
    # write_to_file(m, "relaxed_problem_intersection_cut.mps")
    # input_instance = "relaxed_problem_intersection_cut.mps"
    # presolved_instance = "presolved_instance.mps"
    # postsolve_file = "postsolved_instance.post"
    # presolve_write_from_file(input_instance, postsolve_file, presolved_instance)
    model = read_from_file("presolved_instance.mps") # x_3 deleted and all G, μ variables
    set_optimizer(model, HiGHS.Optimizer)
    objective_value, _, solution, _, termination = optimize_model(model)

    # extract A and b
    A = Matrix(lp_matrix_data(model).A)
    b = lp_matrix_data(model).b_upper
    rows, cols = size(A)

    cs = all_constraints(model, include_variable_in_set_constraints=false)
    var_vec = all_variables(model)
    var_mapping = Dict([(val, idx) for (idx, val) in enumerate(var_vec)])
    
    # extract basis
    nonbasic_cons = filter(cs) do c
        MOI.get(model, MOI.ConstraintBasisStatus(), c) != MOI.BASIC
    end
    nonbasic_cons_idxs = [idx for (idx, val) in enumerate(cs) if val in nonbasic_cons]
    basic_vars = filter(var_vec) do v
        MOI.get(model, MOI.VariableBasisStatus(), v) == MOI.BASIC
    end
    cons_idxs_to_select = vcat(nonbasic_cons_idxs)
    vars_idxs_to_select = [var_mapping[var_name] for var_name in basic_vars]

    A = A[cons_idxs_to_select, :]
    # A = A[:, vars_idxs_to_select]
    b = b[cons_idxs_to_select]

    nonbasic_lower_vars = filter(var_vec) do v
        MOI.get(model, MOI.VariableBasisStatus(), v) == MOI.NONBASIC_AT_LOWER
    end
    nonbasic_upper_vars = filter(var_vec) do v
        MOI.get(model, MOI.VariableBasisStatus(), v) == MOI.NONBASIC_AT_UPPER
    end

    # add rows for var at bounds
    cs_var_in_set = [cs for cs in all_constraints(model, include_variable_in_set_constraints=true) if !(cs in all_constraints(model, include_variable_in_set_constraints=false))]
    for cs_var in cs_var_in_set
        new_row = zeros(cols)
        idx = var_mapping[constraint_object(cs_var).func]
        if constraint_object(cs_var).func in nonbasic_lower_vars
            new_row[idx] = -1
            lb = constraint_object(cs_var).set.lower
            A = vcat(A, new_row')
            b = vcat(b, -lb)
        end
        if constraint_object(cs_var).func in nonbasic_upper_vars
            new_row[idx] = 1
            ub = constraint_object(cs_var).set.upper        
            A = vcat(A, new_row')
            b = vcat(b, ub)
        end
    end
    @test rank(A) == 4
end

@testset "extract basis of model with loop with auxiliary LP" begin
    S = [[1,0,0] [-1,1,0] [0,-1,1] [-1,0,1] [0,0,-1]]
    lb = [0,-10,-10,-10,0]
    ub = [10,30,30,30,10]
    internal_rxn_idxs = [2,3,4]
    m, num_reactions = size(S)
    @show m, num_reactions

    # relaxed model 
    relaxed_problem = build_fba_model_inequalitiy_form(S, lb, ub, set_objective=false, optimizer=HiGHS.Optimizer)
    @objective(relaxed_problem, Max, relaxed_problem[:x][2])
    build_relaxed_problem_nullspace_inequality_form(relaxed_problem, S, internal_rxn_idxs)
    # print(relaxed_problem)
    objective_value, _, solution_relaxed_problem, _, termination = optimize_model(relaxed_problem)
    @show solution_relaxed_problem
    solution_x_relaxed_problem = value.(relaxed_problem[:x])
    solution_G_relaxed_problem = value.(relaxed_problem[:G])
    @show solution_x_relaxed_problem, solution_G_relaxed_problem
    basis_relaxed_problem = extract_basis(relaxed_problem)
    @show basis_relaxed_problem.basic_vars # we need nxn basis

    # TODO: solution of subproblems should match solution of relaxed model
    # extract basis inequality form
    m = build_fba_model_inequalitiy_form(S, lb, ub, optimizer=HiGHS.Optimizer)
    @objective(m, Max, m[:x][2])
    # set_start_value.(m[:x], solution_x)
    # print(m)
    objective_value, _, solution_x, _, termination = optimize_model(m)
    @show solution_x
    basis = extract_basis(m) 
    # @show basis.A
    # extract_basis_1(m)
    @test rank(basis.A) == num_reactions

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

    # intersection cut on relaxed problem
    C1_x_idx, C1_G_idx, C2_x_idx, C2_G_idx = get_reaction_index(solution_x, solution_G, internal_rxn_idxs)
    @show C1_x_idx, C1_G_idx
    intersection_cut(relaxed_problem, A, b, vcat(solution_x, solution_G), C1_x_idx=C1_x_idx, C1_G_idx=C1_G_idx)
    objective_value, _, solution, _, termination = optimize_model(relaxed_problem)
    @show solution

    # ensure that solution we derive cut on is no longer feasible
    key_vector = vcat(relaxed_problem[:x], relaxed_problem[:G])
    value_vector = vcat(solution_x, solution_G)
    point = Dict(key_vector .=> value_vector)
    report = primal_feasibility_report(relaxed_problem, point, atol=0.001)
    @test !isempty(report)

    basis_relaxed_problem = extract_basis(relaxed_problem)
    @show basis_relaxed_problem.basic_vars # TODO: How to derive basis if some G variables in basis of relaxed problem?

    # C2 intersection cut on relaxed problem
    @show C2_x_idx, C2_G_idx
    intersection_cut(relaxed_problem, A, b, vcat(solution_x, solution_G), C2_x_idx=C2_x_idx, C2_G_idx=C2_G_idx)
    objective_value, _, solution, _, termination = optimize_model(relaxed_problem)
    @show solution
    report = primal_feasibility_report(relaxed_problem, point, atol=0.001)
    @test !isempty(report)
end

@testset "intersection cut of model with loop with auxiliary LP" begin
    S = [[1,0,0] [-1,1,0] [0,-1,1] [-1,0,1] [0,0,-1]]
    lb = [0,-10,-10,-10,0]
    ub = [10,30,30,30,10]
    internal_rxn_idxs = [2,3,4]

    println("C1 cut")
    intersection_cut_relaxed_problem("", verbose=true, C2=false, C1=true, max_cuts=10, set_silent=true, S=S, lb=lb, ub=ub, internal_rxn_idxs=internal_rxn_idxs)

    println("")
    println("C2 cut")
    intersection_cut_relaxed_problem("", verbose=true, C2=true, C1=false, max_cuts=10, set_silent=true, S=S, lb=lb, ub=ub, internal_rxn_idxs=internal_rxn_idxs)
end

# ### EAXMPLE 2 with variable bound at optimal solution
# @testset "variable bound at optimal solution" begin
#     println("")
#     println("EXAMPLE 2 - active variable bound")
#     println("--------------------------------------------------------")

#     m = Model(HiGHS.Optimizer)
#     @variable(m, 0 <= x <= 3)
#     @variable(m, 0.5 <= y <= 3.5) # has to be greater than 3.5, otherwise c_2 not active

#     @constraint(m, c_1, y <= 1.5x + 0.5) # c_1
#     @constraint(m, c_2, y <= -0.5x + 4.5) # c_2
#     n = length(all_variables(m))

#     @objective(m, Max, y)

#     var_vec = all_variables(m)
#     cs = all_constraints(m, include_variable_in_set_constraints=false)

#     # optimal solution of relaxed LP
#     set_attribute(m, MOI.Silent(), true)
#     optimize!(m)
#     solution_before = [value(var) for var in all_variables(m)]
#     @show solution_before

#     # extract basis
#     basis = extract_basis(m, var_vec, cs)

#     @show det(Matrix(basis.A)) # invertible as det != 0
#     @test isapprox(inv(basis.A) * basis.b, solution_before) #Matrix(B)\b

#     # extreme rays are not normalized
#     # TODO: compute intersection
#     λ = [0.9, 0.5]
#     # compute point on ray
#     for (idx, col) in enumerate(eachcol(inv(Matrix(basis.A))))
#         @show col
#         @show solution_before + λ[idx] * (- col)
#     end

#     # add intersection cut
#     weighted_rows = []
#     for (idx, row) in enumerate(eachrow(basis.A))
#         push!(weighted_rows, (row' * vcat(x,y) - basis.b[idx]) / λ[idx])
#     end
#     @constraint(m, sum(weighted_rows) <= -1)

#     # optimal solution of relaxed LP after adding intersection cut
#     optimize!(m)
#     solution_after = [value(var) for var in all_variables(m)]
#     @show solution_after
#     @test solution_after[2] <= solution_before[2]
# end

# @testset "computing intersection with axes" begin
#     println("")
#     println("axes intersection")
#     println("--------------------------------------------------------")

#     m = Model(HiGHS.Optimizer)
#     @variable(m, -3 <= x <= 3)
#     @variable(m, -4 <= y <= 4)

#     @constraint(m, c_1, y <= 0.268x + 2.62) # c_1
#     @constraint(m, c_2, y <= -0.5x + 4.5) # c_2
#     n = length(all_variables(m))

#     @objective(m, Max, y)

#     var_vec = all_variables(m)
#     cs = all_constraints(m, include_variable_in_set_constraints=false)

#     # optimal solution of relaxed LP
#     set_attribute(m, MOI.Silent(), true)
#     optimize!(m)
#     solution_before = [value(var) for var in all_variables(m)]
#     @show solution_before
#     objective_value_before = MOI.get(m, MOI.ObjectiveValue())

#     # extract basis
#     basis = extract_basis(m, var_vec, cs) 

#     @show det(Matrix(basis.A)) # invertible as det != 0
#     @test isapprox(inv(Matrix(basis.A)) * basis.b, solution_before) #Matrix(B)\b

#     # extreme rays are not normalized
#     # compute intersection
#     # S-free set is all points above x_axis and right of y-axis
#     λ = ones(n) * Inf
#     for (idx, col) in enumerate(eachcol(inv(Matrix(basis.A))))
#         A = hcat(-col, [0, 1]) # y-axis
#         b = - solution_before
#         # @show A\b # [lambda var]
#         # @show (A\b)[1]
#         if (A\b)[1] > 0
#             λ[idx] = (A\b)[1]
#         end 

#         A = hcat(-col, [1, 0]) # x-axis
#         b = - solution_before
#         if (A\b)[1] > 0
#             λ[idx] = min((A\b)[1], λ[idx])
#         end
#     end 
#     intersecting_points = []
#     for (idx, col) in enumerate(eachcol(inv(Matrix(basis.A))))
#         push!(intersecting_points, solution_before + λ[idx] * (- col))
#     end
#     @show intersecting_points
#     @test intersecting_points[1][2] == 0 && intersecting_points[2][1] == 0

#     # add intersection cut
#     weighted_rows = []
#     for (idx, row) in enumerate(eachrow(basis.A))
#         push!(weighted_rows, (row' * vcat(x,y) - basis.b[idx]) / λ[idx])
#     end
#     @constraint(m, sum(weighted_rows) <= -1)

#     # optimal solution of relaxed LP after adding intersection cut
#     optimize!(m)
#     solution_after = [value(var) for var in all_variables(m)]
#     @show solution_after
#     @test round.(solution_before, digits=4) != round.(solution_after, digits=4)
#     @test solution_after == solution_after
#     objective_value_after = MOI.get(m, MOI.ObjectiveValue())
#     @show objective_value_after, objective_value_before
#     @test objective_value_after <= objective_value_before
# end

# @testset "simplified thermodynamic constraint" begin
#     println("")
#     println("thermodynamic constraint")
#     println("--------------------------------------------------------")

#     m = Model(HiGHS.Optimizer)
#     @variable(m, -3 <= x <= 3)
#     @variable(m, -4 <= y <= 4)

#     @constraint(m, c_1, y <= 0.268x + 2.62) # c_1
#     @constraint(m, c_2, y <= -0.5x + 4.5) # c_2
#     n = length(all_variables(m))

#     @objective(m, Max, y)

#     var_vec = all_variables(m)
#     cs = all_constraints(m, include_variable_in_set_constraints=false)

#     # optimal solution of relaxed LP
#     set_attribute(m, MOI.Silent(), true)
#     optimize!(m)
#     solution_before = [value(var) for var in all_variables(m)]
#     @show solution_before
#     objective_value_before = MOI.get(m, MOI.ObjectiveValue())

#     # extract basis
#     basis = extract_basis(m, var_vec, cs) 

#     @show det(Matrix(basis.A)) # invertible as det != 0
#     @test isapprox(inv(Matrix(basis.A)) * basis.b, solution_before) #Matrix(B)\b

#     # compute intersection
#     λ = ones(n) * Inf
#     for (idx, col) in enumerate(eachcol(inv(Matrix(basis.A))))
#         A = hcat(-col, [0, -1]) # y-axis -1
#         b = - solution_before + [-1, 0]
#         if (A\b)[1] > 0
#             λ[idx] = (A\b)[1]
#         end 

#         A = hcat(-col, [-1, 0]) # x-axis
#         b = - solution_before
#         if (A\b)[1] > 0
#             λ[idx] = min((A\b)[1], λ[idx])
#         end
#     end 
#     intersecting_points = []
#     for (idx, col) in enumerate(eachcol(inv(Matrix(basis.A))))
#         push!(intersecting_points, solution_before + λ[idx] * (- col))
#     end
#     @show intersecting_points
#     @test isapprox(intersecting_points[1][2], 0) && isapprox(intersecting_points[2][1], -1)

#     # add intersection cut
#     weighted_rows = []
#     for (idx, row) in enumerate(eachrow(basis.A))
#         push!(weighted_rows, (row' * vcat(x,y) - basis.b[idx]) / λ[idx])
#     end
#     @constraint(m, sum(weighted_rows) <= -1)

#     # optimal solution of relaxed LP after adding intersection cut
#     optimize!(m)
#     solution_after = [value(var) for var in all_variables(m)]
#     @show solution_after
#     @test round.(solution_before, digits=4) != round.(solution_after, digits=4)
#     @test solution_after == solution_after
#     objective_value_after = MOI.get(m, MOI.ObjectiveValue())
#     @show objective_value_after, objective_value_before
#     @test objective_value_after <= objective_value_before
# end

# @testset "artificial molecular network" begin 
#     println("")
#     println("simple graph with infeasible cycle")
#     println("--------------------------------------------------------")

#     S = [[1,0,0] [-1,1,0] [0,-1,1] [-1,0,1] [0,0,-1]]
#     lb = [0,-30,-30,-30,0]
#     ub = [10,30,30,30,10]
#     internal_rxn_idxs = [2,3,4]
#     m, num_reactions = size(S)
#     optimization_model = build_fba_model(S, lb, ub)
#     x = optimization_model[:x]
#     @objective(optimization_model, Max, x[2]+x[3]+x[4])

#     add_loopless_constraints_mu(optimization_model, S, internal_rxn_idxs)
#     objective_value_ll_fba , _, solution_ll_fba, _, _ = optimize_model(optimization_model)
#     solution_dict = Dict(:G => value.(optimization_model[:G]), :μ => value.(optimization_model[:μ]), :x => value.(optimization_model[:x]))
#     @show objective_value_ll_fba
#     @show solution_dict[:x], solution_dict[:G]
#     # @assert iszero.(solution_dict[:x]) != ones(length(solution_dict[:x]))

#     solution_intersection_cut = intersection_cut_relaxed_problem(""; min=false, max=true, compare_to_ll_fba_solution=true, S=S, lb=lb, ub=ub, internal_rxn_idxs=internal_rxn_idxs, solution_dict=solution_dict)
#     @show solution_intersection_cut

#     println("")
#     println("artificial molecular network")
#     println("--------------------------------------------------------")

#     S = [[1,0,0,0] [-1,1,0,0] [0,-1,1,0] [-1,0,1,0] [0,0,-1,0] [-1,0,0,1] [0,0,1,-1]]
#     lb = [0,-10,-10,-10,0,0,0]
#     ub = [20,30,30,30,20,10,10]
#     m, num_reactions = size(S)
#     @show m, num_reactions

#     # test thermodynamic feasible fba with nullspace formulation
#     m = build_fba_model(S, lb, ub, set_objective=true, optimizer=HiGHS.Optimizer)
#     internal_rxn_idxs = [2,3,4,6,7]

#     # print(m)
#     @show length(internal_rxn_idxs)

#     build_relaxed_problem(m, S, internal_rxn_idxs)

#     var_vec = all_variables(m)
#     var_mapping = Dict([(val, idx) for (idx, val) in enumerate(var_vec)])
#     cs = all_constraints(m, include_variable_in_set_constraints=false)

#     objective_value, _, solution, _, termination = optimize_model(m)
#     @test termination == MOI.OPTIMAL
#     @show objective_value
#     solution_x = [value(var) for var in m[:x]]
#     solution_x_internal = solution_x[internal_rxn_idxs]
#     solution_G = [value(var) for var in m[:G]]
#     solution = [value(var) for var in var_vec]

#     infeasible_idxs = []
#     for (idx, val) in enumerate(solution_G)
#         # TODO: some violating cases might not be caught due to ambiguous if clause
#         if solution_x_internal[idx] < 1e-10
#             if solution_G[idx] < 1 + 1e-10  
#                 push!(infeasible_idxs, idx)
#             end
#         elseif solution_x_internal[idx] > -1e-10
#             if solution_G[idx] > -1 + 1e-10
#                 push!(infeasible_idxs, idx)
#             end
#         else
#             @error "val not binary"
#         end
#     end
#     @show infeasible_idxs
#     @show findmax(solution_x_internal)
#     _, max_idx = findmax(solution_x_internal)

#     # TODO: correct ?
#     # extract basis
#     basis = extract_basis(m, var_vec, cs) 

#     @test !isapprox(det(Matrix(basis.A)), 0) # invertible as det != 0
#     basic_var_names = [var_name for var_name in union(basis.basic_vars, basis.nonbasic_lower_vars, basis.nonbasic_upper_vars)]
#     basic_var_idxs = [var_mapping[var_name] for var_name in basic_var_names]
#     solution_basic_vars = solution[sort(basic_var_idxs)]
#     @test isapprox((inv(Matrix(basis.A)) * basis.b), solution_basic_vars)

#     # TODO: not correct
#     # compute intersection
#     max_idx_basic_G = findfirst(x -> x==m[:G][max_idx], basic_var_names)
#     max_idx_basic_v = findfirst(v -> v==m[:x][max_idx], basic_var_names)
#     n = size(basis.A)[1]
#     λ = ones(n) * Inf
#     for (idx, col) in enumerate(eachcol(inv(Matrix(basis.A))))
#         A = col[max_idx_basic_G] # r^i = -col
#         b = (1  + solution_basic_vars[max_idx_basic_G]) # y-axis -1
#         if (A\b)[1] > 0
#             λ[idx] = (A\b)[1]
#         end 

#         A = col[max_idx_basic_v] # x-axis
#         b = solution_basic_vars[max_idx_basic_v]
#         if (A\b)[1] > 0
#             # @show idx, (A\b)[1], λ[idx], min((A\b)[1], λ[idx])
#             λ[idx] = min((A\b)[1], λ[idx])
#         end
#     end 

#     intersecting_points = []
#     for (idx, col) in enumerate(eachcol(inv(Matrix(basis.A))))
#         if !isinf(λ[idx])
#             @show solution_basic_vars + λ[idx] * (- col)
#             push!(intersecting_points, solution_basic_vars + λ[idx] * (- col))
#         end
#     end
#     @show intersecting_points

#     @show solution_basic_vars[max_idx_basic_G], solution_basic_vars[max_idx_basic_v]
#     v_vals = [i[max_idx_basic_v] for i in intersecting_points]
#     G_vals = [i[max_idx_basic_G] for i in intersecting_points]
#     # @test [v_vals[i] == 0 || G_vals[i] == -1 for i in 1:length(v_vals)] == ones(length(v_vals))

#     # add intersection cut
#     weighted_rows = []
#     for (idx, row) in enumerate(eachrow(basis.A))
#         if !isinf(λ[idx])
#             push!(weighted_rows, (row' * all_variables(m)[basic_var_idxs] - basis.b[idx]) / λ[idx])
#         end
#     end
#     # @show [Bool(i) for i in [v_vals[i] == 0 || G_vals[i] == -1 for i in 1:length(v_vals)]]
#     @assert [Bool(i) for i in [v_vals[i] == 0 || G_vals[i] == -1 for i in 1:length(v_vals)]] == ones(length(weighted_rows))
#     @constraint(m, sum(weighted_rows) <= -1)

#     # test that solution is cut off 
#     # TODO: check correctly
#     weighted_rows = []
#     for (idx, row) in enumerate(eachrow(basis.A))
#         push!(weighted_rows, (row' * solution[basic_var_idxs] - basis.b[idx]) / λ[idx])
#     end
#     # TODO: why does assert fail but solution not cut off?
#     # @test sum(weighted_rows) <= -1

#     # optimal solution of relaxed LP after adding intersection cut
#     optimize!(m)
#     status = termination_status(m)
#     @test status == MOI.OPTIMAL

#     solution_after = [value(var) for var in all_variables(m)]
#     @show solution_after
#     @test round.(solution, digits=4) != round.(solution_after, digits=4)
#     objective_value_after = MOI.get(m, MOI.ObjectiveValue())
#     @show objective_value_after, objective_value
#     @test objective_value_after <= objective_value
#     @test solution_x[max_idx] != [value(var) for var in m[:x]][internal_rxn_idxs][max_idx] || solution_G[max_idx] != [value(var) for var in m[:G]][max_idx]

#     @test objective_value_after >= 80.0 # objective value of ll FBA
# end 

# @testset "artificial molecular network using functions in intersection.jl" begin 
#     println("")
#     println("artificial molecular network")
#     println("--------------------------------------------------------")

#     S = [[1,0,0,0] [-1,1,0,0] [0,-1,1,0] [-1,0,1,0] [0,0,-1,0] [-1,0,0,1] [0,0,1,-1]]
#     lb = [0,-10,-10,-10,0,0,0]
#     ub = [20,30,30,30,20,10,10]
#     internal_rxn_idxs = [2,3,4,6,7]
#     S_int = Array(S[:, internal_rxn_idxs])
#     num_metabolites, num_reactions = size(S)
#     @show num_metabolites, num_reactions

#     # ll FBA reference 
#     optimization_model = build_fba_model(S, lb, ub, set_objective=true)
#     # print(model)
#     add_loopless_constraints_mu(optimization_model, S, internal_rxn_idxs)
#     objective_value_ll_fba , _, solution_ll_fba, _, _ = optimize_model(optimization_model)
#     solution_ll_fba_G = value.(optimization_model[:G])
#     solution_ll_fba_μ = value.(optimization_model[:μ])
#     solution_ll_fba_x = value.(optimization_model[:x])
#     @show solution_ll_fba_x, solution_ll_fba_G, solution_ll_fba_μ
#     solution_dict = Dict(:G => solution_ll_fba_G, :μ => solution_ll_fba_μ, :x => solution_ll_fba_x)
#     # C1 cuts
#     println("")
#     println("C1 cuts")
#     intersection_cut_relaxed_problem(""; min=false, max=true, compare_to_ll_fba_solution=true, S=S, lb=lb, ub=ub, internal_rxn_idxs=internal_rxn_idxs, solution_dict=solution_dict)
#     # # C2 cuts
#     # println("")
#     # println("C2 cuts")
#     # intersection_cut_relaxed_problem(""; max=false, min=true, compare_to_ll_fba_solution=true, S=S, lb=lb, ub=ub, internal_rxn_idxs=internal_rxn_idxs, solution_dict=solution_dict)

#     # test thermodynamic feasible fba with nullspace formulation
#     m = build_fba_model(S, lb, ub, set_objective=true, optimizer=HiGHS.Optimizer)
#     build_relaxed_problem(m, S, internal_rxn_idxs)

#     # ll FBA reference point
#     key_vector = vcat(m[:x], m[:G], m[:μ])
#     value_vector = vcat(solution_dict[:x], solution_dict[:G], solution_dict[:μ])
#     point_ll_fba = Dict(key_vector .=> value_vector)

#     println("")
#     println("solution before cut")
#     objective_value, _, solution, _, termination = optimize_model(m)
#     @assert termination == MOI.OPTIMAL
#     @show objective_value
#     solution_x = [value(var) for var in m[:x]]
#     solution_G = [value(var) for var in m[:G]]
#     @show solution_x, solution_G

#     # # add C_2 cut    
#     # println("")
#     # println("C2 cut")
#     # max_x_idx, max_G_idx, min_x_idx, min_G_idx = get_reaction_index(m, internal_rxn_idxs)
#     # @show max_x_idx, max_G_idx, min_x_idx, min_G_idx
#     # intersection_cut(m, min_x_idx=min_x_idx, min_G_idx=min_G_idx, ll_fba_solution=value_vector)
#     # intersection_cut(m, max_x_idx=max_x_idx, max_G_idx=max_G_idx)

#     # # optimal solution of relaxed LP after adding intersection cut
#     # optimize!(m)
#     # status = termination_status(m)
#     # @assert status == MOI.OPTIMAL

#     # solution_after = [value(var) for var in all_variables(m)]
#     # @assert round.(solution, digits=4) != round.(solution_after, digits=4)
#     # objective_value_after = MOI.get(m, MOI.ObjectiveValue())
#     # @show objective_value_after, objective_value
#     # @assert objective_value_after <= objective_value
#     # @assert solution_x[min_x_idx] != [value(var) for var in m[:x]][min_x_idx] || solution_G[min_G_idx] != [value(var) for var in m[:G]][min_G_idx]
#     # @show [value(var) for var in m[:x]][min_x_idx], [value(var) for var in m[:G]][min_G_idx]
#     # @test objective_value_after >= objective_value_ll_fba
#     # solution_x = [value(var) for var in m[:x]]
#     # solution_G = [value(var) for var in m[:G]]
#     # @show solution_x, solution_G

#     # # verify the optimal solution of ll FBA is valid after adding cuts
#     # report = primal_feasibility_report(m, point_ll_fba, atol=0.001)
#     # @test isempty(report)

#     # add C_1 cut
#     println("")
#     println("C1 cut")
#     max_x_idx, max_G_idx, min_x_idx, min_G_idx = get_reaction_index(m, internal_rxn_idxs)
#     @show max_x_idx, max_G_idx, min_x_idx, min_G_idx
#     intersection_cut(m, max_x_idx=max_x_idx, max_G_idx=max_G_idx, ll_fba_solution=value_vector)
#     # optimal solution of relaxed LP after adding intersection cut
#     optimize!(m)
#     status = termination_status(m)
#     @assert status == MOI.OPTIMAL

#     solution_after = [value(var) for var in all_variables(m)]
#     @assert round.(solution, digits=4) != round.(solution_after, digits=4)
#     objective_value_after = MOI.get(m, MOI.ObjectiveValue())
#     @show objective_value_after, objective_value
#     @assert objective_value_after <= objective_value
#     @assert solution_x[max_x_idx] != [value(var) for var in m[:x]][max_x_idx] || solution_G[max_G_idx] != [value(var) for var in m[:G]][max_G_idx]
#     @show [value(var) for var in m[:x]][max_x_idx], [value(var) for var in m[:G]][max_G_idx]
#     @test objective_value_after >= objective_value_ll_fba
#     solution_x = [value(var) for var in m[:x]]
#     solution_G = [value(var) for var in m[:G]]

#     # verify the optimal solution of ll FBA is valid after adding cuts
#     report = primal_feasibility_report(m, point_ll_fba, atol=0.001)
#     @test isempty(report)

#     # add C_1 cut
#     println("")
#     println("C1 cut")
#     max_x_idx, max_G_idx, min_x_idx, min_G_idx = get_reaction_index(m, internal_rxn_idxs)
#     @show max_x_idx, max_G_idx, min_x_idx, min_G_idx
#     intersection_cut(m, max_x_idx=max_x_idx, max_G_idx=max_G_idx, ll_fba_solution=value_vector)
#     # optimal solution of relaxed LP after adding intersection cut
#     optimize!(m)
#     status = termination_status(m)
#     @assert status == MOI.OPTIMAL

#     solution_after = [value(var) for var in all_variables(m)]
#     @assert round.(solution, digits=4) != round.(solution_after, digits=4)
#     objective_value_after = MOI.get(m, MOI.ObjectiveValue())
#     @show objective_value_after, objective_value
#     @assert objective_value_after <= objective_value
#     @assert solution_x[max_x_idx] != [value(var) for var in m[:x]][max_x_idx] || solution_G[max_G_idx] != [value(var) for var in m[:G]][max_G_idx]
#     @show [value(var) for var in m[:x]][max_x_idx], [value(var) for var in m[:G]][max_G_idx]
#     @test objective_value_after >= objective_value_ll_fba
#     solution_x = [value(var) for var in m[:x]]
#     solution_G = [value(var) for var in m[:G]]
#     solution_μ = [value(var) for var in m[:μ]]

#     max_x_idx, max_G_idx, min_x_idx, min_G_idx = get_reaction_index(m, internal_rxn_idxs)
#     @show max_x_idx, max_G_idx, min_x_idx, min_G_idx
#     # @test isempty(filter(!isnan, [max_x_idx, max_G_idx, min_x_idx, min_G_idx])) # no more reaction to branch on

#     @test (solution_μ' * S_int == solution_G')

#     non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution_x) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
#     non_zero_flux_directions = [solution_x[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
#     thermo_feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
#     # @test thermo_feasible

#     # verify the optimal solution of ll FBA is valid after adding cuts
#     # report = primal_feasibility_report(m, point_ll_fba, atol=0.001)
#     # @test isempty(report)
# end

# @testset "another simple model" begin 
#     S = [[1,0,0,0,0,0] [-1,1,0,0,1,0] [0,-1,1,0,0,0] [0,0,-1,0,0,0] [-1,0,0,1,0,0] [0,0,0,-1,1,0] [0,0,1,0,-1,1] [0,0,0,0,0,-1] [0,0,0,-1,0,1]] # example nw with reaction from F to D
#     lb = [0,0,0,0,0,0,0,0,-10]
#     ub = [10,10,10,10,10,10,10,10,10]

#     internal_rxn_idxs = [2,3,5,6,7,9]
#     S_int = Array(S[:, internal_rxn_idxs])
#     num_metabolites, num_reactions = size(S)
#     @show num_metabolites, num_reactions

#     # ll FBA reference 
#     optimization_model = build_fba_model(S, lb, ub, set_objective=true)
#     @objective(optimization_model, Min, optimization_model[:x][9])
#     # print(optimization_model)
#     add_loopless_constraints_mu(optimization_model, S, internal_rxn_idxs)
#     objective_value_ll_fba , _, solution_ll_fba, _, _ = optimize_model(optimization_model)
#     solution_dict = Dict(:G => value.(optimization_model[:G]), :μ => value.(optimization_model[:μ]), :x => value.(optimization_model[:x]))

#     # test thermodynamic feasible fba with nullspace formulation
#     m = build_fba_model(S, lb, ub, set_objective=false, optimizer=HiGHS.Optimizer)
#     @objective(m, Min, m[:x][9])
#     build_relaxed_problem(m, S, internal_rxn_idxs)

#     # ll FBA reference point
#     key_vector = vcat(m[:x], m[:G], m[:μ])
#     value_vector = vcat(solution_dict[:x], solution_dict[:G], solution_dict[:μ])
#     point_ll_fba = Dict(key_vector .=> value_vector)

#     objective_value, _, solution, _, termination = optimize_model(m)
#     @assert termination == MOI.OPTIMAL
#     @show objective_value
#     solution_x = [value(var) for var in m[:x]]
#     solution_G = [value(var) for var in m[:G]]
#     @show solution_x, solution_G

#     # add C_1 cut    
#     println("")
#     println("C1 cut")
#     max_x_idx, max_G_idx, min_x_idx, min_G_idx = get_reaction_index(m, internal_rxn_idxs)
#     @show max_x_idx, max_G_idx, min_x_idx, min_G_idx
#     # intersection_cut(m, min_x_idx=min_x_idx, min_G_idx=min_G_idx)
#     intersection_cut(m, max_x_idx=max_x_idx, max_G_idx=max_G_idx)

#     # optimal solution of relaxed LP after adding intersection cut
#     optimize!(m)
#     status = termination_status(m)
#     @assert status == MOI.OPTIMAL

#     solution_after = [value(var) for var in all_variables(m)]
#     @assert round.(solution, digits=4) != round.(solution_after, digits=4)
#     objective_value_after = MOI.get(m, MOI.ObjectiveValue())
#     @show objective_value_after, objective_value
#     @assert objective_value_after >= objective_value
#     @assert solution_x[max_x_idx] != [value(var) for var in m[:x]][max_x_idx] || solution_G[max_G_idx] != [value(var) for var in m[:G]][max_G_idx]
#     @test objective_value_after >= objective_value_ll_fba
#     solution_x = [value(var) for var in m[:x]]
#     solution_G = [value(var) for var in m[:G]]
# end

# @testset "molecular network" begin
#     println("")
#     println("e coli core")
#     println("--------------------------------------------------------")

#     organism = "e_coli_core"
#     molecular_model = load_model("../molecular_models/" * organism * ".json")
#     S = stoichiometry(molecular_model)
#     m, num_reactions = size(S)
#     lb, ub = bounds(molecular_model)
#     internal_rxn_idxs = internal_reactions(molecular_model)

#     # print_model(molecular_model, organism)
#     @show length(internal_rxn_idxs)

#     m = make_optimization_model(molecular_model, HiGHS.Optimizer)
#     build_relaxed_problem(m, S, internal_rxn_idxs)

#     var_vec = all_variables(m)
#     var_mapping = Dict([(val, idx) for (idx, val) in enumerate(var_vec)])
#     cs = all_constraints(m, include_variable_in_set_constraints=false)

#     objective_value, _, solution, _, termination = optimize_model(m)
#     @test termination == MOI.OPTIMAL
#     @show objective_value
#     solution_x = [value(var) for var in m[:x]]
#     solution_x_internal = solution_x[internal_rxn_idxs]
#     solution_G = [value(var) for var in m[:G]]
#     solution = [value(var) for var in var_vec]

#     infeasible_idxs = []
#     for (idx, val) in enumerate(solution_G)
#         # TODO: some violating cases might not be caught due to ambiguous if clause
#         if solution_x_internal[idx] < 1e-10
#             if solution_G[idx] < 1 + 1e-10  
#                 push!(infeasible_idxs, idx)
#             end
#         elseif solution_x_internal[idx] > -1e-10
#             if solution_G[idx] > -1 + 1e-10
#                 push!(infeasible_idxs, idx)
#             end
#         else
#             @error "val not binary"
#         end
#     end
#     @show infeasible_idxs
#     @show findmax(solution_x_internal)
#     _, max_idx = findmax(solution_x_internal)

#     # TODO: correct ?
#     # extract basis
#     basis = extract_basis(m, var_vec, cs) 

#     @test !isapprox(det(Matrix(basis.A)), 0) # invertible as det != 0
#     basic_var_names = [var_name for var_name in union(basis.basic_vars, basis.nonbasic_lower_vars, basis.nonbasic_upper_vars)]
#     basic_var_idxs = [var_mapping[var_name] for var_name in basic_var_names]
#     solution_basic_vars = solution[sort(basic_var_idxs)]
#     @test isapprox((inv(Matrix(basis.A)) * basis.b), solution_basic_vars)

#     # TODO: not correct
#     # compute intersection
#     max_idx_basic_G = findfirst(x -> x==m[:G][max_idx], basic_var_names)
#     max_idx_basic_v = findfirst(v -> v==m[:x][max_idx], basic_var_names)
#     n = size(basis.A)[1]
#     λ = ones(n) * Inf
#     for (idx, col) in enumerate(eachcol(inv(Matrix(basis.A))))
#         A = col[max_idx_basic_G] # r^i = -col
#         b = (1  + solution_basic_vars[max_idx_basic_G]) # y-axis -1
#         if (A\b)[1] > 0
#             λ[idx] = (A\b)[1]
#         end 

#         A = col[max_idx_basic_v] # x-axis
#         b = solution_basic_vars[max_idx_basic_v]
#         if (A\b)[1] > 0
#             # @show idx, (A\b)[1], λ[idx], min((A\b)[1], λ[idx])
#             λ[idx] = min((A\b)[1], λ[idx])
#         end
#     end

#     intersecting_points = []
#     for (idx, col) in enumerate(eachcol(inv(Matrix(basis.A))))
#         if !isinf(λ[idx])
#             # @show solution_basic_vars + λ[idx] * (- col)
#             push!(intersecting_points, solution_basic_vars + λ[idx] * (- col))
#         end
#     end
#     # @show intersecting_points

#     # @show solution_basic_vars[max_idx_basic_G], solution_basic_vars[max_idx_basic_v]
#     v_vals = [i[max_idx_basic_v] for i in intersecting_points]
#     G_vals = [i[max_idx_basic_G] for i in intersecting_points]
#     # @test [v_vals[i] == 0 || G_vals[i] == -1 for i in 1:length(v_vals)] == ones(length(v_vals))

#     # add intersection cut
#     weighted_rows = []
#     for (idx, row) in enumerate(eachrow(basis.A))
#         if !isinf(λ[idx])
#             push!(weighted_rows, (row' * all_variables(m)[basic_var_idxs] - basis.b[idx]) / λ[idx])
#         end
#     end

#     @assert [Bool(i) for i in [isapprox(v_vals[i], 0, atol=1e-8) || isapprox(G_vals[i], -1, atol=1e-8) for i in 1:length(v_vals)]] == ones(length(weighted_rows))
#     @constraint(m, sum(weighted_rows) <= -1)

#     # test that solution is cut off 
#     # TODO: check correctly
#     weighted_rows = []
#     for (idx, row) in enumerate(eachrow(basis.A))
#         push!(weighted_rows, (row' * solution[basic_var_idxs] - basis.b[idx]) / λ[idx])
#     end
#     # TODO: why does assert fail but solution not cut off?
#     # @test sum(weighted_rows) <= -1

#     # optimal solution of relaxed LP after adding intersection cut
#     optimize!(m)
#     status = termination_status(m)
#     @test status == MOI.OPTIMAL

#     solution_after = [value(var) for var in all_variables(m)]
#     # @show solution_after
#     @test round.(solution, digits=4) != round.(solution_after, digits=4)
#     objective_value_after = MOI.get(m, MOI.ObjectiveValue())
#     @show objective_value_after, objective_value
#     @test objective_value_after <= objective_value
#     @test solution_x[max_idx] != [value(var) for var in m[:x]][internal_rxn_idxs][max_idx] || solution_G[max_idx] != [value(var) for var in m[:G]][max_idx]

#     @test objective_value_after >= 0.8739215 # objective value of ll FBA
# end 

# @testset "molecular network using functions in intersection.jl" begin
#     println("")
#     println("e coli core")
#     println("--------------------------------------------------------")

#     organism = "e_coli_core"
#     molecular_model = load_model("../molecular_models/" * organism * ".json")
#     S = stoichiometry(molecular_model)
#     m, _ = size(S)
#     internal_rxn_idxs = internal_reactions(molecular_model)
#     # print_model(molecular_model, organism)

#     m = make_optimization_model_ineq_form(molecular_model, HiGHS.Optimizer)
#     build_relaxed_problem_inequality_form(m, S, internal_rxn_idxs)

#     objective_value, _, solution, _, termination = optimize_model(m)
#     @assert termination == MOI.OPTIMAL
#     @show objective_value
#     solution_x = [value(var) for var in m[:x]]
#     solution_G = [value(var) for var in m[:G]]

#     # C_1
#     println("")
#     println("C1 cut")
#     max_x_idx, max_G_idx, min_x_idx, min_G_idx = get_reaction_index(m, internal_rxn_idxs)
#     intersection_cut(m, max_x_idx=max_x_idx, max_G_idx=max_G_idx)

#     # optimal solution of relaxed LP after adding intersection cut
#     optimize!(m)
#     status = termination_status(m)
#     @assert status == MOI.OPTIMAL

#     solution_after = [value(var) for var in all_variables(m)]
#     @assert round.(solution, digits=4) != round.(solution_after, digits=4)
#     objective_value_after = MOI.get(m, MOI.ObjectiveValue())
#     @show objective_value_after, objective_value
#     @assert objective_value_after <= objective_value
#     @assert solution_x[max_x_idx] != [value(var) for var in m[:x]][max_x_idx] || solution_G[max_G_idx] != [value(var) for var in m[:G]][max_G_idx]
#     @test objective_value_after >= 0.8739215 # objective value of ll FBA

#     # C_1
#     println("")
#     println("C1 cut")
#     max_x_idx, max_G_idx, min_x_idx, min_G_idx = get_reaction_index(m, internal_rxn_idxs)
#     intersection_cut(m, max_x_idx=max_x_idx, max_G_idx=max_G_idx)

#     # optimal solution of relaxed LP after adding intersection cut
#     optimize!(m)
#     status = termination_status(m)
#     @assert status == MOI.OPTIMAL

#     solution_after = [value(var) for var in all_variables(m)]
#     @assert round.(solution, digits=4) != round.(solution_after, digits=4)
#     objective_value_after = MOI.get(m, MOI.ObjectiveValue())
#     @show objective_value_after, objective_value
#     @assert objective_value_after <= objective_value
#     @assert solution_x[max_x_idx] != [value(var) for var in m[:x]][max_x_idx] || solution_G[max_G_idx] != [value(var) for var in m[:G]][max_G_idx]
#     @test objective_value_after >= 0.8739215 # objective value of ll FBA

#     # C_2 
#     println("")
#     println("C2 cut")
#     # check l.314 in src/intersection.jl
#     max_x_idx, max_G_idx, min_x_idx, min_G_idx = get_reaction_index(m, internal_rxn_idxs)
#     intersection_cut(m, min_x_idx=min_x_idx, min_G_idx=min_G_idx)

#     # optimal solution of relaxed LP after adding intersection cut
#     optimize!(m)
#     status = termination_status(m)
#     @assert status == MOI.OPTIMAL

#     solution_after = [value(var) for var in all_variables(m)]
#     @assert round.(solution, digits=4) != round.(solution_after, digits=4)
#     objective_value_after = MOI.get(m, MOI.ObjectiveValue())
#     @show objective_value_after, objective_value
#     @assert objective_value_after <= objective_value
#     @assert solution_x[max_x_idx] != [value(var) for var in m[:x]][max_x_idx] || solution_G[max_G_idx] != [value(var) for var in m[:G]][max_G_idx]
#     @test objective_value_after >= 0.8739215 # objective value of ll FBA
# end 

# @testset "self-generated instances" begin
#     println("")
#     println("self-generated instances")
#     println("--------------------------------------------------------")

#     println("EXAMPLE 1")
#     S, lb, ub, internal_rxn_idxs = generate_model(4, 5, 2, seed=2)
#     m = build_fba_model_inequalitiy_form(S, lb, ub, set_objective=true, optimizer=HiGHS.Optimizer)
#     # conic_form = conic_matrix(m) # looks as expected
#     build_relaxed_problem_inequality_form(m, S, internal_rxn_idxs)
#     objective_value, _, solution, _, termination = optimize_model(m)
#     # conic_form = conic_matrix(m) # looks as expected
#     var_vec = all_variables(m)
#     cs = all_constraints(m, include_variable_in_set_constraints=false)
#     basis = extract_basis(m, var_vec, cs)
#     max_x_idx, max_G_idx, min_x_idx, min_G_idx = get_reaction_index(m, internal_rxn_idxs, selected_max_reactions=[], selected_min_reactions=[])
#     intersection_cut(m, max_x_idx=max_x_idx, max_G_idx=max_G_idx)

#     # ll FBA reference 
#     optimization_model = build_fba_model(S, lb, ub, set_objective=true)
#     # print(optimization_model)
#     add_loopless_constraints_mu(optimization_model, S, internal_rxn_idxs)
#     objective_value_ll_fba , _, solution_ll_fba, _, _ = optimize_model(optimization_model)
#     solution_dict = Dict(:G => value.(optimization_model[:G]), :μ => value.(optimization_model[:μ]), :x => value.(optimization_model[:x]))
#     @show objective_value_ll_fba
#     @show solution_dict[:x], solution_dict[:G]
#     # @assert iszero.(solution_dict[:x]) != ones(length(solution_dict[:x]))
#     solution = intersection_cut_relaxed_problem(""; min=false, max=true, compare_to_ll_fba_solution=true, S=S, lb=lb, ub=ub, internal_rxn_idxs=internal_rxn_idxs, solution_dict=solution_dict)

#     println("")
#     println("EXAMPLE 2")
#     # S, lb, ub, internal_rxn_idxs = generate_model(7, 10, 2, seed=1)
#     # S, lb, ub, internal_rxn_idxs = generate_model(6, 5, 2, seed=1)
#     S, lb, ub, internal_rxn_idxs = generate_model(7, 8, 2, seed=1)
   
#     m = build_fba_model_inequalitiy_form(S, lb, ub, set_objective=true, optimizer=HiGHS.Optimizer)
#     # conic_form = conic_matrix(m) # looks as expected
#     build_relaxed_problem_inequality_form(m, S, internal_rxn_idxs)
#     objective_value, _, solution, _, termination = optimize_model(m)
#     # conic_form = conic_matrix(m) # looks as expected
#     var_vec = all_variables(m)
#     cs = all_constraints(m, include_variable_in_set_constraints=false)
#     basis = extract_basis(m, var_vec, cs)
#     max_x_idx, max_G_idx, min_x_idx, min_G_idx = get_reaction_index(m, internal_rxn_idxs, selected_max_reactions=[], selected_min_reactions=[])
#     intersection_cut(m, max_x_idx=max_x_idx, max_G_idx=max_G_idx)

#     # ll FBA reference 
#     optimization_model = build_fba_model(S, lb, ub, set_objective=true)
#     # print(optimization_model)
#     add_loopless_constraints_mu(optimization_model, S, internal_rxn_idxs)
#     objective_value_ll_fba , _, solution_ll_fba, _, _ = optimize_model(optimization_model)
#     solution_dict = Dict(:G => value.(optimization_model[:G]), :μ => value.(optimization_model[:μ]), :x => value.(optimization_model[:x]))
#     @show objective_value_ll_fba
#     @show solution_dict[:x], solution_dict[:G]
#     # @assert iszero.(solution_dict[:x]) != ones(length(solution_dict[:x]))
#     solution = intersection_cut_relaxed_problem(""; min=false, max=true, compare_to_ll_fba_solution=true, S=S, lb=lb, ub=ub, internal_rxn_idxs=internal_rxn_idxs, solution_dict=solution_dict)
# end