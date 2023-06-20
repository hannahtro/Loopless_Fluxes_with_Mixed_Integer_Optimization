using Dates 

include("optimization_model.jl")
include("loopless_constraints.jl")

"""
add an additional cycle to block until the solution is thermodynamically feasible
"""
function no_good_cuts(model, internal_rxn_idxs, S; time_limit=1800)
    x = model[:x]
    m, num_reactions = size(S)

    # add indicator variables 
    a = @variable(model, a[1:length(internal_rxn_idxs)], Bin)
    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        # add indicator 
        @constraint(model, a[cidx] => {x[ridx] - eps() >= 0})
        @constraint(model, !a[cidx] => {x[ridx] + eps() <= 0})
    end
    # @objective(model, Max, 0)
    # println("here")

    start_time = time()
    # println("here")
    objective_value, dual_bound, solution, _, termination = optimize_model(model)
    # @show solution
    solutions = [round.(solution, digits=5)]
    cuts = []
    # print(model)
    iter = 1
    while !thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S) && time()-start_time < time_limit
        @show iter
        # @show thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)
        solution_a = solution[num_reactions+1:end]
        @assert round.(solution_a) == solution_a

        Z = []
        O = []
        for (idx, ridx) in enumerate(internal_rxn_idxs)
            if solution_a[idx] > 0 
                push!(O,idx)
            else 
                push!(Z,idx)
            end
        end 
        # @show Z,O

        cut = @constraint(model, sum(a[O]) + sum([1-a[i] for i in Z]) <= length(internal_rxn_idxs)-1)
        @assert !(cut in cuts)
        push!(cuts,[cut])
        # print(model)

        objective_value, dual_bound, solution, _, termination = optimize_model(model)
        solution = round.(solution, digits=5)
        solution_a = solution[num_reactions+1:end]

        @assert solutions[end][num_reactions+1:end] != solution_a
        @assert sum(solution_a[O]) + sum([1-solution_a[i] for i in Z]) <= length(internal_rxn_idxs)-1
        @assert !(solution in solutions)
        push!(solutions,solution)
        @show thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)
        iter += 1
    end
    # @show cuts
    # @show solutions
    end_time = time()
    time_taken = end_time - start_time
    @show time_taken

    return objective_value, dual_bound, solution, time_taken, termination, iter
end

function no_good_cuts_data(organism; time_limit=1800, csv=true)
    model = deserialize("../data/" * organism * ".js")
    print_model(model, "organism")

    S = stoichiometry(model)
    lb, ub = bounds(model)
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(model)) if
        !is_boundary(reaction_stoichiometry(model, rid))
    ]

    model = build_model(S, lb, ub)
    objective_value, dual_bound, solution, time, termination, iter = no_good_cuts(model, internal_rxn_idxs, S, time_limit=time_limit)

    thermo_feasible = thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)

    df = DataFrame(
        objective_value=objective_value, 
        dual_bound=dual_bound,
        solution=solution, 
        time=time, 
        termination=termination,
        time_limit=time_limit, 
        thermo_feasible=thermo_feasible,
        iter=iter)

    type = "no_good_cuts"
    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")
    if csv
        CSV.write(file_name, df, append=false, writeheader=true)
    end
end