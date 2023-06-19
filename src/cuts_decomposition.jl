include("optimization_model.jl")
include("loopless_constraints.jl")


function no_good_cuts(model, internal_rxn_idxs, S)
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

    _, _, solution, _, _ = optimize_model(model)
    # @show solution
    solutions = [round.(solution, digits=5)]
    cuts = []

    # print(model)
    iter = 1
    while !thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)
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

        _, _, solution, _, _ = optimize_model(model)
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
    return solution
end
