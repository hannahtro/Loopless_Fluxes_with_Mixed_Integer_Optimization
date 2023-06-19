include("optimization_model.jl")
include("loopless_constraints.jl")


function no_good_cuts(model, internal_rxn_idxs, S)
    x = model[:x]

    # add indicator variables 
    a = @variable(model, a[1:length(internal_rxn_idxs)], Bin)
    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        # add indicator 
        @constraint(model, a[cidx] => {x[ridx] - eps() >= 0})
        @constraint(model, !a[cidx] => {x[ridx] + eps() <= 0})
    end
    # @objective(model, Max, 0)

    _, _, solution, _, _ = optimize_model(model)
    @show solution

    while !thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)
        # @show thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)

        Z = []
        O = []
        for (idx, ridx) in enumerate(internal_rxn_idxs)
            if solution[ridx] > 0 
                push!(Z,idx)
            else 
                push!(O,idx)
            end
        end 
        # @show Z,O

        @constraint(model, sum(a[Z]) + sum([1-a[i] for i in O]) <= length(internal_rxn_idxs)-1)

        _, _, solution, _, _ = optimize_model(model)

        @show solution
        @show thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)
    end
    return solution
end
