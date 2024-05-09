using DataFrames, JSON, CSV

include("../src/cycle_free_flux.jl")
include("../src/loopless_fba.jl")

function measure_cycle_search(organism; cff_time_limit=600)
    # load model
    molecular_model = load_model("../molecular_models/" * organism * ".json")
    S = stoichiometry(molecular_model)
    lb, ub = bounds(molecular_model)
    internal_rxn_idxs = internal_reactions(molecular_model)   
    println("-----------------------------------")

    start_time = time()
    cycles = cycle_free_fva(S, internal_rxn_idxs, lb, ub, optimizer=HiGHS.Optimizer, time_limit=cff_time_limit)
    # cycles_to_block = [cycle for cycle in cycles if !thermo_feasible_mu(cycle[1], cycle[2], S)]
    end_time = time()

    @show end_time - start_time
    @show length(cycles)
    unbounded_cycles = [cycle[1] for cycle in cycles]
    flux_directions = [cycle[2] for cycle in cycles]

    # filter infeasible cycles
    infeasible = [thermo_feasible_mu(cycle, flux_directions[idx], S) ? 0 : 1 for (idx, cycle) in enumerate(unbounded_cycles)]
    idx_infeasible = findall(x->x==1, infeasible)
    unbounded_cycles = unbounded_cycles[idx_infeasible]
    flux_directions = flux_directions[idx_infeasible]
    @show length(unbounded_cycles)

    df = DataFrame(
        organism=organism, 
        time_limit=cff_time_limit, 
        time=end_time - start_time, 
        cycles=length(unbounded_cycles)
    )

    file_name = joinpath(@__DIR__,"../experiments/csv/measure_cff.csv")
    if !isfile(file_name)
        CSV.write(file_name, df, append=true, writeheader=true)
    else 
        CSV.write(file_name, df, append=true)    
    end
end