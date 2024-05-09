include("../src/fba.jl")
include("../src/cycle_free_flux.jl")
include("../src/loopless_fba.jl")

println(ARGS[1])
organism = ARGS[1]
seed = parse(Int64, ARGS[2])

@show seed

type = "ll_fba"
try 
    loopless_fba_data(organism; time_limit=1800, type="loopless_fba", nullspace_formulation=false, json=true, scip_seed=seed)
catch e 
    println(e)
    file = organism * "_" * type
    open(file * ".txt","a") do io
        println(io, e)
    end
end

