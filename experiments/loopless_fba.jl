#using Gurobi, HiGHS
using Dates 

include("../src/loopless_fba.jl")

organism = ARGS[1]
@show organism
time_limit = parse(Int64, ARGS[2])
json = parse(Bool, ARGS[3])
yeast = parse(Bool, ARGS[4])
mean = parse(Float64, ARGS[5])

@show organism, mean, time_limit, json, yeast

start_time = time()

type = "loopless_fba"
try 
    loopless_fba_data(organism, time_limit=time_limit, yeast=yeast, nullspace_formulation=false, json=json, enzyme_data=true, μ_mean=mean, scip_tol=1.0e-8)
catch e 
    println(e)
    file = organism * "_" * type * "_" * string(mean)
    open(file * ".txt","a") do io
        println(io, e)
    end
end

end_time = time()
time_taken = end_time - start_time 
@show time_taken
