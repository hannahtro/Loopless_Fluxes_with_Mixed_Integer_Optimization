using HiGHS
using Dates 
using Gurobi

include("../src/loopless_fba.jl")

organism = ARGS[1]
@show organism
time_limit = parse(Int64, ARGS[2])
fast = parse(Bool, ARGS[3])
json = parse(Bool, ARGS[4])
yeast = parse(Bool, ARGS[5])
@show time_limit, json, yeast

start_time = time()

type = "loopless_fba"
try 
    loopless_fba_data(organism, time_limit=time_limit, nullspace_formulation=false, json=json, yeast=yeast, optimizer=Gurobi.Optimizer)
catch e 
    println(e)
    file = organism * "_" * type
    open(file * ".txt","a") do io
        println(io, e)
    end
end

end_time = time()
time_taken = end_time - start_time 
@show time_taken
