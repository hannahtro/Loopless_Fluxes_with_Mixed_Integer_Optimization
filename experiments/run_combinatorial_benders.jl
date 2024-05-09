include("../src/cuts_decomposition.jl")
include("../src/constraint_handler.jl")

println(ARGS[1])
time_limit = parse(Int64, ARGS[2])
fast = parse(Bool, ARGS[3])
json = parse(Bool, ARGS[4])
yeast = parse(Bool, ARGS[5])
@show time_limit, fast, json, yeast

type = "cb_fast"
try 
    combinatorial_benders_data(ARGS[1], time_limit=time_limit, fast=fast, json=json, yeast=yeast)
catch e 
    println(e)
    file = organism * "_" * type
    open(file * ".txt","a") do io
        println(io, e)
    end
end

