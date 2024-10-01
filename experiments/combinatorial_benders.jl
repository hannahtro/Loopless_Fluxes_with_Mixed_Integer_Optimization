include("../src/cuts_decomposition.jl")

@show ENV["GRB_LICENSE_FILE"]
println(ARGS[1])
organism = ARGS[1]
time_limit = parse(Int64, ARGS[2])
fast = parse(Bool, ARGS[3])
json = parse(Bool, ARGS[4])
yeast = parse(Bool, ARGS[5])
mis = parse(Float64, ARGS[6])
max_density = parse(Int64, ARGS[7])
max_cuts = parse(Float64, ARGS[8])
@show time_limit, fast, json, yeast

type = "cb_fast_big_m"
try 
    combinatorial_benders_data(organism, time_limit=time_limit, fast=fast, json=json, yeast=yeast, indicator=true, big_m=false, multiple_mis=mis, distinct_cuts=true)
catch e 
    println(e)
    file = organism * "_" * type
    open(file * ".txt","a") do io
        println(io, e)
    end
end

