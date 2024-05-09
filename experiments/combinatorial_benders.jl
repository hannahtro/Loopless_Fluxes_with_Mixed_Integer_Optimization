include("../src/cuts_decomposition.jl")

# @show ENV["GRB_LICENSE_FILE"]
println(ARGS[1])
organism = ARGS[1]
time_limit = parse(Int64, ARGS[2])
fast = parse(Bool, ARGS[3])
json = parse(Bool, ARGS[4])
μ_mean = parse(Float64, ARGS[5])
mis = parse(Int64, ARGS[6])
@show time_limit, fast, json
type = "cb_fast_big_m"
try 
    combinatorial_benders_data(organism, fast=true, json=true, big_m=true, indicator=false, enzyme_data=true, seed=10, mis_solver=HiGHS.Optimizer, scip_tol=1.0e-8, multiple_mis=0, μ_mean=μ_mean)
catch e 
    println(e)
    file = organism * "_" * type * "_" * string(μ_mean)
    open(file * ".txt","a") do io
        println(io, e)
    end
end

