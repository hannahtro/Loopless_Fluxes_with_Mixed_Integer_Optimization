using GLPK
include("../src/disjunctive_programming.jl")

println(ARGS[1])
organism = ARGS[1]
time_limit = parse(Int64, ARGS[2])
gdp_method = ARGS[3]

type = "dp"
try 
    dp_data(organism; gdp_method=gdp_method, big_m_constant=1000, optimizer=HiGHS.Optimizer)
catch e 
    println(e)
    file = organism * "_" * type * "_" * gdp_method 
    open(file * ".txt","a") do io
        println(io, e)
    end
end
