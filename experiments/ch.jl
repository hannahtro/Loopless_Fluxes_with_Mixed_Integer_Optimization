using Dates 

include("../src/constraint_handler.jl")

organism = ARGS[1]
@show organism
time_limit = parse(Int64, ARGS[2])
json = parse(Bool, ARGS[3])
multiple_mis = parse(Int64, ARGS[4])
@show time_limit, json

start_time = time()

type = "ch"
try 
    constraint_handler_data(organism; time_limit=time_limit, multiple_mis=multiple_mis)
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