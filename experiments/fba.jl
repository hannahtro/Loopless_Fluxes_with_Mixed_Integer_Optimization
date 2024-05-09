include("../src/fba.jl")

@show ENV["GRB_LICENSE_FILE"]
println(ARGS[1])
organism = ARGS[1]
time_limit = parse(Int64, ARGS[2])
fast = parse(Bool, ARGS[3])
json = parse(Bool, ARGS[4])
yeast = parse(Bool, ARGS[5])
@show time_limit, fast, json, yeast

type = "fba"
try 
    get_fba_data(organism, save_lp=false, yeast=yeast, json=true, optimizer=SCIP.Optimizer)
catch e 
    println(e)
    file = organism * "_" * type
    open(file * ".txt","a") do io
        println(io, e)
    end
end
