
using Test 

include("../src/set_primal.jl")


@testset "set optimal solution as primal in loopless FBA" begin

organism = "iAF692"

objective_value, time, nodes = loopless_fba_data(organism, time_limit=1800)
objective_value_primal, time_primal, nodes_primal = loopless_fba_set_primal(organism, time_limit=1800)

@test isapprox(objective_value_primal,objective_value)
@test time_primal < time
@test nodes_primal < nodes
end