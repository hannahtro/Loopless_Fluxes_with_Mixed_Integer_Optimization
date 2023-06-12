include("../src/loopless_fba.jl")

@show VERSION

# organism = "iAF692"
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, shortest_cycles=true, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, shortest_cycles=true, block_limit=50, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, block_limit=50, nullspace_formulation=false)

# loopless_fba_data(organism, time_limit=1800, nullspace_formulation=true, nullspace_formulation=false)
# loopless_fba_data(organism, time_limit=1800, nullspace_formulation=false, nullspace_formulation=false)


organism = "iJR904"
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, shortest_cycles=true, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, shortest_cycles=true, block_limit=50, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, block_limit=50, nullspace_formulation=false)


# loopless_fba_data(organism, time_limit=1800, nullspace_formulation=false)

# loopless_indicator_fba_data(organism, time_limit=1800, nullspace_formulation=false)
loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=10, nullspace_formulation=false)

loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, same_objective=false, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, same_objective=false, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, same_objective=false, nullspace_formulation=false)

loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, vector_formulation=true, shortest_cycles=true, nullspace_formulation=false)

loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, nullspace_formulation=false)


organism = "iAF692"
loopless_fba_data(organism, time_limit=1800, nullspace_formulation=false)

loopless_indicator_fba_data(organism, time_limit=1800, nullspace_formulation=false)
loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=10, nullspace_formulation=false)

loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, same_objective=false, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, same_objective=false, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, same_objective=false, nullspace_formulation=false)

loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, vector_formulation=true, shortest_cycles=true, nullspace_formulation=false)

loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, nullspace_formulation=false)


organism = "iML151"
loopless_fba_data(organism, time_limit=1800, nullspace_formulation=false)

loopless_indicator_fba_data(organism, time_limit=1800, nullspace_formulation=false)
loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=10, nullspace_formulation=false)

loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, same_objective=false, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, same_objective=false, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, same_objective=false, nullspace_formulation=false)

loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, vector_formulation=true, shortest_cycles=true, nullspace_formulation=false)

loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, nullspace_formulation=false)
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, nullspace_formulation=false)
