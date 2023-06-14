include("../src/loopless_fba.jl")

@show VERSION

organism = "iAF692"
# iAF692_loopless_fba_blocked_100_1800_10000
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, nullspace_formulation=false)
# iAF692_loopless_fba_blocked_shortest_cycles_100_1800_10000
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, shortest_cycles=true, nullspace_formulation=false)
# iAF692_loopless_fba_blocked_shortest_cycles_50_1800_10000
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, shortest_cycles=true, block_limit=50, nullspace_formulation=false)
# iAF692_loopless_fba_blocked_50_1800_10000
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, block_limit=50, nullspace_formulation=false)

# iAF692_loopless_fba_1800
# loopless_fba_data(organism, time_limit=1800, nullspace_formulation=true)
# # iAF692_loopless_fba_mu_1800
# loopless_fba_data(organism, time_limit=1800, nullspace_formulation=false)

# loopless_indicator_fba_data(organism, time_limit=1800, nullspace_formulation=false)
# loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=50, nullspace_formulation=false)
# loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=100, nullspace_formulation=false)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, same_objective=false, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, same_objective=false, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, same_objective=false, nullspace_formulation=false)

# # iAF692_loopless_fba_blocked_shortest_cycles_100_1800_50
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, vector_formulation=true, shortest_cycles=true, nullspace_formulation=false)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, nullspace_formulation=false)


organism = "iJR904"
# iJR904_loopless_fba_blocked_100_1800_10000
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, nullspace_formulation=false)
# iJR904_loopless_fba_blocked_shortest_cycles_100_1800_10000
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, shortest_cycles=true, nullspace_formulation=false)
# iJR904_loopless_fba_blocked_shortest_cycles_50_1800_10000
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, shortest_cycles=true, block_limit=50, nullspace_formulation=false)
# iJR904_loopless_fba_blocked_shortest_cycles_100_1800_50
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, block_limit=50, nullspace_formulation=false)

# # iJR904_loopless_fba_1800
loopless_fba_data(organism, time_limit=1800, nullspace_formulation=false)
loopless_fba_data(organism, time_limit=1800, nullspace_formulation=true)

# iJR904_loopless_indicator_fba_mu_1800
loopless_indicator_fba_data(organism, time_limit=1800, nullspace_formulation=false)
# iJR904_loopless_indicator_fba_blocked_mu_1800_10_same_objective
loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=50, nullspace_formulation=false)
loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=100, nullspace_formulation=false)

# iJR904_loopless_fba_blocked_1800_50
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, nullspace_formulation=false)
# iJR904_loopless_fba_blocked_100_1800_100
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, same_objective=false, nullspace_formulation=false)
# iJR904_loopless_fba_blocked_100_1800_200
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, same_objective=false, nullspace_formulation=false)
# iJR904_loopless_fba_blocked_100_1800_500
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, same_objective=false, nullspace_formulation=false)

# iJR904_loopless_fba_blocked_shortest_cycles_100_1800_50 TODO: same file twice? 
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, vector_formulation=true, shortest_cycles=true, nullspace_formulation=false)

# iJR904_loopless_fba_blocked_1800_50_same_objective
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, nullspace_formulation=false)
# iJR904_loopless_fba_blocked_100_1800_100_same_objective
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, nullspace_formulation=false)
# iJR904_loopless_fba_blocked_100_1800_200_same_objective
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, nullspace_formulation=false)
# iJR904_loopless_fba_blocked_100_1800_500_same_objective
loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, nullspace_formulation=false)


# organism = "iML1515"
# loopless_fba_data(organism, time_limit=1800, nullspace_formulation=false)
# loopless_fba_data(organism, time_limit=1800, nullspace_formulation=true)

# loopless_indicator_fba_data(organism, time_limit=1800, nullspace_formulation=false)
# loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=50, nullspace_formulation=false)
# loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=100, nullspace_formulation=false)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, same_objective=false, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, same_objective=false, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, same_objective=false, nullspace_formulation=false)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, vector_formulation=true, shortest_cycles=true, nullspace_formulation=false)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, nullspace_formulation=false)