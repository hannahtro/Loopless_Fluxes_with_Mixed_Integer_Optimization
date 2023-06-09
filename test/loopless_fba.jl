include("../src/loopless_fba.jl")

@show VERSION

organism = "iAF692"
# loopless_fba_data(organism, time_limit=600, nullspace_formulation=true)
loopless_fba_data(organism, time_limit=600, nullspace_formulation=false)

organism = "iJR904"
# loopless_fba_data(organism, time_limit=600, nullspace_formulation=true)
loopless_fba_data(organism, time_limit=600, nullspace_formulation=false)

# loopless_indicator_fba_data(organism, time_limit=1800)
# loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=10)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, same_objective=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, same_objective=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, same_objective=false)

# loopless_fba_blocked_data(organism, time_limit=600, ceiling=50, same_objective=false, vector_formulation=true, smallest_cycles=true)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500)


# organism = "iAF692"
# loopless_fba_data(organism, time_limit=1800)

# loopless_indicator_fba_data(organism, time_limit=1800)
# loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=10)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, same_objective=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, same_objective=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, same_objective=false)

# loopless_fba_blocked_data(organism, time_limit=600, ceiling=50, same_objective=false, vector_formulation=true, smallest_cycles=true)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500)


# organism = "iML151"
# loopless_fba_data(organism, time_limit=1800)

# loopless_indicator_fba_data(organism, time_limit=1800)
# loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=10)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, same_objective=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, same_objective=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, same_objective=false)

# loopless_fba_blocked_data(organism, time_limit=600, ceiling=50, same_objective=false, vector_formulation=true, smallest_cycles=true)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500)
