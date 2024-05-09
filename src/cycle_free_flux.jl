using HiGHS, JuMP

include("fba.jl")

function cycle_free_flux(flux, objective_value, S, internal_rxn_idxs, objective_func_vars, objective_func_coeffs; optimizer=HiGHS.Optimizer, constrain_var=NaN, zero_rxns=[])
    m, n = size(S)
    exchange_rxn_idxs = symdiff(1:n, internal_rxn_idxs)
    pos_flux = [idx for (idx, val) in enumerate(flux) if val > 0.0001]
    neg_flux = [idx for (idx, val) in enumerate(flux) if val < -0.0001]

    optimization_model = Model(optimizer)
    model = direct_model(optimization_model.moi_backend)
    x = @variable(model, x[1:n])
    
    # mass balance
    @constraint(model, mb, S * x .== 0) 

    # fix objective value and value of exchange reactions to FBA solution
    @constraint(model, objective_func_coeffs' * x[objective_func_vars] == objective_value)
    @constraint(model, x[exchange_rxn_idxs] == flux[exchange_rxn_idxs])

    # fix flux direction
    @constraint(model, 0 .<= x[pos_flux] .<= flux[pos_flux])
    @constraint(model, 0 .>= x[neg_flux] .>= flux[neg_flux])

    # constrain flux through reaction I
    if !isnan(constrain_var)
        @constraint(model, x[constrain_var]==flux[constrain_var])
    end

    # minimze flux
    @objective(model, Min, sum(x[pos_flux]) - sum(x[neg_flux]))

    objective_val, _, vars, time_fba, termination_fba = optimize_model(optimization_model, print_objective=false)
    return objective_val, vars
end

function cycle_free_fva(S, internal_rxn_idxs, lb, ub; optimizer=HiGHS.Optimizer, time_limit=1800)
    start_time = time()
    current_time = time()
    cycles = []
    for i in 1:length(lb)
        if current_time - start_time < time_limit
            # FBA solution for maximal flux though reaction i
            model = build_fba_model(S, lb, ub; optimizer=optimizer)
            x = model[:x]
            @objective(model, Max, x[i])
            objective_fba, _, solution_fba, _, _ = optimize_model(model, print_objective=false)
            objective_func_vars = [i]
            objective_func_coeffs = [1]
            # println("")
            # @show solution_fba
            # @show i, solution_fba[i]

            # perform CycleFreeFlux
            objective_val_cff, solution_cff = cycle_free_flux(solution_fba, objective_fba, S, internal_rxn_idxs, objective_func_vars, objective_func_coeffs)
            # @show solution_cff

            # @show solution_fba, solution_cff
            # TODO: skip zero fluxes ?
            # if solutions differ, we know that the FBA solution contains a loop
            if solution_fba != solution_cff
                # CycleFreeFlux with fixed reaction i to value in FBA solution 
                # objective_val_cff_2, solution_cff_2 = cycle_free_flux(solution_fba, objective_fba, S, internal_rxn_idxs, objective_func_vars, objective_func_coeffs; constrain_var=i) # TODO: add constraint 
                
                # @show solution_fba, solution_cff, solution_cff_2
                cycle = [idx for (idx,i) in enumerate(round.(solution_fba - solution_cff)) if !isapprox(i, 0, atol=1e-3)]
                dir = [i > 0 ? 1 : 0 for i in solution_fba[cycle]]
                # objective_val_cff_i, solution_cff_i = cycle_free_flux(solution_fba, objective_fba, S, internal_rxn_idxs, objective_func_vars, objective_func_coeffs, constrain_var=i)#, zero_rxns=zero_rxns)
                # @show solution_cff_i
                # cycle = solution_cff_i - solution_cff 
                if !isempty(cycle)
                    push!(cycles, [cycle, dir])
                end
            end
            current_time = time()
        end
    end
    return unique!(cycles)
end
