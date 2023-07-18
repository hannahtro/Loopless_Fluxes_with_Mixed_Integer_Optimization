using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP
using LinearAlgebra
using Boscia, FrankWolfe
# import MathOptInterface
# const MOI = MathOptInterface

#TODO: why do we get loops in max_flow?
"""
build FBA model
"""
function build_fba_model(S_transform, lb_transform, ub_transform; set_objective=false, optimizer=SCIP.Optimizer)
    # make optimization model
    optimization_model = Model(optimizer)
    _, n = size(S_transform)
    # @show size(S_transform)
    # @show size(lb_transform)
    # @show size(ub_transform)

    @variable(optimization_model, x[1:n])
    @constraint(optimization_model, mb, S_transform * x .== 0) # mass balance #TODO set coefficients to -1/1?
    @constraint(optimization_model, lbs, lb_transform .<= x) # lower bounds
    @constraint(optimization_model, ubs, x .<= ub_transform) # upper bounds
    # @show optimization_model

    if set_objective
        @objective(optimization_model, Max, sum(x))
    end 
    
    return optimization_model
end

"""
build FBA model using MOI interface
"""
function build_fba_indicator_model_moi(S_transform, lb_transform, ub_transform, internal_rxn_idxs; set_objective=false, optimizer=SCIP.Optimizer)
    # # make optimization model
    # optimization_model = Model(optimizer)
    # moi_model = direct_model(optimization_model.moi_backend)
    m, n = size(S_transform)
    # # @show size(S_transform)
    # # @show size(lb_transform)
    # # @show size(ub_transform)

    # @variable(moi_model, x[1:n])
    # @constraint(moi_model, mb, S_transform * x .== 0) # mass balance #TODO set coefficients to -1/1?
    # @constraint(moi_model, lbs, lb_transform .<= x) # lower bounds
    # @constraint(moi_model, ubs, x .<= ub_transform) # upper bounds
    # # @show optimization_model
    
    # return moi_model
    moi_model = SCIP.Optimizer()
    x = MOI.add_variables(moi_model, n)
    # MOI.add_constraints(moi_model, S_transform * x .== 0) # mass balance
    S_transform = Float64.(S_transform)
    for i in 1:n
        MOI.add_constraints(moi_model, x[i], MOI.LessThan(Float64(ub_transform[i]))) # upper bounds
        MOI.add_constraints(moi_model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(-1.0, x[i])], 0.0), MOI.LessThan(Float64(-lb_transform[i]))) # lower bounds
    end
    for i in 1:m
        MOI.add_constraints(moi_model, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(S_transform[i,:], x), 0.0), MOI.LessThan(0.0))
        MOI.add_constraints(moi_model, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(-S_transform[i,:], x), 0.0), MOI.LessThan(0.0))
    end

    if set_objective
        MOI.set(moi_model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), 1.0 * x[1])    
    end 

    # add indicator variables 
    a = []
    for i in 1:length(internal_rxn_idxs)
        var = MOI.add_constrained_variable(moi_model, MOI.ZeroOne())
        push!(a,var)
    end 
    a = [i[1] for i in a]
    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        # add indicator 
        f = MOI.VectorAffineFunction(
            [    
                MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, a[cidx])),
                MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(-1.0, x[ridx])),
            ],
            [0.0, 0.0]
        )
        s = MOI.Indicator{MOI.ACTIVATE_ON_ONE}(MOI.LessThan(-eps()))
        MOI.add_constraint(moi_model, f, s)
        # @constraint(moi_model, a[cidx] => {x[ridx] - eps() >= 0})
    end
    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        # add indicator 
        f = MOI.VectorAffineFunction(
            [    
                MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, a[cidx])),
                MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(1.0, x[ridx])),
            ],
            [0.0, 0.0]
        )
        s = MOI.Indicator{MOI.ACTIVATE_ON_ZERO}(MOI.LessThan(-eps()))
        MOI.add_constraint(moi_model, f, s)
        # @constraint(moi_model, !a[cidx] => {x[ridx] + eps() <= 0})
    end
    # print(moi_model)
    return moi_model, a, x
end

"""
print information on COBREXA model
"""
function print_model(model, name="MODEL")
    println("")
    println(name)
    println("----------------------------------")
    println("number of metabolites : ", length(model.metabolites))
    println("number of reactions : ", length(model.reactions))
    println("number of genes : ", length(model.genes))
    # @show model.annotations
    # @show model.notes
    println("objective function: ", model.objective)
    # @show molecular_model.reactions
    # @show molecular_model.metabolites
    println("")
end

"""
optimize model and print process,
returns objective value, solution, time taken and termination status
"""
function optimize_model(model, type="FBA"; time_limit = Inf, print_objective=false, silent=true, mute=true)
    if !mute 
        println("")
        println(type)
        println("----------------------------------")
    end
    if print_objective
        println("objective function : ", objective_function(model))
    end
    if !isinf(time_limit)
        set_time_limit_sec(model, time_limit)
    end
    if silent
        set_attribute(model, MOI.Silent(), true)
    else 
        set_attribute(model, MOI.Silent(), false)
    end
    optimize!(model)
    status = termination_status(model)
    time = solve_time(model)

    # @show solution_summary(model)
    # @show status
    if has_values(model)
        primal_objective_value = MOI.get(model, MOI.ObjectiveValue())
        dual_objective_value = MOI.get(model, MOI.ObjectiveBound())
        if !mute
            println("objective value : ", round(primal_objective_value, digits=2))
            println("")
        end
        solution = [value(var) for var in all_variables(model)]
    else 
        primal_objective_value = NaN
        dual_objective_value = NaN
        solution = NaN
    end

    return primal_objective_value, dual_objective_value, solution, time, status
end

function internal_reactions(molecular_model)
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]
    return internal_rxn_idxs
end