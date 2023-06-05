using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP
using LinearAlgebra
using Boscia, FrankWolfe

#TODO: why do we get loops in max_flow?
"""
build FBA model
"""
function build_model(S_transform, lb_transform, ub_transform; optimizer=SCIP.Optimizer)
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

    @objective(optimization_model, Max, sum(x))

    return optimization_model
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

