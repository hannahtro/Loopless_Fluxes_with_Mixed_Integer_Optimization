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
# TODO: ensure order of variables after copying
function build_fba_indicator_model_moi(S_transform, lb_transform, ub_transform, internal_rxn_idxs; set_objective=false, optimizer=SCIP.Optimizer, time_limit=Inf, objective_func_vars=Nothing, objective_func_coeffs=Nothing, silent=true)
    # make optimization model
    optimization_model = Model(optimizer)
    model = direct_model(optimization_model.moi_backend)

    m, n = size(S_transform)

    x = @variable(model, x[1:n])
    @constraint(model, mb, S_transform * x .== 0) # mass balance
    @constraint(model, lbs, lb_transform .<= x) # lower bounds
    @constraint(model, ubs, x .<= ub_transform) # upper bounds
    
    if set_objective
        if objective_func_vars != Nothing
            # @show objective_func_coeffs
            # @show objective_func_vars
            @objective(model, Max, objective_func_coeffs' * objective_func_vars)        
        else 
            @objective(model, Max, sum(x))
        end
    end 

    a = build_master_problem(model, internal_rxn_idxs)

    o_inner = SCIP.Optimizer(; presolving_maxrounds=0)
    o = MOI.Bridges.full_bridge_optimizer(o_inner, Float64)
    MOI.copy_to(o, model) # adds complementary v variables
    if silent
        MOI.set(o, MOI.Silent(), true)
    else 
        MOI.set(o, MOI.Silent(), false)
    end
    if !isinf(time_limit)
        MOI.set(o, MOI.TimeLimitSec(), time_limit)
    end

    # print(o)
    # print(o_inner)
    # @show MOI.get(o, MOI.NumberOfVariables())
    # @show MOI.get(o, MOI.ListOfVariableIndices())
    # a = [i.index for i in a]
    # x = [i.index for i in x]
    # append!(x, a)
    # x = MOI.get(o,  MOI.ListOfVariableIndices())
    # @show MOI.get(o, MOI.ListOfConstraintTypesPresent())
    # name = MOI.get(o, MOI.VariableName(), MOI.VariableIndex(1))
    # @show name
    # @show MOI.get(o, MOI.VariableIndex, name)
    # name = MOI.get(o, MOI.VariableName(), MOI.VariableIndex(length(internal_rxn_idxs)+1))
    # @show name
    # @show MOI.get(o, MOI.VariableIndex, name)

    var_names = [MOI.get(o_inner, MOI.VariableName(), MOI.VariableIndex(i)) for i in 1:n+length(internal_rxn_idxs)]
    # @show var_names
    # @show MOI.get(o, MOI.ZeroOne())
    binary_vars = [MOI.VariableIndex(i) for i in 1:length(internal_rxn_idxs)]
    flux_vars = [MOI.VariableIndex(i) for i in length(internal_rxn_idxs)+1:MOI.get(o, MOI.NumberOfVariables())]
    @assert length(flux_vars) == n
    # @show binary_vars
    # @show flux_vars
    return o_inner, binary_vars, flux_vars
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
