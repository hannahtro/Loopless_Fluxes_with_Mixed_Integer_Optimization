using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP
using LinearAlgebra
using Boscia, FrankWolfe

"""
"""
function moma(model, x, reference_flux)
    L = - reference_flux
    Q = I(length(x))
    @objective(model, Min, 1/2 * x' * Q * x + L' * x)
    # @show objective_function(model)
    # @show typeof(objective_function(model))
    objective_value_primal, solution, time, status = optimize_model(model, "loopless MOMA"; time_limit=time_limit)

    return objective_value_primal, solution, time, status
end

"""
"""
function moma_boscia(model, x, reference_flux, type="loopless MOMA in Boscia"; time_limit=Inf)
    println("")
    println(type)
    println("----------------------------------")
    
    L = - reference_flux
    Q = I(length(x))

    function f(x)
        length_var = size(Q)[1]
        f_x = 1/2 * x[1:length_var]' * Q * x[1:length_var] + L' * x[1:length_var]
    end
    # @show f(ones(length(x)))
    
    function grad!(storage, x)
        length_var = size(Q)[1]
        storage[1:length_var] = Q * x[1:length_var] + L
        storage[length_var+1:length(x)] .= 0
    end
    
    set_time_limit_sec(model, time_limit)
    set_objective_sense(model, FEASIBILITY_SENSE)
    moi_model = backend(model)
    lmo = FrankWolfe.MathOptLMO(moi_model)
    x, _, result = Boscia.solve(f, grad!, lmo, verbose=true, time_limit=time_limit) 
    objective_value_primal = f(x)
    println("objective value : ", round(objective_value_primal, digits=2))   
    println("")

    return objective_value_primal, x, result[:total_time_in_sec], result[:status]
end