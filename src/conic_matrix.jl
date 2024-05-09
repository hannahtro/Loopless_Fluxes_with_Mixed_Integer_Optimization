using JuMP
using MathOptInterface

function conic_matrix(model::GenericModel{T}) where {T}
    columns = Dict(var => i for (i, var) in enumerate(all_variables(model)))
    I_bounds, J_bounds, V_bounds = Int[], Int[], T[]
    bound_constraints = ConstraintRef[]
    affine_constraints = ConstraintRef[]
    n = length(columns)
    bound_constraint_number = length(all_constraints(model, VariableRef, MathOptInterface.GreaterThan{Float64})) + length(all_constraints(model, VariableRef, MathOptInterface.LessThan{Float64}))
    c_u = fill(typemax(T), bound_constraint_number) # comes from bound constraints
    r_l, r_u = T[], T[]
    I, J, V = Int[], Int[], T[]
    for (F, S) in list_of_constraint_types(model)
        _fill_conic_form(
            model,
            columns,
            bound_constraints,
            affine_constraints,
            F,
            S,
            c_u,
            r_l,
            r_u,
            I,
            J,
            V,
            I_bounds,
            J_bounds, 
            V_bounds
        )
    end
    rows_to_add = sort(unique(I))[end]
    I_bounds = I_bounds .+ rows_to_add 
    # @show r_u, c_u
    columns = columns
    if isinf.(c_u) == ones(length(c_u))
        b = r_u
    else
        b = vcat(r_u, c_u)
    end 
    A = SparseArrays.sparse(append!(I, I_bounds), append!(J, J_bounds), append!(V, V_bounds))
    bounds = bound_constraints
    constraints = affine_constraints
    if size(A)[1] != length(b)
        @show size(A)[1], length(b)
    end
    @assert size(A)[1] == length(b)
    return (
        columns = columns,
        b = b,
        A = A,
        bounds = bounds,
        constraints = constraints,
    )
end

function _fill_conic_form(
    model::GenericModel{T},
    x::Dict{GenericVariableRef{T},Int},
    bound_constraints::Vector{ConstraintRef},
    ::Vector{ConstraintRef},
    F::Type{GenericVariableRef{T}},
    S::Type,
    c_u::Vector{T},
    ::Vector{T},
    ::Vector{T},
    ::Vector{Int},
    ::Vector{Int},
    ::Vector{T},
    I_bounds::Vector{Int},
    J_bounds::Vector{Int},
    V_bounds::Vector{T},
) where {T}
    for c in all_constraints(model, F, S)
        push!(bound_constraints, c)
        row = length(I_bounds) + 1
        c_obj = constraint_object(c)
        var = c_obj.func
        i = x[c_obj.func]
        set = MOI.Interval(c_obj.set)
        # @show set.lower, set.upper
        # @show c_u[row]
        if S == MOI.GreaterThan{Float64}
            c_u[row] = min(c_u[row], set.lower)
            val = -1
        elseif S == MOI.LessThan{Float64}
            c_u[row] = min(c_u[row], set.upper)
            val = 1
        end
        # @show c_u[row]
        push!(I_bounds, row)
        push!(J_bounds, findall(i->i==c_obj.func, all_variables(model))[1]) # get index
        push!(V_bounds, val)  
    end
    return
end

function _fill_conic_form(
    model::GenericModel{T},
    x::Dict{GenericVariableRef{T},Int},
    ::Vector{ConstraintRef},
    affine_constraints::Vector{ConstraintRef},
    F::Type{<:GenericAffExpr},
    S::Type,
    ::Vector{T},
    r_l::Vector{T},
    r_u::Vector{T},
    I::Vector{Int},
    J::Vector{Int},
    V::Vector{T},
    ::Vector{Int},
    ::Vector{Int},
    ::Vector{T},
) where {T}
    for c in all_constraints(model, F, S)
        push!(affine_constraints, c)
        c_obj = constraint_object(c)
        @assert iszero(c_obj.func.constant)
        row = length(r_l) + 1
        set = MOI.Interval(c_obj.set)
        # @show set.upper
        push!(r_l, set.lower)
        push!(r_u, set.upper)
        for (var, coef) in c_obj.func.terms
            push!(I, row)
            push!(J, x[var])
            push!(V, coef)
        end
    end
    return
end
