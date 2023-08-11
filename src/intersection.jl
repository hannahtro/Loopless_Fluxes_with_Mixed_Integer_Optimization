using SparseArrays, LinearAlgebra

using JuMP, HiGHS, Clp
using MatrixOptInterface

# m = read_from_file(joinpath(@__DIR__, "..", "csv/models/loopless_fba_relaxed_nullspace_iML1515.lp"))
# m = read_from_file(joinpath(@__DIR__, "..", "csv/models/loopless_fba_relaxed_iSB619.lp"))
m = read_from_file(joinpath(@__DIR__, "..", "lp_models/loopless_fba_relaxed_nullspace_iML1515.lp"))

var_vec = all_variables(m)
cs = all_constraints(m, include_variable_in_set_constraints=false)

var_to_delete = []
for cons in cs
    func = MOI.get(backend(m), MOI.ConstraintFunction(), cons.index)
    @assert func.constant == 0
    set = MOI.get(backend(m), MOI.ConstraintSet(), cons.index)
    if set isa MOI.EqualTo
        @assert MOI.constant(set) == 0
        if length(func.terms) == 1
            push!(var_to_delete, VariableRef(m, func.terms[1].variable))
        end
    elseif length(func.terms) == 1
        @assert set isa MOI.LessThan
        val = MOI.constant(set)
        coeff = func.terms[1].coefficient
        var = func.terms[1].variable
        @assert coeff in (-1, 1)
        if coeff < 0
            set_lower_bound(VariableRef(m, var), -val)
        else
            set_upper_bound(VariableRef(m, var), val)
        end
        JuMP.delete(m, cons)
    end
end

JuMP.delete.(m, unique(var_to_delete))

for v in all_variables(m)
    if occursin("G_", name(v))
        if !has_lower_bound(v)
            set_lower_bound(v, -100)
        end
        if !has_upper_bound(v)
            set_upper_bound(v, 100)
        end
    end
end

set_optimizer(m, HiGHS.Optimizer)
optimize!(m)

cs = all_constraints(m, include_variable_in_set_constraints=false)

basic_cons = filter(cs) do c
    MOI.get(m, MOI.ConstraintBasisStatus(), c) == MOI.BASIC
end

nonbasic_cons = filter(cs) do c
    MOI.get(m, MOI.ConstraintBasisStatus(), c) != MOI.BASIC
end

var_vec = all_variables(m)
basic_vars = filter(var_vec) do v
    MOI.get(m, MOI.VariableBasisStatus(), v) == MOI.BASIC
end
nonbasic_lower_vars = filter(var_vec) do v
    MOI.get(m, MOI.VariableBasisStatus(), v) == MOI.NONBASIC_AT_LOWER
end
nonbasic_upper_vars = filter(var_vec) do v
    MOI.get(m, MOI.VariableBasisStatus(), v) == MOI.NONBASIC_AT_UPPER
end


basis_matrix = zeros(length(var_vec), length(var_vec))
nbasic_cons_names = name.(nonbasic_cons)
basis_rhs = Float64[]
for row_idx in eachindex(nonbasic_cons)
    cons = nonbasic_cons[row_idx]
    func = MOI.get(backend(m), MOI.ConstraintFunction(), cons.index)
    @assert func.constant == 0
    set = MOI.get(backend(m), MOI.ConstraintSet(), cons.index)
    push!(is_gt_cons, set isa MOI.GreaterThan)
    row_multiplier = if set isa MOI.GreaterThan
        -1
    else
        1
    end
    push!(basis_rhs, row_multiplier * MOI.constant(set))
    for term in func.terms
        col_idx = findfirst(v -> v.index == term.variable, var_vec)
        basis_matrix[row_idx, col_idx] = row_multiplier * term.coefficient
    end
end
for idx in eachindex(nonbasic_lower_vars)
    push!(basis_rhs, -lower_bound(nonbasic_lower_vars[idx]))
    col_idx = findfirst(v -> v.index == nonbasic_lower_vars[idx].index , var_vec)
    basis_matrix[idx + length(basic_vars),  col_idx] = -1
    if col_idx == 297
        @show idx + length(basic_vars)
    end
end
for idx in eachindex(nonbasic_upper_vars)
    push!(basis_rhs, upper_bound(nonbasic_upper_vars[idx]))
    col_idx = findfirst(v -> v.index == nonbasic_upper_vars[idx].index, var_vec)
    basis_matrix[idx + length(nonbasic_lower_vars) + length(basic_vars),  col_idx] = 1
end

@assert norm(basis_matrix \ basis_rhs - JuMP.value.(var_vec)) <= 1e-9
