using SparseArrays, LinearAlgebra

using JuMP, HiGHS, Clp
using MatrixOptInterface

# m = read_from_file(joinpath(@__DIR__, "..", "csv/models/loopless_fba_relaxed_nullspace_iML1515.lp"))
m = read_from_file(joinpath(@__DIR__, "..", "csv/models/loopless_fba_relaxed_iSB619.lp"))

# @constraint(m, Gvars .<= 100)
# @constraint(m, -Gvars .<= 100)

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
    end
end

JuMP.delete.(m, unique(var_to_delete))


# we prefer AffineFunction in LessThan
# for g in Gvars
#     set_lower_bound(g, -100)
#     set_upper_bound(g, 100)
# end

## used to convert constraints to variable bounds
# for cons in all_constraints(m, include_variable_in_set_constraints=false)
#     func = MOI.get(backend(m), MOI.ConstraintFunction(), cons.index)
#     set = MOI.get(backend(m), MOI.ConstraintSet(), cons.index)
#     if !(set isa MOI.EqualTo) && length(func.terms) == 1
#         vname = MOI.get(backend(m), MOI.VariableName(), func.terms[1].variable)
#         jump_v = variable_by_name(m, vname)
#         val = MOI.constant(set)
#         if set isa MOI.GreaterThan
#             set_lower_bound(jump_v, val)
#         elseif set isa MOI.LessThan
#             set_upper_bound(jump_v, val)
#         end
#         JuMP.delete(m, cons)
#     end
# end

## Converting everything to {Ax <= b} form.
inner_optimizer = MOI.Utilities.Model{Float64}()
optimizer = MOI.Bridges.Constraint.SplitInterval{Float64}(MOI.Bridges.Constraint.GreaterToLess{Float64}(inner_optimizer))
MOI.copy_to(optimizer, m)

model = Model()
MOI.copy_to(model, optimizer.model.model)

m = model
set_optimizer(m, HiGHS.Optimizer)
MOI.set(model, MOI.RawOptimizerAttribute("presolve"), "off")
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
nonbasic_vars = filter(var_vec) do v
    MOI.get(m, MOI.VariableBasisStatus(), v) != MOI.BASIC
end

# xvars_coupled = [JuMP.variable_by_name() ]
# are our variables sorted and contiguous?
@assert sort(getproperty.(getproperty.(var_vec, :index), :value)) == 1:length(var_vec)

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
        basis_matrix[row_idx, term.variable.value] = row_multiplier * term.coefficient
    end
end
