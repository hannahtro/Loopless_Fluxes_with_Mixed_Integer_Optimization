include("optimization_model.jl")
include("loopless_constraints.jl")

S = [[1,0,0] [-1,1,0] [0,-1,1] [-1,0,1] [0,0,-1]]
lb = [0,-10,-10,-10,0]
ub = [10,30,30,30,10]
m, num_reactions = size(S)
@show m, num_reactions

master_problem = build_model(S, lb, ub)

# add indicator variables 

x = master_problem[:x]
internal_rxn_idxs = [2,3,4]
a = @variable(master_problem, a[1:length(internal_rxn_idxs)])

for (cidx, ridx) in enumerate(internal_rxn_idxs)
    # add indicator 
    @constraint(master_problem, a[cidx] => {x[ridx] - eps() >= 0})
    @constraint(master_problem, !a[cidx] => {x[ridx] + eps() <= 0})
end

# @objective(master_problem, Max, 0)

_, _, solution, _, _ = optimize_model(master_problem)

@show solution
@show thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)

Z = []
O = []
for (idx, ridx) in enumerate(internal_rxn_idxs)
    if solution[ridx] > 0 
        push!(Z,idx)
    else 
        push!(O,idx)
    end
end 

@show Z,O

@constraint(master_problem, sum(a[Z]) + sum([1-a[i] for i in O]) <= length(internal_rxn_idxs)-1)

_, _, solution, _, _ = optimize_model(master_problem)

@show solution
@show thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)
