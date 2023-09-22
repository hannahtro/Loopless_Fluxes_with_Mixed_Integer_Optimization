using COBREXA, Serialization
using SCIP, JuMP
using LinearAlgebra
using SparseArrays

function check_symmetry(ridx_1, ridx_2)
    coeffs_1 = sort([vals[idx] for (idx, val) in enumerate(rows) if idx == ridx_1])
    coeffs_2 = sort([vals[idx] for (idx, val) in enumerate(rows) if idx == ridx_1])
    lb_1 = lb[ridx_1]
    lb_2 = lb[ridx_2]
    ub_1 = ub[ridx_1]
    ub_2 = ub[ridx_2]
    return (coeffs_1 == coeffs_2) && (lb_1 == lb_2) & (ub_1 == ub_2)
end

function merge_symmetric_pairs(symmetric_rxns)
    symmetric_sets = []
    while !(isempty(symmetric_rxns))
        temp_set = deepcopy(symmetric_rxns[1])
        # @show temp_set
        to_delete = [1]
        for (idx, val) in enumerate(symmetric_rxns[2:end])
            if !(isempty(intersect(temp_set, val)))
                push!(to_delete, idx+1)
                # @show temp_set, val
                append!(temp_set, val)
            end
        end
        # @show symmetric_rxns
        to_delete = unique(sort(to_delete))
        deleteat!(symmetric_rxns, to_delete)
        # @show symmetric_rxns
        # @show symmetric_sets, temp_set
        temp_set = unique(sort(temp_set))
        push!(symmetric_sets, temp_set)
    end
    return symmetric_sets
end

function get_symmetric_sets(m, num_metabolites)
    symmetric = []

    for i in 5:m
        idxs = [idx for (idx, val) in enumerate(num_metabolites) if val==i]
        if !isempty(idxs)

            # get symmetric pairs
            symmetric_temp = []
            for idx_1 in idxs
                for idx_2 in idxs
                    if idx_1 != idx_2 && !(sort([idx_1, idx_2]) in symmetric_temp)
                        if check_symmetry(idx_1, idx_2)
                            push!(symmetric_temp, sort([idx_1, idx_2]))
                        end
                    end 
                end 
            end

            # merge symmetric pairs
            symmetric_temp_merged = merge_symmetric_pairs(symmetric_temp)
            
            if !isempty(symmetric_temp_merged)
                # @show symmetric_temp
                push!(symmetric, symmetric_temp_merged)
            end
        end
    end
    return symmetric
end

organism = "iAF692"
optimizer = SCIP.Optimizer

molecular_model = load_model("../molecular_models/" * organism * ".json")
S = stoichiometry(molecular_model)
m, num_reactions = size(S)

# @show length([val for val in S.nzval if !(val >= -1 && val <= 1)])
# @show length([val for val in S.nzval if !(val >= -2 && val <= 2)])
# @show length([val for val in S.nzval if !(val >= -3 && val <= 3)])
# @show length([val for val in S.nzval if !(val >= -10 && val <= 10)])
# @show length([val for val in S.nzval if !(val >= -35 && val <= 35)])

lb, ub = bounds(molecular_model)

rows = [] # row
cols = [] # column
vals = []
for (x,y,v) in zip(findnz(S)...)
    push!(rows, x)
    push!(cols, y)
    push!(vals, v)
end

num_metabolites = [length([idx for (idx, val) in enumerate(rows) if val == row]) for row in 1:m] # per reaction

symmetric = get_symmetric_sets(m, num_metabolites)
@show symmetric
