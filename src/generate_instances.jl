using Random
using StatsBase: sample

function generate_model(metabolites=10, reactions=12, exchange=2; seed=1)
    Random.seed!(seed)

    @assert reactions >= 5
    @assert metabolites >= 3 # to have loop
    @assert exchange >= 2 # at least one exchange reaction has to uptake/ consume
    S = zeros(metabolites, reactions)

    # set exchange reactions
    exchange_reaction_idxs = sample(1:reactions, exchange; replace=false)
    exchange_metabolite_idxs = sample(1:metabolites, exchange; replace=false)
    for i in 1:exchange
        if i % 2 == 0
            S[exchange_metabolite_idxs[i], exchange_reaction_idxs[i]] = 1
        else 
            S[exchange_metabolite_idxs[i], exchange_reaction_idxs[i]] = -1
        end
    end

    # set internal reactions
    internal_rxn_idxs = [i for i in 1:reactions if !(i in exchange_reaction_idxs)]
    @show internal_rxn_idxs
    for r in 1:reactions
        if (r in internal_rxn_idxs[4:end])
            m_idxs = sample(1:metabolites, rand(2:min(5, metabolites)); replace=false)
            r_coeffs = [rand(-3:3) for m_idx in m_idxs]
            r_coeffs[1] = 1 
            r_coeffs[end] = -1
            for (i, m_idx) in enumerate(m_idxs)
                S[m_idx, r] = r_coeffs[i]
            end
        end 
    end

    # add internal loop
    loop_idxs = sample(1:metabolites, 3; replace=false)
    S[loop_idxs[1], internal_rxn_idxs[1]] = -1
    S[loop_idxs[2], internal_rxn_idxs[1]] = 1
    S[loop_idxs[2], internal_rxn_idxs[2]] = -1
    S[loop_idxs[3], internal_rxn_idxs[2]] = 1
    S[loop_idxs[3], internal_rxn_idxs[3]] = -1
    S[loop_idxs[1], internal_rxn_idxs[3]] = 1

    # fix direction of exchange directions
    lb = -30 * ones(reactions)
    ub = 30 * ones(reactions)
    for r in 1:reactions
        if !(r in internal_rxn_idxs)
            if sum(S[:,r]) == 1
                lb[r] = 0
            elseif sum(S[:,r]) == -1
                ub[r] = 0
            end
        end 
    end
    return S, lb, ub, internal_rxn_idxs
end

