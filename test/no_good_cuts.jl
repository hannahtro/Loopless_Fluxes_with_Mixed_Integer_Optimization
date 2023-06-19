using Test 

include("../src/cuts_decomposition.jl")

@testset "simple model" begin
    S = [[1,0,0] [-1,1,0] [0,-1,1] [-1,0,1] [0,0,-1]]
    lb = [0,-10,-10,-10,0]
    ub = [10,30,30,30,10]
    m, num_reactions = size(S)
    @show m, num_reactions

    model = build_model(S, lb, ub)
    internal_rxn_idxs = [2,3,4]

    solution = no_good_cuts(model, internal_rxn_idxs, S)

    @test thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)
end