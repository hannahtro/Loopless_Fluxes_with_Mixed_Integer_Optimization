using DelimitedFiles
using COBREXA, COBREXA.Everything
using SCIP 

organisms = readdlm("ecModel_names.txt", '\t', String, '\n')

organism_data = []
for organism in organisms
    molecular_model = load_model("../molecular_models/ecModels/Classical/emodel_" * organism * "_classical.mat")
    S = stoichiometry(molecular_model)
    m, num_reactions = size(S)
    push!(organism_data, (organism, m, num_reactions))
end

# num metabolites
sort!(organism_data, by = organism_data -> organism_data[2])
@show organism_data[1:10]
min_m = organism_data[1:5]
# sort!(organism_data, by = organism_data -> organism_data[2], rev=true)
# @show organism_data[1:10]
println("")

# num reactions
sort!(organism_data, by = organism_data -> organism_data[3])
@show organism_data[1:10]
min_r = organism_data[1:5]
# sort!(organism_data, by = organism_data -> organism_data[3], rev=true)
# @show organism_data[1:10]
println("")

sort!(organism_data, by = organism_data -> organism_data[2] + organism_data[3])
@show organism_data[1:10]
min_sum = organism_data[1:5]
println("")

@show unique(vcat(min_m, min_r, min_sum))
smallest_organisms = [i[1] for i in unique(vcat(min_m, min_r, min_sum))]
writedlm("ecModel_small_model_names.txt", smallest_organisms)

