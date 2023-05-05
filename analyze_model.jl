using COBREXA, Serialization, COBREXA.Everything

model = deserialize("data/ec_e_coli_core.js")
println("E COLI COBRA MODEL")
println("number of metabolites : ", length(model.metabolites))
println("number of reactions : ", length(model.reactions))
println("number of genes : ", length(model.genes))
# @show model.annotations
# @show model.notes
println("objective : ", model.objective)

println("")
println("sMOMENT data")
m_rids = deserialize("data/metabolic_rids.js")
println("number of metabolic reaction ids : ", length(m_rids))

t_rids = deserialize("data/transport_rids.js")
println("number of transport reaction ids : ", length(t_rids))

println("")
println("GECKO data")
m_gids = deserialize("data/metabolic_gids.js")
println("number of metabolic gene product ids : ", length(m_gids))

t_gids = deserialize("data/transport_gids.js")
println("number of transport gene product ids : ", length(t_gids))

fieldnames(ObjectModel)
m = load_model(ObjectModel, "data/e_coli_core.json")
@show stoichiometry(model)
@show reactions(model)
@show bounds(model)