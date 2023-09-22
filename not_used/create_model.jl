using COBREXA, Serialization

model = ObjectModel()

metabolite_list = [Metabolite(string("m_", num)) for num = 1:4]
add_metabolites!(model, metabolite_list)


r_1 = ReactionBidirectional("r_1", Dict("m_1" => 1.0))
r_2 = ReactionForward("r_2", Dict("m_1" => -2.0, "m_2" => 1.0, "m_3" => 1.0))
r_3 = ReactionBidirectional("r_3", Dict("m_2" => -1.0))
r_4 = ReactionBidirectional("r_4", Dict("m_3" => -1.0))
r_5 = ReactionForward("r_5", Dict("m_4" => -1.0, "m_1" => 1.0, "m_3" => 3.0))
r_6 = ReactionForward("r_6", Dict("m_4" => 1.0))

add_reactions!(model, [r_1, r_2, r_3, r_4, r_5, r_6])

@show metabolites(model)
@show reactions(model)

# save model
serialize("test_model.js", model)
save_model(model, "test_model.json")