using COBREXA, Serialization, COBREXA.Everything

function json_to_js(organism="iML1515")
    # could not extract json.gz
    molecular_model = load_model(ObjectModel, "../data/" * organism * ".json")
    serialize("../data/" * organism * ".js", molecular_model)
end



