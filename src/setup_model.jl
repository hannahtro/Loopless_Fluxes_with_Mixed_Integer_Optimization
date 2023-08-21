using COBREXA, Serialization, COBREXA.Everything
using DelimitedFiles 

"""
store json as js file
TODO: COBREXA can read json directly, therefore serialization not needed 
"""
function json_to_js(organism="iML1515")
    # could not extract json.gz
    molecular_model = load_model(ObjectModel, "../molecular_models/" * organism * ".json")
    serialize("../molecular_models/" * organism * ".js", molecular_model)
end


"""
get names of ecModels from https://zenodo.org/record/6438262#.Y3NLqi8w0dU
"""
function scrape_model_names(path = "../molecular_models/ecModels/ecModels/Classical/")
    files = readdir(path; join=true)
    files = replace.(files, r"../molecular_models/ecModels/ecModels/Classical/emodel_" => "")
    files = replace.(files, r"_classical.mat" => "")
    writedlm("../molecular_models/ecModel_names.txt", files)
end
