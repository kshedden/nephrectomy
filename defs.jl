using LightXML

Region = Array{Float64,2}

# Standardize names of phenotypes
trans = Dict{String,String}(
    "GGS" => "FGGS",
    "Imploding Phenotype" => "Imploding",
    "Questionable / Flagged" => "Unclassifiable",
    "Questionble / Flagged" => "Unclassifiable",
    "Normal / No Phenotype" => "Normal",
    "Normal / No Phenoytpe" => "Normal",
    "SGS" => "FSGS",
    "Ischemic- Like" => "Ischemic",
    "Ischemic-like" => "Ischemic",
    "Ischemic-Like" => "Ischemic",
    "Description" => "All_glomeruli",
    "All Glomeruli" => "All_glomeruli",
)

atypical_glom_types = ["FGGS", "FSGS", "BSPC", "Ischemic", "Imploding", "Empty BC"]
glom_types = vcat(["All_glomeruli", "Normal"], atypical_glom_types)

boundary_types = ["Capsule", "CMJ", "Cortex", "Tissue"]

colors = Dict{String,String}(
    "All_glomeruli" => "grey",
    "FGGS" => "blue",
    "FSGS" => "magenta",
    "BSPC" => "green",
    "Ischemic" => "cyan",
    "Normal" => "black",
    "Imploding" => "red",
    "Empty BC" => "lime",
    "Capsule" => "purple",
    "CMJ" => "orange",
    "Tissue" => "yellow",
    "Atypical" => "blue",
    "Cortex" => "violet",
)
