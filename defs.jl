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

# Variables to analyze
avn = [
    :Female,
    :NonWhite,
    :CURRENT_AGE,
    #:htn_code_pre,
    #:dm_code_pre,
    #:bmi,
    #:htn_sbp,
    :Lab_Result,
    :htn,
    #:dm_lab_pre,
    :dm_med,
    :Smoking,
    #:bl_cr_prenx,
    #:bl_egfr_prenx,
    #:bl_ckd_prenx,
    #:kdigo,
    :"Hyper-Bin",
    :"T2D-Bin",
    :"eGFR Min4 Mean",
    :"I/L Ratio",
    :">100%",
    :">50%",
    :NormalGlomPercent,
    :AbnormalPercent,
    :"AKI-Bin",
    :GGS,
    :Imploding,
    :GlomVol,
    :"% CA Glomerular",
    :"Nephron Density (All)",
    :"Nephron Density (Normal)",
    :"GV/Podocyte",
    :PodoDensity,
    :MesIndex,
    :PodoPerGlom,
    :MesVol,
    :PodoVol,
    :FIA,
    :SGS,
    :Ischemic,
    :vGlepp,
]

atypical_glom_types = ["FGGS", "FSGS", "BSPC", "Ischemic", "Imploding"]
glom_types = vcat(["All_glomeruli", "Normal"], atypical_glom_types)

colors = Dict{String,String}(
    "All_glomeruli" => "grey",
    "FGGS" => "blue",
    "FSGS" => "magenta",
    "BSPC" => "green",
    "Ischemic" => "cyan",
    "Normal" => "black",
    "Imploding" => "red",
    "Capsule" => "purple",
    "CMJ" => "orange",
    "Tissue" => "yellow",
    "Atypical" => "blue",
    "Cortex" => "violet",
)
