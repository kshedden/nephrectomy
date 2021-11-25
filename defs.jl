using LightXML

# Path to the annotation data
pa = "/home/kshedden/data/Markus_Bitzer/Annotations"

# All annotation files
fi = readdir(pa)
fi = [x for x in fi if endswith(x, ".xml")]

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
    "Description" => "All Glomeruli",
)

# Variables to analyze
avn = [
    :Female,
    :NonWhite,
    :Age,
    :htn_code_pre,
    :dm_code_pre,
    :bmi,
    :htn_sbp,
    :htn,
    :dm_lab_pre,
    :dm,
    :bl_cr_prenx,
    :bl_egfr_prenx,
    :bl_ckd_prenx,
    :kdigo,
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

colors = Dict{String,String}(
    "All Glomeruli" => "grey",
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
)


# Read one annotation file.  Returns a dictionary mapping glomerular phenotypes
# to arrays of bounding boxes.
function read_annot(fn::String)

    xd = parse_file(joinpath(pa, fn))
    xr = root(xd)

    rd = Dict{String,Array{Region}}()

    for annot in xr["Annotation"]
        atp = attribute(annot["Attributes"][1]["Attribute"][1], "Name")
        if haskey(trans, atp)
            atp = trans[atp]
        end
        for region in annot["Regions"]
            for region1 in region["Region"]
                for vertex in region1["Vertices"]
                    vx = []
                    for vertex1 in vertex["Vertex"]
                        x = attribute(vertex1, "X")
                        y = attribute(vertex1, "Y")
                        push!(vx, [parse(Float64, x), parse(Float64, y)])
                    end
                    vx = hcat(vx...)
                    if !haskey(rd, atp)
                        rd[atp] = []
                    end
                    push!(rd[atp], vx)
                end
            end
        end
    end

    free(xd)
    return rd
end

# Collapse all atypical gloms into one classe labeled "Atypical"
function condense(a)

    b = Dict{String,Array{Array{Float64,2},1}}()
    b["All Glomeruli"] = a["All Glomeruli"]
    b["Atypical"] = Float64[]
    for x in ["All Glomeruli", "Tissue", "Capsule", "CMJ"]
        if !haskey(b, x)
            b[x] = Float64[]
        end
        if haskey(a, x)
            push!(b[x], a[x]...)
        end
    end
    for x in ["FSGS", "FGGS", "Ischemic", "Imploding"]
        if haskey(a, x)
            push!(b["Atypical"], a[x]...)
        end
    end

    return b

end
