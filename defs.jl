using LightXML

# Path to the annotation data
pa = "/nfs/kshedden/Bitzer/Annotations"

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
