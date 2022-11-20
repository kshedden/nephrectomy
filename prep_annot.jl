using TarIterators, LightXML, CodecZlib, PolygonOps, StaticArrays
using Statistics, Serialization, Printf, JSON, LinearAlgebra

# Path to annotations file
pa = "/home/kshedden/data/Markus_Bitzer/Annotations"

include("defs.jl")

function parse_annot(raw_xml::String)

    xd = parse_string(raw_xml)
    xr = root(xd)

    # Map from nephrectomy id to features
    rd = Dict{String,Vector{Any}}()

    # Loop over the annotation layers.  There should be 14 of these.
    for annot in xr["Annotation"]

        # Annotation type, translate if needed
        atr = annot["Attributes"]
        @assert length(atr) == 1
        atr = atr[1]
        atr = atr["Attribute"]
        if length(atr) != 1
            continue # No name for this layer
        end
        atp = attribute(atr[1], "Name")
        if haskey(trans, atp)
            atp = trans[atp]
        end

        # Get all the points in the current annotation region.
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
    rd = create_normal!(rd)
    return rd
end

# Create a "normal" category consisting of glomeruli that are in "All_glomeruli"
# but not in any atypical class.
function create_normal!(rd)

    if !haskey(rd, "All_glomeruli")
        return rd
    end

    # All atypical glom centroids
    atp = Vector{SVector{2,Float64}}()
    for k in keys(rd)
        if k in atypical_glom_types
            for x in rd[k]
                px = mean(x[1, :])
                py = mean(x[2, :])
                push!(atp, SVector{2,Float64}(px, py))
            end
        end
    end

    nrml = []
    for x in rd["All_glomeruli"]
        px = mean(x[1, :])
        py = mean(x[2, :])
        z = SVector{2,Float64}(px, py)
        di = Inf
        for u in atp
            d = norm(u - z)
            di = min(d, di)
        end
        if di > 200
            push!(nrml, x)
        end
    end
    rd["Normal"] = nrml

    return rd
end

function build_annotations(fn)

    annots = Dict{String,Any}()

    open(GzipDecompressorStream, fn) do io
        ti = TarIterator(io, :file)

        for (h, iox) in ti
            p = h.path
            pp = splitext(p)
            @assert length(pp) == 2 && pp[2] == ".xml"
            println("ID=", pp[1])

            x = read(iox, String)
            y = parse_annot(x)
            annots[pp[1]] = y
        end
    end

    return annots
end

function main()

    fn = joinpath(pa, "annotations.tar.gz")
    annots = build_annotations(fn)

    println(@sprintf("%d samples processed\n", length(annots)))
    open(
        GzipCompressorStream,
        "/home/kshedden/data/Markus_Bitzer/Annotations/annotations.ser.gz",
        "w",
    ) do io
        serialize(io, annots)
    end
    open(
        GzipCompressorStream,
        "/home/kshedden/data/Markus_Bitzer/Annotations/annotations.json.gz",
        "w",
    ) do io
        JSON.print(io, annots)
    end

    return annots
end

annots = main()
