using TarIterators, LightXML, GZip, PolygonOps, StaticArrays
using Statistics, Serialization

# Path to annotations file
fn = "/home/kshedden/data/Markus_Bitzer/Annotations/annotations.tar.gz"

include("defs.jl")

function parse_annot(raw_xml::String)

    xd = parse_string(raw_xml)
    xr = root(xd)

    # Map from nephrectomy id to features
    rd = Dict{String,Vector{Any}}()

    # Loop over the annotation layers
    for annot in xr["Annotation"]

        # Annotation type, translate if needed
        atp = attribute(annot["Attributes"][1]["Attribute"][1], "Name")
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
    rd = dedup(rd)
    return rd
end

# Create a "normal" category consisting of glomeruli that are in "All_glomeruli"
# but not in any atypical class.
function dedup(rd)

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
            d = sqrt(sum(abs2, u - z))
            di = d < di ? d : di
        end
        if di > 200
            push!(nrml, x)
        end
    end
    rd["Normal"] = nrml

    return rd
end

function find_components(anno)

    components = Dict{String,Any}()

    # The bounding polygon for each component
    polys = []
    if haskey(anno, "Cortex")
        for (k, ti) in enumerate(anno["Cortex"])
            tix = SVector{2,Float64}[]
            for i = 1:size(ti, 2)
                push!(tix, ti[:, i])
            end
            push!(tix, tix[1]) # Close the loop
            push!(polys, tix)
        end
    end
    m = length(polys)

    # Determine the component to which each glomerulus
    # belongs.
    for g in glom_types
        if !(g in keys(anno))
            continue
        end

        gloms = anno[g]
        cmp = Vector{Union{Int,Nothing}}()

        for (i, glom) in enumerate(gloms)
            # The centroid of the glomerulus
            c = mean(glom, dims = 2)

            ix = [j for (j, pl) in enumerate(polys) if inpolygon(c, pl) != 0]
            if length(ix) == 1
                push!(cmp, ix[1])
            else
                push!(cmp, nothing)
            end
        end
        components[g] = cmp
    end

    return components
end

annots = Dict{String,Any}()
GZip.open(fn) do io
    ti = TarIterator(io, :file)
    for (h, iox) in ti
        p = h.path
        pp = splitext(p)
        println(pp[1])
        @assert length(pp) == 2 && pp[2] == ".xml"
        x = read(io, String)
        y = parse_annot(x)
        cmp = find_components(y)
        for (k, v) in cmp
            y["$(k)_components"] = v
        end
        annots[pp[1]] = y
    end
end

GZip.open("annotations.ser.gz", "w") do io
    serialize(io, annots)
end
