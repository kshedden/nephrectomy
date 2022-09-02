using Serialization, CodecZlib, StaticArrays, Statistics

include("defs.jl")

# Read the annotation information into a dictionary.
# annots maps each sample's id to its annotation data.
pa = "/home/kshedden/data/Markus_Bitzer/Annotations"
annots = open(GzipDecompressorStream, joinpath(pa, "annotations.ser.gz")) do io
    deserialize(io)
end

# Reduce each glom bounding box to its centroid
function glom_centroids(annots)
    annotsx = Dict()
    for (k, g) in annots
        h = Dict{String,Any}()
        for q in keys(g)
            if q in glom_types
                u = StaticVector{2}[]
                for x in g[q]
                    z = mean(x, dims = 2)
                    push!(u, SVector{2,Float64}(z[1], z[2]))
                end
                h[q] = u
            else
                h[q] = g[q]
            end
        end
        annotsx[k] = h
    end
    return annotsx
end

# Collapse all atypical gloms into one class labeled "Atypical"
function condense(a)

    b = Dict{String,Any}()

    # Copy these without change
    for x in [
        "All_glomeruli",
        "All_glomeruli_components",
        "Tissue",
        "Capsule",
        "Cortex",
        "CMJ",
        "Normal",
        "Normal_components",
    ]
        if haskey(a, x)
            b[x] = a[x]
        else
            b[x] = Float64[]
        end
    end

    # Combine these glomerular types into one category
    b["Atypical"] = SVector{2}[]
    b["Atypical_components"] = []
    for x in atypical_glom_types
        if haskey(a, x)
            push!(b["Atypical"], a[x]...)
            push!(b["Atypical_components"], a["$(x)_components"]...)
        end
    end

    return b
end

# Retain only the glomeruli contained in the major connected component of each nephrectomy.
function major_components(annots)

    b = Dict()
    for k in keys(annots)
        b[k] = Dict()

        # Find the largest component
        c = annots[k]["All_glomeruli_components"]
        h = Dict()
        for x in c
            y = get(h, x, 0)
            y += 1
            h[x] = y
        end
        m = maximum(values(h))
        km = first([h[x] == m for x in keys(h)])

        for ky in keys(annots[k])
            if ky in glom_types && haskey(annots[k], "$(ky)_components")
                ii = [i for (i, c) in enumerate(annots[k]["$(ky)_components"]) if c == km]
                b[k][ky] = annots[k][ky][ii]
                b[k]["$(ky)_components"] = annots[k][ky][ii]
            else
                b[k][ky] = annots[k][ky]
            end
        end
    end

    return b
end
