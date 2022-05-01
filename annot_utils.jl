using Serialization, CodecZlib

include("defs.jl")

# Read the annotations information into a dictionary.
# annots maps each sample's id to its annotation data.
annots = open(GzipDecompressorStream, "annotations.ser.gz") do io
    deserialize(io)
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
    b["Atypical"] = []
    b["Atypical_components"] = []
    for x in ["FSGS", "FGGS", "Ischemic", "Imploding"]
        if haskey(a, x)
            push!(b["Atypical"], a[x]...)
            push!(b["Atypical_components"], a["$(x)_components"]...)
        end
    end

    return b
end
