using CodecZlib, CSV, DataFrames

# Clinical variables
clin = open(
    GzipDecompressorStream,
    "/home/kshedden/data/Markus_Bitzer/Spatial_IDs.csv.gz",
) do io
    CSV.read(io, DataFrame)
end

# Return a phenotype vector y for variable 'vname', and the corresponding
# array of normalized distance quantiles.
function get_response(vname, scid, pcq)

    println(vname)

    cm = Dict()
    for r in eachrow(clin)
        if !ismissing(r[vname])
            cm[r.Scanner_ID] = r[vname]
        end
    end

    y, x, ids = [], [], []
    for (j, id) in enumerate(scid)
        if haskey(cm, id)
            push!(ids, id)
            push!(y, cm[id])
            push!(x, pcq[j, :])
        end
    end
    x = hcat(x...)'

    return y, x, ids
end
