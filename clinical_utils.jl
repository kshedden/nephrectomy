using CodecZlib, CSV, DataFrames

# Clinical variables
pa = "/home/kshedden/data/Markus_Bitzer"
clin = open(GzipDecompressorStream, joinpath(pa, "Spatial_IDs.csv.gz")) do io
    CSV.read(io, DataFrame)
end

egfr = open(
    GzipDecompressorStream,
    joinpath(pa, "eGFR  - updated formula.xlsx - Sheet1.csv.gz"),
) do io
    CSV.read(io, DataFrame)
end
egfr = egfr[:, [:Precise_ID, :Baseline_eGFR]]
egfr = rename(egfr, :Baseline_eGFR => :Revised_eGFR)

clin = leftjoin(clin, egfr, on = :Precise_ID)

# Return a phenotype vector y for variable 'vname', and the corresponding
# array of normalized distance quantiles.
function get_response(vname, scid, pcq)

    # Map from scanner id to clinical variable value
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
    x = copy(hcat(x...)')

    return y, x, ids
end
