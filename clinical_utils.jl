using GZip, CSV, DataFrames

df = GZip.open("clinical_var.csv.gz") do io
    CSV.read(io, DataFrame)
end

# The TCP id's mix hyphens and underscores, so standardize to use underscores.
df[!, :TCP_ID] = [replace(x, "-"=>"_") for x in df[:, :TCP_ID]]

# Map between id's used in the clinical and annotation files
# idm maps TCP ID to Scanner ID
# idmr maps Scanner ID to TCP ID
idm, idmr = GZip.open("id_map.csv.gz") do io
    dx = CSV.read(io, DataFrame)
    idm = Dict{String,Int}()
    idmr = Dict{Int,String}()
    for r in eachrow(dx)
        if !ismissing(r.Scanner_ID)
            tid = replace(r.TCP_ID, "-"=>"_")
            sid = Int(r.Scanner_ID)
            idm[tid] = sid
            idmr[sid] = tid
        end
    end
    idm, idmr
end

# Scanner id's
df[:, :sid] = [haskey(idm, y) ? idm[y] : missing for y in df[:, :TCP_ID]]

# Map from TCP id to row number in df
tcp_rownum = Dict{String,Int}()
for (i, r) in enumerate(eachrow(df))
    tcp_rownum[r.TCP_ID] = i
end

# Map from annotation id's to row indices in the clinical data.
sid_rownum = Dict{Int,Int}()
for (k,v) in idmr
    if haskey(tcp_rownum, v)
        sid_rownum[k] = tcp_rownum[v]
    end
end
