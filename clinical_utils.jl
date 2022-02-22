using CodecZlib, CSV, DataFrames

df1 = open("/home/kshedden/data/Markus_Bitzer/clinical_var.csv.gz") do io
    CSV.read(GzipDecoderStream(io), DataFrame)
end

df2 = open("/home/kshedden/data/Markus_Bitzer/nephrectomy_characteristics.csv.gz") do io
    CSV.read(GzipDecoderStream(io), DataFrame)
end

# The TCP id's mix hyphens and underscores, so standardize to use underscores.
df1[!, :TCP_ID] = [replace(x, "-" => "_") for x in df1[:, :TCP_ID]]
df2[!, :TCP_id] = [replace(x, "-" => "_") for x in df2[:, :TCP_id]]

df = outerjoin(df1, df2, on = :TCP_ID => :TCP_id)

df[:, :Female] = [!ismissing(x) && x == "F" ? 1 : 0 for x in df[:, :SEX]]
df[:, :NonWhite] =
    [!ismissing(x) && x != "White or Caucasian" ? 1 : 0 for x in df[:, :RACE]]

# Map between id's used in the clinical and annotation files
# idm maps TCP ID to Scanner ID
# idmr maps Scanner ID to TCP ID
idm, idmr = open("id_map.csv.gz") do io
    dx = CSV.read(GzipDecoderStream(io), DataFrame)
    idm = Dict{String,Int}()
    idmr = Dict{Int,String}()
    for r in eachrow(dx)
        if !ismissing(r.Scanner_ID)
            tid = replace(r.TCP_ID, "-" => "_")
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
for (k, v) in idmr
    if haskey(tcp_rownum, v)
        sid_rownum[k] = tcp_rownum[v]
    end
end
