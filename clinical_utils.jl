using CodecZlib, CSV, DataFrames

# Clinical variables
df1 = open("/home/kshedden/data/Markus_Bitzer/clinical_var_021022.csv.gz") do io
    CSV.read(GzipDecompressorStream(io), DataFrame)
end
df1 = rename(df1, "GENDER" => "SEX")
df1 = rename(df1, "CURRENT_AGE" => "AGE")

df1[:, :Smoking] = replace(
    df1[:, :SMOKING_STATUS],
    "Never Assessed" => Missing,
    "Never Smoker" => 0,
    "Passive Smoke Exposure -" => 1,
    "Former Smoker" => 2,
    "Current Some Day Smoker" => 3,
    "Current Every Day Smoker" => 4,
)

# Nephrectomy characteristics
df2 = open("/home/kshedden/data/Markus_Bitzer/nephrectomy_characteristics.csv.gz") do io
    CSV.read(GzipDecompressorStream(io), DataFrame)
end
df2 = rename(df2, "TCP_id" => "TCP_ID")

# The TCP id's mix hyphens and underscores, so standardize to use underscores.
df1[!, :TCP_ID] = [replace(x, "-" => "_") for x in df1[:, :TCP_ID]]
df2[!, :TCP_ID] = [replace(x, "-" => "_") for x in df2[:, :TCP_ID]]

df = outerjoin(df1, df2, on = :TCP_ID)

# Recode a few variables
df[:, :Female] = replace(df[:, :SEX] .== "Female", true => 1, false => 0)
df[:, :NonWhite] = replace(df[:, :RACE] .!= "White or Caucasian", true => 1, false => 0)

# Map between id's used in the clinical and annotation files
# idm maps TCP ID to Scanner ID
# idmr maps Scanner ID to TCP ID
idm, idmr = open("id_map.csv.gz") do io
    dx = CSV.read(GzipDecompressorStream(io), DataFrame)
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
