using CodecZlib, CSV, DataFrames, LinearAlgebra
using Printf, StaticArrays, Statistics

include("defs.jl")
include("annot_utils.jl")
include("clinical_utils.jl")

annotsx = glom_centroids(annots)
annotsx = major_components(annotsx)

# Create a column index for each glom type
gti = Dict()
for (k, v) in enumerate(glom_types)
    gti[v] = k
end

function l2_depth(v1::StaticVector{2}, v2::Vector{StaticVector{2}})
    d = mean([norm(v1 - x) for x in v2])
    return 1e6 / (1 + d)
end

function tukey_depth(v1::StaticVector{2}, v2::Vector{StaticVector{2}}; steps::Int = 100)

    angle = range(0, pi * (1 - 1 / steps), length = steps)
    vv = [v - v1 for v in v2]

    mn = length(v2)
    for a in angle
        m = 0
        for v in vv
            b = cos(a) * v[1] + sin(a) * v[2]
            m += b > 0 ? 1 : 0
        end
        mn = min(mn, m, length(v2) - m)
    end

    return mn / length(v2)
end

# Calculate the median depth for each glom type in nephrectomy 'neph' using the given depth function.
function run_depth(neph, depth)

    # Compute depth relative to these points.
    rf = neph["Normal"]

    r = Vector{Union{Missing,Float64}}(undef, length(glom_types))
    r .= missing

    if length(rf) < 5
        return r
    end

    for k in keys(neph)
        if k in glom_types
            di = [depth(z, rf) for z in neph[k]]
            r[gti[k]] = length(di) > 0 ? median(di) : missing
        end
    end
    return r
end

# Get the median depth for each glom type within each nephrectomy
function get_depths()
    idx, res = [], []
    for (k, v) in annotsx
        push!(idx, k)
        r = run_depth(v, l2_depth)
        push!(res, r)
    end

    depths = hcat(res...)'
    depths = DataFrame(depths, glom_types)
    depths[:, :Scanner_ID] = [parse(Int, x) for x in idx]

    clin[!, :Race] = [ismissing(x) ? missing : Int(x) for x in clin[:, :Race]]
    for x in unique(clin[:, :Race])
        clin[:, Symbol("Race$(x)")] = clin[:, :Race] .== x
    end
    depths = leftjoin(depths, clin, on = :Scanner_ID)

    return depths
end

# Statistically compare the median depths for every pair of glom types.
function depth_analysis(depths)
    rslt = (Groups = String[], Z = Float64[], N = Int[])
    for j1 = 1:length(glom_types)
        for j2 = 1:j1-1
            if "All_glomeruli" in [glom_types[j1], glom_types[j2]]
                continue
            end
            x = depths[:, [j1, j2]]
            x = x[completecases(x), :]
            if size(x, 1) == 0
                continue
            end
            di = x[:, 2] - x[:, 1]
            z = sqrt(length(di)) * mean(di) / std(di)
            push!(rslt.Groups, @sprintf("%s - %s", glom_types[j2], glom_types[j1]))
            push!(rslt.Z, z)
            push!(rslt.N, length(di))
        end
    end
    rslt = DataFrame(rslt)
    rslt = sort(rslt, :Z)
    return rslt
end

function clinical_analysis(depths)
    crslt =
        (Glomtype = String[], Clinvar = String[], r = Float64[], Z = Float64[], N = Int[])
    for c in names(clin)[6:end]
        if c in ["Race"]
            continue
        end
        for gt in glom_types
            if gt == "All_glomeruli"
                continue
            end
            xx = depths[:, [gt, c]]
            xx = xx[completecases(xx), :]
            if size(xx, 1) < 30
                continue
            end
            r = cor(xx[:, 1], xx[:, 2])
            rz = sqrt(size(xx, 1) - 3) * r
            push!(crslt.Glomtype, gt)
            push!(crslt.Clinvar, c)
            push!(crslt.r, r)
            push!(crslt.Z, rz)
            push!(crslt.N, size(xx, 1))
        end
    end
    crslt = DataFrame(crslt)
    crslt = sort(crslt, :Z)
    return crslt
end

depths = get_depths()
drslt = depth_analysis(depths)
crslt = clinical_analysis(depths)
