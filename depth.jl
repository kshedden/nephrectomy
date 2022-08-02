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

# Statistically compare the median depths for every pair of glom types.
function depth_analysis(depths)
    rslt = (Groups = String[], Z = Float64[], N = Int[])

    # Loop over all distinct pairs of glomerulus types, excluding "All_glomeruli"
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

            k1, k2, z = if z > 0
                j1, j2, z
            else
                j2, j1, -z
            end

            push!(rslt.Groups, @sprintf("%s - %s", glom_types[k2], glom_types[k1]))
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

function boundary_depth(v1::StaticVector{2}, v2::Vector{StaticVector{2}})

    dmin::Float64 = Inf

    for c in v2
        dmin = min(dmin, norm(v1 - c))
    end
    return dmin
end

# Calculate the median depth for each glom type in nephrectomy 'neph' using the given depth function.
function run_depth(neph, depth, ref)

    r = Vector{Union{Missing,Float64}}(undef, length(glom_types))
    r .= missing

    if length(ref) < 5
        return r
    end

    for k in keys(neph)
        if k in glom_types
            di = [depth(z, ref) for z in neph[k]]
            r[gti[k]] = length(di) > 0 ? median(di) : missing
        end
    end
    return r
end

# Get the median depth for each glom type within each nephrectomy
function get_depths(depth, rkey)
    idx, res = [], []
    for (k, neph) in annotsx
        push!(idx, k)

        # Compute depth relative to these points.
        ref = get(neph, rkey, [])
        if rkey in boundary_types
            # Convert vector of arrays to vetor of 2-vectors.
            v = Vector{StaticVector{2}}()
            for b in ref
                for c in eachcol(b)
                    push!(v, SVector{2}(c[1], c[2]))
                end
            end
            ref = v
        end

        r = run_depth(neph, depth, ref)
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

function dfclean(df)
    x = split(string(df), "\n")
    x = vcat(x[2], x[4:end])
    x = join(x, "\n")
    return x
end

out = open("depth.txt", "w")

rkeys = ["Normal", "Normal", "Capsule", "CMJ"]

dn = ["L2 depth", "Tukey halfspace depth", "Distance to capsule", "Distance to CMJ"]

for (jd, depth) in enumerate([l2_depth, tukey_depth, boundary_depth, boundary_depth])
    depths = get_depths(depth, rkeys[jd])
    drslt = depth_analysis(depths)
    crslt = clinical_analysis(depths)
    write(out, @sprintf("=== %s ===\n\n", dn[jd]))
    write(out, "Differences in median depth based on glomerulus type:\n")
    write(out, dfclean(drslt))
    write(out, "\n\nAssociations between depth and clinical variables:\n")
    write(out, dfclean(crslt))
    write(out, "\n\n")
end

close(out)
