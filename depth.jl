using DataFrames, LinearAlgebra, Printf, StaticArrays, Statistics

include("defs.jl")
include("annot_utils.jl")
include("clinical_utils.jl")

# Statistically compare the median depths for every pair of glom types.
function depth_analysis(depths, gt)
    rslt = (Groups = String[], Z = Float64[], N = Int[])

    # Loop over all distinct pairs of glomerulus types, excluding "All_glomeruli"
    for j1 = 1:length(gt)
        for j2 = 1:j1-1
            if occursin("All_glomeruli", gt[j1]) || occursin("All_glomeruli", gt[j2])
                continue
            end
            if !((gt[j1] in names(depths)) && (gt[j2] in names(depths)))
                continue
            end
            x = depths[:, [gt[j1], gt[j2]]]
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

            push!(rslt.Groups, @sprintf("%s - %s", gt[k2], gt[k1]))
            push!(rslt.Z, z)
            push!(rslt.N, length(di))
        end
    end
    rslt = DataFrame(rslt)
    rslt = sort(rslt, :Z)
    return rslt
end

function clinical_analysis(depths, gt)
    crslt =
        (Glomtype = String[], Clinvar = String[], r = Float64[], Z = Float64[], N = Int[])

    for c in names(clin)[6:end]
        if occursin("Race", c)
            continue
        end
        for gg in gt

            if !(gg in names(depths))
                continue
            end

            xx = depths[:, [gg, c]]
            xx = xx[completecases(xx), :]
            if size(xx, 1) < 30
                continue
            end
            r = cor(xx[:, 1], xx[:, 2])
            rz = sqrt(size(xx, 1) - 3) * r
            push!(crslt.Glomtype, gg)
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

# Calculate a quantile of the depths for each glom type in nephrectomy 'neph' using the given depth function.
function calc_depth_quantile(neph, depth, ref; pp = 0.5)

    r = Dict{String,Union{Missing,Float64}}()

    # Not enough reference gloms to compute depths against
    if length(ref) < 5
        return r
    end

    for k in keys(neph)
        if k in glom_types || k == "Atypical"
            di = [depth(z, ref) for z in neph[k]]
            if length(di) > 0
                r[k] = length(di) > 0 ? quantile(di, pp) : missing
            end
        end
    end

    return r
end

function make_df(res, pp)

    # Glom types
    ky = union([keys(r) for r in res]...)
    ky = Vector([k for k in ky])
    sort!(ky)
    kyi = Dict([k => i for (i, k) in enumerate(ky)])

    n = length(res)
    p = length(ky)
    x = Matrix{Union{Missing,Float64}}(missing, n, p)

    for (j, r) in enumerate(res)
        for (k, v) in r
            x[j, kyi[k]] = v
        end
    end

    return DataFrame(x, [@sprintf("%s_%.2f", k, pp) for k in ky])
end

# Get a quantile of the depths for the gloms of each type within each nephrectomy
function get_depth_quantile(depth, rkey, annots; pp = 0.5)
    idx, res = [], []

    # Loop over nephrectomy samples
    for (k, neph) in annots
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

        r = calc_depth_quantile(neph, depth, ref; pp = pp)
        push!(res, r)
    end

    depths = make_df(res, pp)
    depths[:, :Scanner_ID] = [parse(Int, x) for x in idx]

    return depths
end

function build_depths(pp, depthfun, annots)
    dd = nothing
    for p in pp
        dd1 = get_depth_quantile(depthfun, "Normal", annots; pp = p)
        if isnothing(dd)
            dd = dd1
        else
            dd = leftjoin(dd, dd1, on = :Scanner_ID)
        end
    end
    return dd
end

function dfclean(df)
    x = split(string(df), "\n")
    x = vcat(x[2], x[4:end])
    x = join(x, "\n")
    return x
end

function main(annots)
    out = open("depth.txt", "w")

    pp = 0.5
    gt = [@sprintf("%s_%.2f", g, pp) for g in glom_types]

    rkeys = ["Normal"]#, "Normal", "Capsule", "CMJ"]

    dn = ["L2 depth", "Tukey halfspace depth"]#, "Distance to capsule", "Distance to CMJ"]
    dpf = [l2_depth, tukey_depth] #, boundary_depth, boundary_depth]

    for (jd, depthfun) in enumerate(dpf)

        depths = build_depths([0.5], depthfun, annots)

        clin[!, :Race] = [ismissing(x) ? missing : Int(x) for x in clin[:, :Race]]
        for x in unique(clin[:, :Race])
            clin[:, Symbol("Race$(x)")] = clin[:, :Race] .== x
        end
        depths = leftjoin(depths, clin, on = :Scanner_ID)

        drslt = depth_analysis(depths, gt)
        crslt = clinical_analysis(depths, gt)
        write(out, @sprintf("=== %s ===\n\n", dn[jd]))
        write(out, "Differences in median depth based on glomerulus type:\n")
        write(out, dfclean(drslt))
        write(out, "\n\nAssociations between depth and clinical variables:\n")
        write(out, dfclean(crslt))
        write(out, "\n\n")
    end

    close(out)
end

annotsx = glom_centroids(annots)
annotsx = major_components(annotsx)
#annotsc = Dict(k => condense(v) for (k, v) in annotsx)

main(annotsx)
