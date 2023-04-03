# Statistically compare a collection of depth summaries between every
# pair of glomerulus types in the vector 'gt'.
function depth_analysis(depths, gt::Vector{String})

    # The results will go here
    rslt = (Groups = String[], Z = Float64[], N = Int[])

    # Loop over all distinct pairs of glomerulus types, excluding "All_glomeruli"
    for j1 = 1:length(gt)
        for j2 = 1:j1-1

            # Don't use the all glom class.
            if occursin("All_glomeruli", gt[j1]) || occursin("All_glomeruli", gt[j2])
                continue
            end

            # Some classes may be missing
            if !((gt[j1] in names(depths)) && (gt[j2] in names(depths)))
                continue
            end

            # Extract the depths to be compared
            x = depths[:, [gt[j1], gt[j2]]]
            x = x[completecases(x), :]
            if size(x, 1) == 0
                continue
            end

            # Paired difference Z-score
            di = x[:, 2] - x[:, 1]
            z = sqrt(length(di)) * mean(di) / std(di)

            # Reorder to make the Z-scores positive.
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

function clinical_analysis(depths, gt::Vector{String})
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

            # Require data for at least 30 people to take a correlation.
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

# Calculate the L2 depth of the vector 'v1' relative to the vectors
# in 'v2'.
function l2_depth(v1::StaticVector{2}, v2::Vector{StaticVector{2}})
    d = mean([norm(v1 - x) for x in v2])
    return 1e6 / (1 + d)
end

# Calculate the spatial depth.
function spatial_depth(v1::StaticVector{2}, v2::Vector{StaticVector{2}})

    u = zeros(2)
    s = zeros(2)

    n = 0
    for v in v2
        u .= v - v1
        nu = norm(u)
        if nu > 1e-10
            u ./= nu
            n += 1
        end
        s .+= u
    end
    s ./= n
    return 1 - norm(s)
end

# Calculate the Tukey depth of the vector 'v1' relative to the vectors in 'v2'.
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

# Calculate a quantile of the depths for each glom type in nephrectomy 'neph' using
# the given depth function.
function calc_depth_quantile(neph, depthfun, ref, gt::Vector{String}; pp = 0.5)

    # Map glom types to a quantile of the depths for that glom type.
    r = Dict()

    # Not enough reference gloms to compute depths against
    if length(ref) < 10
        return r
    end

    for k in gt

        if !haskey(neph, k)
            continue
        end

        di = [depthfun(z, ref) for z in neph[k]]
        m = length(di)
        if m > 0
            r[k] = m > 0 ? [m, quantile(di, pp)] : [0, missing]
        end
    end

    return r
end

# Create a dataframe from a vector of dictionaries.
function make_df(res, pp::Float64)

    # Glom types
    ky = union([keys(r) for r in res]...)
    ky = [k for k in ky]
    sort!(ky)
    kyi = Dict([k => i for (i, k) in enumerate(ky)])

    n = length(res)
    p = length(ky)
    x = Matrix{Union{Missing,Float64}}(missing, n, 2*p)

    for (j, r) in enumerate(res)
        for (k, v) in r
            x[j, 2*kyi[k]-1] = v[1]
            x[j, 2*kyi[k]] = v[2]
        end
    end

    # Make names for the dataframe columns
    na = []
    for k in ky
        push!(na, @sprintf("%s_n", k))
        push!(na, @sprintf("%s_q%02.0f", k, 100 * pp))
    end

    return DataFrame(x, na)
end

# Get the pp'th quantile of the depths for the glomeruli of each type
# within each nephrectomy
function get_depth_quantile(depthfun, ref::String, annots, gt; pp::Float64 = 0.5)
    idx, res = [], []

    # Loop over nephrectomy samples
    for (k, neph) in annots

        # Sample id
        push!(idx, k)

        # Compute depth relative to these points.
        refdata = get(neph, ref, [])
        if ref in boundary_types
            # Convert vector of arrays to vector of 2-vectors.
            v = Vector{StaticVector{2}}()
            for b in ref
                for c in eachcol(b)
                    push!(v, SVector{2}(c[1], c[2]))
                end
            end
            ref = v
        end

        r = calc_depth_quantile(neph, depthfun, refdata, gt; pp = pp)
        push!(res, r)
    end

    depths = make_df(res, pp)
    depths[:, :ID] = idx
    depths[:, :Scanner_ID] = [parse(Int, first(split(x, "_"))) for x in idx]

    return depths
end

# Create a matrix containing the pp'th quantile of depths
# calculated using 'depthfun' relative to the given reference
# class, using the data in 'annots'.
function build_depths(
    pp::Vector{Float64},
    depthfun,
    annots,
    gt::Vector{String};
    ref::String = "Normal",
)
    dd = nothing
    for p in pp
        dd1 = get_depth_quantile(depthfun, ref, annots, gt; pp = p)
        if isnothing(dd)
            dd = dd1
        else
            dd = leftjoin(dd, dd1, on = :Scanner_ID)
        end
    end
    return dd
end
