using Statistics

# Calculate pair-correlations within and between
# the glomerular categories in the dictionary 'a'.
# If 'use' is present, exclude categories not listed
# in it.  Only include distances in the same tissue
# component.
function pair_corr(a; use::Vector{String} = String[])

    dit = Dict{Tuple{String,String},Vector{Float64}}()

    ptp = collect(keys(a))
    ptp = [x for x in ptp if !endswith(x, "_components")]
    if length(use) > 0
        ptp = [x for x in ptp if x in use]
    end

    # Loop over all distinct pairs of types
    m = length(ptp)
    nreject = 0
    for jj in product(1:m, 1:m)

        j1, j2 = jj[1], jj[2]
        if j2 > j1
            continue
        end

        # Category names for the two categories being assessed
        # for pair correlation.
        q1, q2 = ptp[j1], ptp[j2]

        # Create a string key for the pair of types
        qq = tuple(q1, q2)
        if !haskey(dit, qq)
            dit[qq] = Float64[]
        end

        # The component label for each glom.
        cmp1 = a["$(q1)_components"]
        cmp2 = a["$(q2)_components"]

        # The centroids for each category.
        u1, u2 = a[q1], a[q2]

        for (k1, v1) in enumerate(u1)
            for (k2, v2) in enumerate(u2)

                # When comparing within a category, only assess each pair of gloms
                # once.
                if (j1 == j2) && (k2 >= k1)
                    break
                end

                # Only compare glom pairs in the same component.
                # Component 'nothing' is a glom that was not bounded
                # by any cortex loop.
                if cmp1[k1] != cmp2[k2] || isnothing(cmp1[k1]) || isnothing(cmp2[k2])
                    continue
                end

                # Distance between two glomeruli
                d = norm(v1 - v2)

                # If the gloms are too close, it is probably the same
                # glom listed twice so exclude.
                if d > 100
                    push!(dit[qq], d)
                else
                    nreject += 1
                end
            end
        end
    end

    if nreject > 0
        println("$(nreject) glom pairs were too close.")
    end

    # Reduce each list of correlations to a set of quantiles
    m = 20
    pp = collect(range(1 / m, 1 - 1 / m, length = m))
    for (k, v) in dit
        if length(v) >= 20
            dit[k] = [quantile(v, p) for p in pp]
        else
            dit[k] = [] # too little data to work with
        end
    end

    return dit
end

# Get the normalized pairwise distance quantiles.
function get_normalized_paircorr(annots)

    # We are interested in the pairwise distances between these
    # two types
    k1 = ("Atypical", "Atypical")

    # Normalize to the pairwise distances between these two types
    k2 = ("Normal", "Normal")

    pc = Dict{String,Vector{Float64}}()
    x, xn, xd, idx = [], [], [], []
    for neph_id in keys(annots)

        # Get the scanner ID from the file name
        fni = parse(Int, neph_id)

        a = annots[neph_id]
        dit = pair_corr(a, use = ["Normal", "Atypical"])

        if !(
            haskey(dit, k1) &&
            (length(dit[k1]) == 20) &&
            haskey(dit, k2) &&
            (length(dit[k2]) == 20)
        )
            println("skipping $(neph_id)")
            continue
        end

        push!(idx, fni)
        push!(x, log.(dit[k1] ./ dit[k2]))
        push!(xn, dit[k1])
        push!(xd, dit[k2])
    end

    x = copy(hcat(x...)')
    xn = copy(hcat(xn...)')
    xd = copy(hcat(xd...)')

    return tuple(idx, x, xn, xd)
end
