using Statistics

# Calculate pair-correlations within and between
# the glomerular categories in the dictionary 'a'.
# If 'use' is present, exclude categories not listed
# in it.
function pair_corr(a; use::Array{String,1}=[])

    dit = Dict{Tuple{String,String},Array{Float64,1}}()

    ptp = [x for x in keys(a) if x in use]

    # Loop over all distinct pairs of types
    m = length(ptp)
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
            dit[qq] = []
        end

        # The bounding boxes for each category.
        u1, u2 = a[q1], a[q2]

        for (k1, v1) in enumerate(u1)

            # Reduce each bounding box in category 1 to its centroid
            x1 = mean(v1[1, :])
            y1 = mean(v1[2, :])

            for (k2, v2) in enumerate(u2)

                # When comparing within a category, only assess each pair of gloms
                # once.
                if (j1 == j2) && (k2 >= k1)
                    break
                end

                # Reduce each bounding box in category 2 to its centroid
                x2 = mean(v2[1, :])
                y2 = mean(v2[2, :])

                # Distance between two glomeruli
                d = sqrt((x1-x2)^2 + (y1-y2)^2)

                # If the gloms are too close, it is probably the same
                # glom listed twice so exclude.
                if d > 10
                    push!(dit[qq], d)
                end
            end
        end
    end

    # Reduce each list of correlations to a set of quantiles
    m = 20
    pp = collect(range(1/m, 1-1/m, length=m))
    for (k, v) in dit
        if length(v) >= 20
            dit[k] = [quantile(v, p) for p in pp]
        else
            dit[k] = [] # too little data to work with
        end
    end

    return dit

end
