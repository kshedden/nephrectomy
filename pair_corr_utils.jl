using Statistics

# Calculate pair-correlations within and between
# the glomerular categories in the dictionary 'a'.
# If 'use' is present, exclude categories not listed
# in it.
function pair_corr(a; use::Array{String,1} = [])

    dit = Dict{Tuple{String,String},Array{Float64,1}}()

    ptp = collect(keys(a))
    if length(use) > 0
        ptp = [x for x in ptp if x in use]
    end

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
                d = sqrt((x1 - x2)^2 + (y1 - y2)^2)

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
function get_normalized_paircorr()

    # We are interested in the pairwise distances between these
    # two types
    k1 = ("Atypical", "Atypical")

    # Normalize to the pairwise distances between these two types
    k2 = ("All Glomeruli", "All Glomeruli")

    pc = Dict{String,Array{Float64,1}}()
    x, xn, xd, idx = [], [], [], []
    for fn in fi

        # Get the scanner ID from the file name
        fni = replace(fn, ".xml" => "")
        fni = parse(Int, fni)

        # Skip if we can't match this scanner id to a
        # TCP id.
        if !haskey(idmr, fni)
            continue
        end

        println(fn)
        a = read_annot(fn)

        # Condense to typical and atypical groups
        b = Dict{String,Array{Array{Float64,2},1}}()
        b["All Glomeruli"] = a["All Glomeruli"]
        b["Atypical"] = Array{Float64,1}()
        for x in ["FGGS", "Ischemic", "FGGS", "Imploding"]
            if haskey(a, x)
                push!(b["Atypical"], a[x]...)
            end
        end

        dit = pair_corr(b, use = ["All Glomeruli", "Atypical"])

        if (length(dit[k1]) == 0) || (length(dit[k2]) == 0)
            continue
        end

        push!(idx, fni)
        push!(x, log.(dit[k1] ./ dit[k2]))
        push!(xn, dit[k1])
        push!(xd, dit[k2])

    end

    x = hcat(x...)'
    xn = hcat(xn...)'
    xd = hcat(xd...)'

    return tuple(idx, x, xn, xd)

end

# Return a phenotype vector y for variable 'vname', and the corresponding
# array of normalized distance quantiles.
function get_response(vname, idpcq, pcq)

    # Keep track of the subjects that cannot be matched.
    out = open("nomatch.csv", "w")
    y = Union{Float64,Missing}[]
    ids = []
    for id in idpcq
        if haskey(tcp_rownum, idmr[id])
            ri = tcp_rownum[idmr[id]]
            push!(ids, [id, idmr[id]])
            push!(y, df[ri, vname])
        else
            write(out, @sprintf("%s,%s\n", id, idmr[id]))
            push!(y, missing)
        end
    end
    close(out)

    ii = [i for (i, v) in enumerate(y) if !ismissing(v)]
    y = Array{Float64}(y[ii])
    x = Array{Float64,2}(pcq[ii, :])

    ids = ids[ii]
    id_df = DataFrame(:Scanner_id => [x[1] for x in ids], :TCP_id => [x[2] for x in ids])

    return (y, x, id_df)
end
