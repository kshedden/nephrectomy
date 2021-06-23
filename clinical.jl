using GZip, CSV, DataFrames, IterTools, LinearAlgebra, UnicodePlots, Printf

include("defs.jl")
include("pair_corr_utils.jl")
include("clinical_utils.jl")

# Get the normalized pairwise distance quantiles.
function get_paircorr()

    # We are interested in the pairwise distances between these
    # two types
    k1 = ("Atypical", "Atypical")

    # Normalize to the pairwise distances between these two types
    k2 = ("All Glomeruli", "All Glomeruli")

    pc = Dict{String,Array{Float64,1}}()
    x, idx = [], []
    for fn in fi

        # Get the scanner ID from the file name
        fni = replace(fn, ".xml"=>"")
        fni = parse(Int, fni)

        # Skip if we can't match this scanner id to a
        # TCP id.
        if !haskey(idmr, fni)
            continue
        end

        println(fn)
        a = read_annot(fn)

        b = Dict{String,Array{Array{Float64,2},1}}()
        b["All Glomeruli"] = a["All Glomeruli"]
        b["Atypical"] = Array{Float64,1}()
        for x in ["FGGS", "Ischemic", "FGGS", "Imploding"]
            if haskey(a, x)
                push!(b["Atypical"], a[x]...)
            end
        end

        dit = pair_corr(b, use=["All Glomeruli", "Atypical"])

        if (length(dit[k1]) == 0) || (length(dit[k2]) == 0)
            continue
        end

        push!(idx, fni)
        push!(x, log.(dit[k1] ./ dit[k2]))

    end

    x = hcat(x...)'

    return tuple(idx, x)

end

# Pairwise correlation quantiles
idpcq, pcq = get_paircorr()

# Return a phenotype vector y for variable 'vname', and the corresponding
# array of normalized distance quantiles.
function get_response(vname)

    xm = Dict{String,Float64}()
    for (i, r) in enumerate(eachrow(df))
        if !ismissing(r[vname])
            xm[r.TCP_ID] = r[vname]
        end
    end

    # Keep track of the subject that cannot be matched.
    out = open("nomatch.csv", "w")

    y = []
    for id in idpcq
        if haskey(tcp_rownum, idmr[id])
            ri = tcp_rownum[idmr[id]]
            push!(y, df[ri, vname])
        else
            write(out, @sprintf("%s,%s\n", id, idmr[id]))
            push!(y, missing)
        end
    end

    close(out)

    ii = [i for (i,v) in enumerate(y) if !ismissing(v)]
    return tuple(y[ii], pcq[ii, :])

end

function analyze(vname)

    y, x = get_response(vname)

    # PCA
    for j in 1:size(x, 2)
        x[:, j] .-= mean(x[:, j])
    end
    u,s,v = svd(x)

    # Check to see if the PC scores are correlated with
    # the clinical trait.
    for j in 1:3
        println(cor(y, u[:, j]))

        # Plot the PC loading vector
        plt = lineplot(v[:, j])
        println(plt)
    end

end

analyze(:bmi)
