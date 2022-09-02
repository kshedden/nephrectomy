#=
Use Principal Components Regression (PCR) to relate clinical
and spatial characteristics.
=#

using CSV, DataFrames, IterTools, LinearAlgebra, UnicodePlots
using Distributions, Printf, PyPlot, Random

rm("plots", force = true, recursive = true)
mkdir("plots")

include("defs.jl")
include("annot_utils.jl")
include("pair_corr_utils.jl")
include("clinical_utils.jl")

annotsx = glom_centroids(annots)
annotsx = Dict(k => condense(v) for (k, v) in annotsx)

# Only consider samples with a minimum number of total glomeruli
annotsx = Dict(k => v for (k, v) in annotsx if length(v["All_glomeruli"]) > 100)

# Pairwise correlation quantiles
# pcq are the log ratios of atypical/atypical distances versus typical/typical distances
# pcqn are the atypical/atypical distances
# pcqd are the typical/typical distances
scid, pcq, pcqn, pcqd = get_normalized_paircorr(annotsx)

# Permutation test for correlation coefficients
function permcor(y, x; nrep = 1000)

    y = copy(y)
    n = length(y)
    r = cor(y, x)
    z = sqrt(n - 3) * log((1 + r) / (1 - r)) / 2

    zl = zeros(nrep)
    for i = 1:nrep
        shuffle!(y)
        r1 = cor(y, x)
        z1 = sqrt(n - 3) * log((1 + r1) / (1 - r1)) / 2
        zl[i] = z1
    end

    return (r, z, mean(abs.(zl) .> abs(z)))
end


function analyze(vname, ifig, out; ncomp = 4)

    y, x, ids = get_response(vname, scid, pcq)

    if size(x, 1) < 20
        return ifig
    end

    # Keep only the lower quantiles
    x = x[:, 1:10]

    m = size(x, 2)

    # PCA
    for j = 1:size(x, 2)
        x[:, j] .-= mean(x[:, j])
    end
    u, s, v = svd(x)

    # Try to make the loadings mostly positive
    for j = 1:ncomp
        if sum(v[:, j] .< 0) > sum(v[:, j] .> 0)
            v[:, j] .*= -1
            u[:, j] .*= -1
        end
    end

    # Obtain the correlation coefficients between the PC score and
    # the clinical trait.
    c = Float64[length(y),]
    for j = 1:ncomp
        r, z, p = permcor(y, u[:, j])
        push!(c, r, z, p)
    end

    # Plot the PC loading vector
    PyPlot.clf()
    PyPlot.axes([0.13, 0.12, 0.75, 0.8])
    PyPlot.grid(true)
    pr = collect(range(0, 1, length = size(v, 1)))
    pr = pr ./ maximum(pr)
    for j = 1:ncomp
        PyPlot.plot(pr, v[:, j], label = @sprintf("%d", j))
    end
    ha, lb = PyPlot.gca().get_legend_handles_labels()
    leg = PyPlot.figlegend(ha, lb, "center right")
    leg.draw_frame(false)
    PyPlot.title(vname)
    PyPlot.xlabel("Probability point", size = 15)
    PyPlot.ylabel("Loading", size = 15)
    PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
    ifig += 1

    row = [vname, @sprintf("%d", c[1])]
    row = vcat(row, [@sprintf("%.4f", x) for x in c[2:end]])
    write(out, join(row, ","))
    write(out, "\n")

    return ifig
end

function main()
    ncomp = 4
    ifig = 0
    out = open("clinical_pcr_results.csv", "w")
    head = "Variable,N," * join(["R$(j),Z$(j),P$(j)" for j = 1:ncomp], ",")
    write(out, head)
    write(out, "\n")
    for av in names(clin)[6:end]
        ifig = analyze(av, ifig, out; ncomp = ncomp)
    end
    close(out)
    return ifig
end

ifig = main()

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=clinical_pcr_loadings.pdf $f`
run(c)
