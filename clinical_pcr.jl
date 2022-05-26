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

# Pairwise correlation quantiles
# pcq are the log ratios of atypical/atypical distances versus typical/typical distances
# pcqn are the atypical/atypical distances
# pcqd are the typical/typical distances
scid, pcq, pcqn, pcqd = get_normalized_paircorr(annots)

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


function analyze(vname, ifig, out)

    y, x, ids = get_response(vname, scid, pcq)
    m = size(x, 2)
    if size(x, 1) < 20
        return ifig
    end

    # PCA
    for j = 1:size(x, 2)
        x[:, j] .-= mean(x[:, j])
    end
    u, s, v = svd(x)

    # Try to make the loadings mostly positive
    for j = 1:3
        if sum(v[:, j] .< 0) > sum(v[:, j] .> 0)
            v[:, j] .*= -1
            u[:, j] .*= -1
        end
    end

    # Obtain the correlation coefficients between the PC score and
    # the clinical trait.
    c = Float64[length(y),]
    for j = 1:3
        r, z, p = permcor(y, u[:, j])
        push!(c, r, z, p)
    end

    # Plot the PC loading vector
    PyPlot.clf()
    PyPlot.axes([0.13, 0.12, 0.75, 0.8])
    PyPlot.grid(true)
    pr = collect(range(0, 1, length = size(v, 1)))
    pr = pr ./ maximum(pr)
    for j = 1:3
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

    write(out, @sprintf("%s,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", vname, c...))

    return ifig
end

function main()
    ifig = 0
    out = open("clinical_pcr_results.csv", "w")
    write(out, "Variable,N,R1,Z1,P1,R2,Z2,P2,R3,Z3,P3\n")
    for av in names(clin)[6:end]
        ifig = analyze(av, ifig, out)
    end
    close(out)
    return ifig
end

ifig = main()

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=clinical_loadings.pdf $f`
run(c)
