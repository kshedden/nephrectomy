using Printf, Statistics, IterTools, LinearAlgebra, PyPlot
using DataFrames
using Distributions
using CSV

rm("plots", force = true, recursive = true)
mkdir("plots")

include("defs.jl")
include("annot_utils.jl")
include("pair_corr_utils.jl")

out = open("nopheno_analysis.txt", "w")

annotsx = glom_centroids(annots)
annotsx = Dict(k => condense(v) for (k, v) in annotsx)

# Only consider samples with a minimum number of total glomeruli
annotsx = Dict(k => v for (k, v) in annotsx if length(v["All_glomeruli"]) >= 100)

# Pairwise correlation quantiles
# pcq are the log ratios of atypical/atypical distances versus typical/typical distances
# pcqn are the atypical/atypical distances
# pcqd are the typical/typical distances
scid, pcq, pcqn, pcqd = get_normalized_paircorr(annotsx)

m = 20
m2 = m * (m - 1) / 2
pp = collect(range(10 / m2, 1 - 10 / m2, length = m))

function compare_means(ifig)
    n = size(pcq, 1)
    z = sqrt(n) * mean(pcq, dims = 1) ./ std(pcq, dims = 1)

    mn = mean(pcq, dims = 1)[:]
    se = std(pcq, dims = 1)[:] / sqrt(n)
    z = mn./se
    pval = 2*cdf(Normal(0, 1), -abs.(z))

    odf = DataFrame(p=pp, mean=mn, sd=se, z=z, pval=pval, ratio=exp.(mn) .- 1, lcb=exp.(mn-2*se) .- 1, ucb=exp.(mn+2*se) .- 1)
    CSV.write("nopheno_analysis.csv", odf)

    PyPlot.clf()
    PyPlot.axes([0.15, 0.1, 0.8, 0.8])
    PyPlot.grid(true)
    PyPlot.plot(pp, mn, "-", color = "black")
    PyPlot.fill_between(pp, mn - 2 * se, mn + 2 * se, color = "grey")
    PyPlot.xlabel("Probability point", size = 15)
    PyPlot.ylabel(raw"$\log d_p(AA) / d_p(TT)$", size = 15)
    PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
    return ifig + 1
end

function plot_logr_dist(ifig)
    for k = 1:10
        PyPlot.clf()
        PyPlot.hist(pcq[:, k], ec = "black", fc = "white")
        PyPlot.xlabel(raw"$\log d_p(AA)/d_p(TT)$", size = 15)
        PyPlot.ylabel("Frequency", size = 15)
        PyPlot.title(@sprintf("p=%.2f", pp[k]))
        PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
        ifig += 1
    end
    return ifig
end

ifig = 0
ifig = compare_means(ifig)
ifig = plot_logr_dist(ifig)

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=nopheno_analysis.pdf $f`
run(c)

close(out)
