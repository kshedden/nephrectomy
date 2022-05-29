#=
Use dimension reduction regression (DRR) to relate clustering patterns of
atypical and typical glomeruli to clinical endpoints.
=#

using CodecZlib, CSV, DataFrames, IterTools, LinearAlgebra
using UnicodePlots, Printf, PyPlot, Dimred

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

#=
Pairwise correlation quantiles
  pcq are the log ratios of atypical/atypical distances versus typical/typical distances
  pcqn are the atypical/atypical distances
  pcqd are the typical/typical distances
=#
scid, pcq, pcqn, pcqd = get_normalized_paircorr(annotsx)

function analyze(vname, ifig, out)

    y, x, ids = get_response(vname, scid, pcq)

    # Only use the lower quantiles
    x = x[:, 1:10]
    n, m = size(x)

    if size(x, 1) < 20
        return ifig
    end

    # Center the covariates
    for j = 1:m
        x[:, j] .-= mean(x[:, j])
    end

    # Functional PCA
    npc = 5
    fm = zeros(m, m - 2)
    for j = 1:m-2
        fm[j:j+2, j] = [1, -2, 1]
    end
    w = 10.0 # Smoothing penalty weight
    a, load = eigen(cov(x) - w * fm * fm')
    ii = sortperm(a, rev = true)
    a = a[ii]
    load = load[:, ii[1:npc]]
    xx = x * load

    # Use SIR or PHD to estimate the directions
    nslice = 10
    ndir = 3
    drm = (y, x) -> sir(y, x; nslice = nslice, ndir = ndir)
    #drm = (y,x)->phd(y, x; ndir=ndir)
    mf = drm(y, xx)

    # Get the coefficients and flip as needed so
    # that the linear predictors are positively
    # correlated with the trait.
    cf = load * mf.dirs
    for j = 1:size(cf, 2)
        if cor(y, xx * mf.dirs[:, j]) < 0
            cf[:, j] .*= -1
        end
    end

    # Plot the coefficient vector
    PyPlot.clf()
    PyPlot.axes([0.13, 0.12, 0.75, 0.8])
    PyPlot.grid(true)
    pr = collect(range(0, 1, length = 20))
    pr = pr ./ maximum(pr)
    pr = pr[1:10]
    for j = 1:size(cf, 2)
        PyPlot.plot(pr, cf[:, j], label = @sprintf("%d", j))
    end
    ha, lb = PyPlot.gca().get_legend_handles_labels()
    leg = PyPlot.figlegend(ha, lb, "center right")
    leg.draw_frame(false)
    PyPlot.title(vname)
    PyPlot.xlabel("Probability point", size = 15)
    PyPlot.ylabel("Coefficient", size = 15)
    PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
    ifig += 1

    # Significance tests
    st = sir_test(mf)
    p = st.Pvalues

    if st.Pvalues[1] < 0.1
        println(lineplot(cf[:, 1]))
    end

    write(out, @sprintf("%s,%d,%f,%f,%f\n", vname, length(y), p[1], p[2], p[3]))

    return ifig
end

function main()
    ifig = 0
    out = open("clinical_drr_results.csv", "w")
    write(out, "Variable,N,P1,P2,P3\n")

    # Loop through the clinical variables, which are outcomes for the regressions.
    for av in names(clin)[6:end]
        println(av)
        ifig = analyze(av, ifig, out)
    end
    close(out)
    return ifig
end

ifig = main()

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c =
    `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=clinical_drr_loadings.pdf $f`
run(c)
