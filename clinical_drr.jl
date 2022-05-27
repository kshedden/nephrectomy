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


#=
Pairwise correlation quantiles
  pcq are the log ratios of atypical/atypical distances versus typical/typical distances
  pcqn are the atypical/atypical distances
  pcqd are the typical/typical distances
=#
scid, pcq, pcqn, pcqd = get_normalized_paircorr(annots)

function analyze(vname, ifig, out)

    y, x, ids = get_response(vname, scid, pcq)
    n, m = size(x)

    if size(x, 1) < 20
        return
    end

    # Center the covariates
    for j = 1:m
        x[:, j] .-= mean(x[:, j])
    end

    # Functional PCA
    fm = zeros(m, m - 2)
    for j = 1:m-2
        fm[j:j+2, j] = [1, -2, 1]
    end
    w = 100.0 # Smoothing penalty weight
    a, tm = eigen(cov(x) - w * fm * fm')
    ii = sortperm(a, rev = true)
    a = a[ii]
    tm = tm[:, ii]
    tm = tm[:, 1:3]
    xx = x * tm

    # Use SIR or PHD to estimate the directions
    nslice = 10
    ndir = 3
    drm = (y, x) -> sir(y, x; nslice = nslice, ndir = ndir)
    #drm = (y,x)->phd(y, x; ndir=ndir)
    mf = drm(y, xx)

    # Get the coefficients and flip as needed so
    # that the linear predictors are positively
    # correlated with the trait.
    cf = tm * mf.dirs
    for j = 1:size(cf, 2)
        if cor(y, xx * mf.dirs[:, j]) < 0
            cf[:, j] .*= -1
        end
    end

    st = sir_test(mf)
    p = st.Pvalues
    write(out, @sprintf("%s,%d,%f,%f,%f\n", vname, length(y), p[1], p[2], p[3]))

    if st.Pvalues[1] < 0.1
        println(lineplot(cf[:, 1]))
    end

end

function main()
    ifig = 0
    out = open("clinical_results_drr.csv", "w")
    write(out, "Variable,N,P1,P2,P3\n")
    for av in names(clin)[6:end]
        println(av)
        ifig = analyze(av, ifig, out)
    end
    close(out)
    return ifig
end

ifig = main()

#f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
#c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=clinical_loadings.pdf $f`
#run(c)
