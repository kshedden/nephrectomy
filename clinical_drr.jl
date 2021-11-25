#=
Use dimension reduction regression (DRR) to relate clustering patterns of
atypical and typical glomeruli to clinical endpoints.
=#

using GZip, CSV, DataFrames, IterTools, LinearAlgebra, UnicodePlots, Printf, PyPlot, Dimred

rm("plots", force = true, recursive = true)
mkdir("plots")

include("defs.jl")
include("pair_corr_utils.jl")
include("clinical_utils.jl")

#= 
Pairwise correlation quantiles
  pcq are the log ratios of atypical/atypical distances versus typical/typical distances
  pcqn are the atypical/atypical distances
  pcqd are the typical/typical distances
=#
idpcq, pcq, pcqn, pcqd = get_normalized_paircorr()

function analyze(vname, ifig, out)

    y, x, ids = get_response(vname, idpcq, pcq)
    n, m = size(x)

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
    nslice = 20
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

    kt = knockoff_test(y, xx, drm)
    if kt.Pvalues[1] < 0.1
        println(lineplot(cf[:, 1]))
        println(kt.Pvalues[:, 1])
    end

end

function main()
    ifig = 0
    out = open("clinical_results_drr.csv", "w")
    write(out, "Variable,N,R1,R2,R3\n")
    for av in avn
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
