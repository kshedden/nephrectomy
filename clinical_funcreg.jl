using GZip, CSV, DataFrames, IterTools, LinearAlgebra, UnicodePlots, Printf
using PyPlot, Distributions, Random

rm("plots", force = true, recursive = true)
mkdir("plots")

include("defs.jl")
include("annot_utils.jl")
include("pair_corr_utils.jl")
include("clinical_utils.jl")

# Pairwise correlation quantiles
idpcq, pcq, pcqn, pcqd = get_normalized_paircorr(annots)

function funcreg(pm, xtx, x, y)

    b = pm \ (x' * y)

    # Fitted values
    yh = x * b

    # Approximate degrees of freedom
    degf = tr(pm \ xtx)

    # Unexplained variance
    sig2 = sum((yh - y) .^ 2) / (length(y) - degf)

    # Covariance matrix of parameter estimates
    cm = sig2 * (pm \ xtx / pm)

    # Chi-square statistic
    cs = b' * (cm \ b)

    # P-values
    pv = 1 - cdf(Chisq(20), cs)

    return (yh, b, cs, pv)
end

function analyze(vname, ifig, out)

    y, x = get_response(vname, idpcq, pcq)
    n, p = size(x)

    # Center the covariates and standardize the dependent variable.
    for j = 1:p
        x[:, j] .-= mean(x[:, j])
    end
    y = (y .- mean(y)) ./ std(y)

    # The second derivative penalty matrix
    F = zeros(p - 2, p)
    for i = 1:p-2
        F[i, i:i+2] = [1, -2, 1]
    end
    G = F' * F

    println(vname)
    bl = []
    for la in [0.01, 0.1, 1.0, 10.0, 100.0]

        # Functional regression coefficients
        xtx = x' * x
        pm = xtx + la * G

        yh, b, cs, pv = funcreg(pm, xtx, x, y)
        push!(bl, b)

		# Permutation p-value for the global null E[b] = 0.
		y0 = copy(y)
		nrep = 1000
		csa = zeros(nrep)
		for k in 1:nrep
			shuffle!(y0)
	        _, _, cs0, _ = funcreg(pm, xtx, x, y0)
	        csa[k] = cs0
		end
		pv1 = mean(csa .>= cs)

        r = cor(yh, y)
        write(out, @sprintf("%s,%.2f,%.3f,%.3f,%.3f\n", vname, la, pv, pv1, r))
    end

    # Plot the PC loading vector
    PyPlot.clf()
    PyPlot.axes([0.13, 0.12, 0.75, 0.8])
    PyPlot.grid(true)
    pr = collect(range(0, 1, length = size(x, 2)))
    pr = pr ./ maximum(pr)
    cm = PyPlot.cm.get_cmap("jet")
    for (j, b) in enumerate(bl)
        col = cm(j / (length(bl) + 1))
        PyPlot.plot(pr, b, "-", color = col)
    end
    PyPlot.title(vname)
    PyPlot.xlabel("Probability point", size = 15)
    PyPlot.ylabel("Coefficient", size = 15)
    PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
    ifig += 1

    return ifig

end

function main()
    ifig = 0
    out = open("clinical_funcreg_results.csv", "w")
    write(out, "Variable,Lambda,Parameteric-P-value,Permutation-P-value,Correlation\n")
    for av in avn
        ifig = analyze(av, ifig, out)
    end
    close(out)
    return ifig
end

ifig = main()

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=clinical_reg.pdf $f`
run(c)
