using GZip, CSV, DataFrames, IterTools, LinearAlgebra, UnicodePlots, Printf
using PyPlot, Distributions

rm("plots", force = true, recursive = true)
mkdir("plots")

include("defs.jl")
include("annot_utils.jl")
include("pair_corr_utils.jl")
include("clinical_utils.jl")

# Pairwise correlation quantiles
idpcq, pcq, pcqn, pcqd = get_normalized_paircorr(annots)

function analyze(vname, ifig, out)

    y, x = get_response(vname, idpcq, pcq)
    n, p = size(x)

    for j = 1:p
        x[:, j] = x[:, j] .- mean(x[:, j])
    end
    y = (y .- mean(y)) ./ std(y)

    u, s, v = svd(x)
    sd = std(u[:, 1])

    # The second derivative penalty matrix
    F = zeros(p - 2, p)
    for i = 1:p-2
        F[i, i:i+2] = [1, -2, 1]
    end
    G = F' * F

    println(vname)
    bl = []
    for la in [0.1, 1.0, 10.0, 100.0]

        # Functional regression coefficients
        pm = x' * x + la * G
        b = pm \ (x' * y)
        push!(bl, b)

        # Fitted values
        yh = x * b

        # Unexplained variance
        sig2 = mean((yh - y) .^ 2)

        # Covariance matrix of parameter estimates
        cm = sig2 * (pm \ (x' * x) / pm)

        # Chi-square statistic
        cs = b' * (cm \ b)

        # P-values
        pv = 1 - cdf(Chisq(20), cs)
        r = cor(yh, y)
        write(out, @sprintf("%s,%.2f,%.3f,%.3f\n", vname, la, pv, r))
    end

    # Plot the PC loading vector
    PyPlot.clf()
    PyPlot.axes([0.13, 0.12, 0.75, 0.8])
    PyPlot.grid(true)
    pr = collect(range(0, 1, length = size(v, 1)))
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
    out = open("clinical_reg_results.csv", "w")
    write(out, "Variable,Lambda,P-value,Correlation\n")
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
