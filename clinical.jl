using GZip, CSV, DataFrames, IterTools, LinearAlgebra, UnicodePlots, Printf, PyPlot

rm("plots", force = true, recursive = true)
mkdir("plots")

include("defs.jl")
include("pair_corr_utils.jl")
include("clinical_utils.jl")

# Pairwise correlation quantiles
# pcq are the log ratios of atypical/atypical distances versus typical/typical distances
# pcqn are the atypical/atypical distances
# pcqd are the typical/typical distances
idpcq, pcq, pcqn, pcqd = get_normalized_paircorr()

# Age has the most complete coverage so use it to save
# the scores and loadings.
function save_age()

    y, x, ids = get_response(:Age, idpcq, pcq)

    # PCA
    for j = 1:size(x, 2)
        x[:, j] .-= mean(x[:, j])
    end
    u, s, v = svd(x)

    u_df = DataFrame(
        :Scanner_id => ids[:, :Scanner_id],
        :TCP_id => ids[:, :TCP_id],
        :score1 => u[:, 1],
        :score2 => u[:, 2],
        :score3 => u[:, 3],
    )
    v_df = DataFrame(:factor1 => v[:, 1], :factor2 => v[:, 2], :factor3 => v[:, 3])
    CSV.write("age_scores.csv", u_df)
    CSV.write("age_factors.csv", v_df)
end

save_age()


function analyze(vname, ifig, out)

    y, x, ids = get_response(vname, idpcq, pcq)
    m = size(x, 2)

    # PCA
    for j = 1:size(x, 2)
        x[:, j] .-= mean(x[:, j])
    end
    u, s, v = svd(x)

    for j = 1:3
        if sum(v[:, j] .< 0) > sum(v[:, j] .> 0)
            v[:, j] = -v[:, j]
            u[:, j] = -u[:, j]
        end
    end

    # Check to see if the PC scores are correlated with
    # the clinical trait.
    c = Float64[]
    c = push!(c, length(y))
    for j = 1:3
        push!(c, cor(y, u[:, j]))
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

    # For age only, plot the distance quantiles
    if vname == :Age
        pr = collect(range(1 / m, 1 - 1 / m, length = m))
        for (i, id) in enumerate(idpcq)
            PyPlot.clf()
            PyPlot.axes([0.1, 0.15, 0.75, 0.8])
            PyPlot.grid(true)
            PyPlot.title(id)
            PyPlot.plot(pr, pcqn[i, :] / 1000, label = "AA")
            PyPlot.plot(pr, pcqd[i, :] / 1000, label = "TT")
            ha, lb = PyPlot.gca().get_legend_handles_labels()
            leg = PyPlot.figlegend(ha, lb, "center right")
            leg.draw_frame(false)
            PyPlot.xlabel("Probability", size = 15)
            PyPlot.ylabel("Quantile", size = 15)
            PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
            ifig += 1
        end
    end

    println(c)
    write(out, @sprintf("%s,%d,%f,%f,%f\n", vname, c...))

    return ifig

end

function main()
    ifig = 0
    out = open("clinical_results.csv", "w")
    write(out, "Variable,N,R1,R2,R3\n")
    for av in avn
        ifig = analyze(av, ifig, out)
    end
    close(out)
    return ifig
end

ifig = main()

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=clinical_loadings.pdf $f`
run(c)
