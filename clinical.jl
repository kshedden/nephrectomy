using GZip, CSV, DataFrames, IterTools, LinearAlgebra, UnicodePlots, Printf, PyPlot

rm("plots", force = true, recursive = true)
mkdir("plots")

include("defs.jl")
include("pair_corr_utils.jl")
include("clinical_utils.jl")

# Pairwise correlation quantiles
idpcq, pcq = get_normalized_paircorr()


function analyze(vname, ifig, out)

    y, x = get_response(vname, idpcq, pcq)

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

    println(c)
    write(out, @sprintf("%s,%d,%f,%f,%f\n", vname, c...))

    return ifig

end

# Variables to analyze
avn = [
    :Female,
    :NonWhite,
    :Age,
    :htn_code_pre,
    :dm_code_pre,
    :bmi,
    :htn_sbp,
    :htn,
    :dm_lab_pre,
    :dm,
    :bl_cr_prenx,
    :bl_egfr_prenx,
    :bl_ckd_prenx,
    :kdigo,
    :"Hyper-Bin",
    :"T2D-Bin",
    :"eGFR Min4 Mean",
    :"I/L Ratio",
    :">100%",
    :">50%",
    :NormalGlomPercent,
    :AbnormalPercent,
    :"AKI-Bin",
    :GGS,
    :Imploding,
    :GlomVol,
    :"% CA Glomerular",
    :"Nephron Density (All)",
    :"Nephron Density (Normal)",
    :"GV/Podocyte",
    :PodoDensity,
    :MesIndex,
    :PodoPerGlom,
    :MesVol,
    :PodoVol,
    :FIA,
    :SGS,
    :Ischemic,
    :vGlepp,
]

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
