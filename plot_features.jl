using PyPlot, Printf, Statistics

include("defs.jl")
include("annot_utils.jl")
include("clinical_utils.jl")

rm("plots", force = true, recursive = true)
mkdir("plots")

# Plot as points
ptf = ["All_glomeruli", "FGGS", "Ischemic", "FSGS", "BSPC", "Imploding", "Atypical"]

# Plot as paths
paf = ["Tissue", "Capsule", "CMJ", "Cortex"]

function plot_one(neph_id::String, mode::Int, ixp::Int)::Int

    a = annots[neph_id]

    if mode == 1
        a = condense(a)
    end

    fni = parse(Int, neph_id)
    if (mode == 1) && !haskey(sid_rownum, fni)
        return ixp
    end

    PyPlot.clf()
    PyPlot.axes([0.1, 0.1, 0.7, 0.8])
    PyPlot.title("$(fni)")

    for (k, v) in a

        # Features to be plotted with a path
        if k in paf
            for (j, u) in enumerate(v)
                args = (color = colors[k], alpha = 0.5, zorder = 1)
                if j == 1
                    args = merge(args, (label = k,))
                end
                PyPlot.plot(u[1, :], u[2, :], "-"; args...)
            end
        end

        # Features to be plotted with a point
        if k in ptf

            xx, yy = [], []
            for (j,u) in enumerate(v)
                for i in size(u, 2)
                    x = mean(u[1, :])
                    y = mean(u[2, :])
                    push!(xx, x)
                    push!(yy, y)
                end
            end
            PyPlot.plot(
                xx,
                yy,
                "o",
                mfc = "none",
                label = k,
                color = colors[k],
                alpha = 0.5,
                zorder = 2,
            )
        end
    end

    PyPlot.axis("off")

    # Annotations
    if haskey(sid_rownum, fni)
        ii = sid_rownum[fni]
        age = df[ii, :age]
        bmi = df[ii, :bmi]
        sex = df[ii, :SEX]
        htn = df[ii, :htn]
        dm = df[ii, :dm]
        race = df[ii, :RACE]
        la = "age=$(age), sex=$(sex), HTN=$(htn), BMI=$(bmi), DM=$(dm), race=$(race)"
        PyPlot.figtext(0.05, 0.02, la)
    end

    ha, lb = PyPlot.gca().get_legend_handles_labels()
    leg = PyPlot.figlegend(ha, lb, loc = "center right")
    leg.draw_frame(false)

    PyPlot.savefig(@sprintf("plots/%03d.pdf", ixp))
    return ixp + 1

end

function plot_all(mode::Int, outname::String)

    rm("plots", force = true, recursive = true)
    mkdir("plots")

    ixp = 0
    for neph_id in keys(annots)
        ixp = plot_one(neph_id, mode, ixp)
    end

    f = [@sprintf("plots/%03d.pdf", j) for j = 0:ixp-1]
    c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=$(outname) $f`
    run(c)
end

function plot_sorted(mode::Int, level, outname::String)

    u_df = open("age_scores.csv") do io
        CSV.read(io, DataFrame)
    end

    s = Symbol("score$(level)")
    ii = sortperm(u_df[:, s])
    idx = u_df[ii, :Scanner_id]
    fi = ["$(id).xml" for id in idx]
    plot_all(mode, outname)
end

# If mode is 1, collapse all atypical glom types into one category.
level = 1
for mode in [0, 1]
    plot_sorted(mode, level, "nephrectomies_sorted_$(level).pdf")
    plot_all(mode, "nephrectomies$(mode).pdf")
end
