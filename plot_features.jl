using PyPlot, Printf, Statistics

include("defs.jl")
include("clinical_utils.jl")

rm("plots", force = true, recursive = true)
mkdir("plots")

# Plot as points
ptf = ["All Glomeruli", "FGGS", "Ischemic", "FSGS", "BSPC", "Imploding", "Atypical"]

# Plot as paths
paf = ["Tissue", "Capsule", "CMJ"]

# If mode is 1, collapse all atypical glom types into one category.
mode = 0

function plot_one(fn::String, mode::Int, ixp::Int)::Int

    a = read_annot(fn)

    if mode == 1
        a = condense(a)
    end

    fni = replace(fn, ".xml" => "")
    fni = parse(Int, fni)

    if (mode == 1) && !haskey(sid_rownum, fni)
        return
    end

    PyPlot.clf()
    PyPlot.axes([0.1, 0.1, 0.7, 0.8])
    PyPlot.title("$(fni)")

    for (k, v) in a

        # Features to be plotted with a path
        if k in paf
            for (j, u) in enumerate(v)
                if j == 1
                    PyPlot.plot(
                        u[1, :],
                        u[2, :],
                        "-",
                        label = k,
                        color = colors[k],
                        alpha = 0.5,
                        zorder = 1,
                    )
                else
                    PyPlot.plot(
                        u[1, :],
                        u[2, :],
                        "-",
                        color = colors[k],
                        alpha = 0.5,
                        zorder = 1,
                    )
                end
            end
        end

        # Features to be plotted with a point
        if k in ptf
            xx, yy = [], []
            for u in v
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


function plot_all(mode::Int, fi, outname::String)

    rm("plots", force = true, recursive = true)
    mkdir("plots")

    ixp = 0
    for fn in fi
        ixp = plot_one(fn, mode, ixp)
    end

    f = [@sprintf("plots/%03d.pdf", j) for j = 1:ixp-1]
    c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=$(outname) $f`
    run(c)
end

function plot_sorted(mode::Int, level)

	u_df = open("age_scores.csv") do io
		CSV.read(io, DataFrame)
	end

	s = Symbol("score$(level)")
	ii = sortperm(u_df[:, s])
	idx = u_df[ii, :Scanner_id]
	fi = ["$(id).xml" for id in idx]
	plot_all(mode, fi, "nephrectomies_sorted_$(level).pdf")

end

plot_sorted(mode, 1)
#plot_all(mode, fi, "nephrectomies$(mode).pdf")
