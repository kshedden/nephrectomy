using PyPlot, Printf, Statistics

include("defs.jl")
include("annot_utils.jl")
include("clinical_utils.jl")

rm("plots", force = true, recursive = true)
mkdir("plots")

# Plot as points
ptf = [
    "All_glomeruli",
    "Normal",
    "FGGS",
    "Ischemic",
    "FSGS",
    "BSPC",
    "Imploding",
    "Atypical",
    "Empty BC",
]

# Plot as paths
paf = ["Tissue", "Capsule", "CMJ", "Cortex"]

# Make a plot showing the disconneced components of
# each nephrectomy.
function plot_components(neph_id::String, ixp::Int)::Int

    a = annots[neph_id]

    fni = parse(Int, neph_id)

    PyPlot.clf()
    PyPlot.axes([0.1, 0.1, 0.7, 0.8])
    PyPlot.title("$(fni)")

    if haskey(a, "Cortex")
        for (j, u) in enumerate(a["Cortex"])
            PyPlot.plot(u[1, :], u[2, :], "-"; color = "grey", lw = 3)
        end
    end

    colors = Dict(
        nothing => "red",
        1 => "blue",
        2 => "green",
        3 => "orange",
        4 => "yellow",
        5 => "cyan",
    )

    if haskey(a, "All_glomeruli_components")
        cmp = a["All_glomeruli_components"]
        uc = unique(cmp)
        for (j, u) in enumerate(uc)
            xx, yy = Float64[], Float64[]
            for (i, g) in enumerate(a["All_glomeruli"])

                if cmp[i] != u
                    continue
                end

                push!(xx, mean(g[1, :]))
                push!(yy, mean(g[2, :]))
            end
            PyPlot.plot(
                xx,
                yy,
                "o",
                mfc = "none",
                color = colors[u],
                alpha = 0.5,
                zorder = 2,
            )
        end
    end

    PyPlot.axis("off")

    PyPlot.savefig(@sprintf("plots/%03d.pdf", ixp))
    return ixp + 1
end

# Generate one figure for the given nephrectomy id.  If mode = 1,
# show all atypical glomeruli as a single class, otherwise show
# each atypical class distinctly.
function plot_one(neph_id::String, mode::Int, ixp::Int)::Int

    a = annots[neph_id]

    if mode == 1
        a = condense(a)
    end

    fni = parse(Int, neph_id)

    # In mode 1, skip samples with no clinical data.
    if (mode == 1) && !(fni in clin[:, :Scanner_ID])
        return ixp
    end

    PyPlot.clf()
    PyPlot.axes([0.1, 0.1, 0.7, 0.8])
    PyPlot.title("$(fni)")

    # Features to be plotted with a path
    for k in paf
        if haskey(a, k)
            v = a[k]
            for (j, u) in enumerate(v)
                args = (color = colors[k], alpha = 0.5, zorder = 1)
                if j == 1
                    args = merge(args, (label = k,))
                end
                PyPlot.plot(u[1, :], u[2, :], "-"; args...)
            end
        end
    end

    # Features to be plotted with a point
    for k in ptf
        if k == "All_glomeruli"
            continue
        end
        if haskey(a, k)
            v = a[k]
            xx, yy = Float64[], Float64[]
            for u in v
                x = mean(u[1, :])
                y = mean(u[2, :])
                push!(xx, x)
                push!(yy, y)
            end

            PyPlot.plot(
                xx,
                yy,
                "o",
                mfc = "none",
                label = k,
                color = colors[k],
                alpha = 0.5,
                zorder = k == "Normal" ? 2 : 3,
            )
        end
    end

    PyPlot.axis("off")

    # If available, print some clinical information at the bottom
    # of the graph.
    if fni in clin[:, :Scanner_ID]
        ii = findfirst(clin[:, :Scanner_ID] .== fni)
        age = clin[ii, :Age]
        bmi = ismissing(clin[ii, :BMI]) ? "" : @sprintf("%.1f", clin[ii, :BMI])
        sex = clin[ii, :Sex_F_M]
        htn = clin[ii, :Hypertension_No_Yes]
        dm = clin[ii, :Diabetes_No_Yes]
        race = clin[ii, :Race]
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

function plot_components()

    rm("plots", force = true, recursive = true)
    mkdir("plots")

    ixp = 0
    for neph_id in keys(annots)
        ixp = plot_components(neph_id, ixp)
    end

    f = [@sprintf("plots/%03d.pdf", j) for j = 0:ixp-1]
    c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=components.pdf $f`
    run(c)
end

# Plot the tissue islands bound by closed loops of cortex.
plot_components()

# If mode is 1, collapse all atypical glom types into one category.
level = 1
for mode in [0, 1]
    plot_sorted(mode, level, "nephrectomies$(level)_sorted.pdf")
    plot_all(mode, "nephrectomies$(mode).pdf")
end
