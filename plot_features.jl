using PyPlot, Printf, Statistics

include("defs.jl")
include("annot_utils.jl")
include("clinical_utils.jl")

rm("plots", force = true, recursive = true)
mkdir("plots")

annotsx = glom_centroids(annots)

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

# Generate one figure for the given nephrectomy id.  If mode = 1,
# show all atypical glomeruli as a single class, otherwise show
# each atypical class distinctly.
function plot_one(neph_id::String, mode::Int, ixp::Int)::Int

    a = annotsx[neph_id]

    if mode == 1
        a = condense(a)
    end

    fni = parse(Int, first(split(neph_id, "_")))

    # In mode 1, skip samples with no clinical data.
    if (mode == 1) && !(fni in clin[:, :Scanner_ID])
        return ixp
    end

    PyPlot.clf()
    PyPlot.axes([0.1, 0.1, 0.7, 0.8])
    PyPlot.title(neph_id)

    # Features to be plotted with a path
    for k in boundary_types
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
            xx, yy = Float64[], Float64[]
            for u in a[k]
                push!(xx, u[1])
                push!(yy, u[2])
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
    ky = [v for v in keys(annotsx)]
    sort!(ky)
    for neph_id in ky
        ixp = plot_one(neph_id, mode, ixp)
    end

    f = [@sprintf("plots/%03d.pdf", j) for j = 0:ixp-1]
    c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=$(outname) $f`
    run(c)
end

# If mode is 1, collapse all atypical glom types into one category.
for mode in [0, 1]
    plot_all(mode, "nephrectomies$(mode).pdf")
end
