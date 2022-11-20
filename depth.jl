using DataFrames, LinearAlgebra, Printf, StaticArrays, Statistics, PyPlot

# TODO get better matching data
# TODO plot showing difference of means between each glom type and normal

include("defs.jl")
include("annot_utils.jl")
include("clinical_utils.jl")
include("depth_utils.jl")

rm("plots", recursive = true, force = true)
mkdir("plots")

function dfclean(df)
    x = split(string(df), "\n")
    x = vcat(x[2], x[4:end])
    x = join(x, "\n")
    return x
end

function plot_depths(depths, title, ifig)

    zz = Dict()
    for gt in glom_types

        if gt == "Normal" || gt == "All_glomeruli"
            continue
        end
        pp = 0.5
        gtx = @sprintf("%s_q%02.0f", gt, 100 * pp)
        if !(gtx in names(depths))
            continue
        end

        xx = depths[:, ["Normal_q50", gtx]]
        xx = xx[completecases(xx), :]
        xx = sqrt.(xx)
        zz[gt] = xx

        # Scatterplot of atypical depth on typical depth
        PyPlot.clf()
        PyPlot.grid(true)
        PyPlot.title(title)
        mn1 = minimum(xx[:, 1])
        mn2 = minimum(xx[:, 2])
        mx1 = maximum(xx[:, 1])
        mx2 = maximum(xx[:, 2])
        mn = min(mn1, mn2)
        mx = max(mx1, mx2)
        PyPlot.axline([0.9 * mn, 0.9 * mn], [1.1 * mx, 1.1 * mx])
        PyPlot.plot(xx[:, 1], xx[:, 2], "o", mec = "blue", mfc = "none")
        PyPlot.xlim(0.9 * mn, 1.1 * mx)
        PyPlot.ylim(0.9 * mn, 1.1 * mx)
        PyPlot.xlabel("Normal glomeruli median depth", size = 15)
        PyPlot.ylabel(@sprintf("%s glomeruli median depth", gt), size = 15)
        PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
        ifig += 1

        # Mean/difference plot
        PyPlot.clf()
        PyPlot.grid(true)
        PyPlot.title(title)
        mn1 = minimum(xx[:, 1])
        mn2 = minimum(xx[:, 2])
        mx1 = maximum(xx[:, 1])
        mx2 = maximum(xx[:, 2])
        mn = min(mn1, mn2)
        mx = max(mx1, mx2)
        PyPlot.axhline(0, lw = 4, color = "grey")
        PyPlot.plot(
            (xx[:, 1] + xx[:, 2]) / 2,
            xx[:, 2] - xx[:, 1],
            "o",
            mec = "blue",
            mfc = "none",
        )
        PyPlot.xlabel(@sprintf("(Normal + %s)/2", gt), size = 15)
        PyPlot.ylabel(@sprintf("%s - normal median depth", gt), size = 15)
        PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
        ifig += 1

        PyPlot.clf()
        PyPlot.title(title)
        PyPlot.hist(xx[:, 2] - xx[:, 1], ec = "black", fc = "none")
        PyPlot.ylabel("Frequency", size = 15)
        PyPlot.xlabel(@sprintf("%s depth - normal depth", gt), size = 15)
        PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
        ifig += 1

        PyPlot.clf()
        PyPlot.title(title)
        PyPlot.hist(
            (xx[:, 2] - xx[:, 1]) ./ abs.(xx[:, 1] + xx[:, 2]),
            ec = "black",
            fc = "none",
        )
        PyPlot.ylabel("Frequency", size = 15)
        PyPlot.xlabel(
            @sprintf("(%s depth - normal depth) / |%s depth + normal depth|", gt, gt),
            size = 13,
        )
        PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
        ifig += 1
    end

    lab, val = [], []
    for gt in glom_types
        if !haskey(zz, gt)
            continue
        end
        push!(lab, gt)
        xx = zz[gt]
        push!(val, xx[:, 2] - xx[:, 1])
    end
    PyPlot.clf()
    PyPlot.title(title)
    PyPlot.axhline(y = 0, color = "grey")
    PyPlot.boxplot(val, labels = lab)
    PyPlot.ylabel("Depth relative to normal", size = 15)
    PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
    ifig += 1

    vm = mean.(val)
    vs = std.(val) ./ sqrt.(length.(val))
    ii = range(1, length(val))
    PyPlot.clf()
    PyPlot.title(title)
    PyPlot.grid(true, axis = "y")
    for i in eachindex(ii)
        PyPlot.plot([ii[i], ii[i]], [vm[i] - 2 * vs[i], vm[i] + 2 * vs[i]], color = "grey")
    end
    PyPlot.plot(ii, vm, "o", color = "black")
    PyPlot.gca().set_xticks(ii)
    PyPlot.gca().set_xticklabels(lab)
    PyPlot.ylabel("Depth relative to normal", size = 15)
    PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
    ifig += 1

    return ifig
end

function main(annots)
    out = open("depth.txt", "w")

    ifig = 0
    pp = 0.5
    gt = [@sprintf("%s_q%02.0f", g, 100 * pp) for g in glom_types]

    rkeys = ["Normal", "Normal"]#, "Capsule", "CMJ"]

    dn = ["L2 depth", "Spatial depth", "Tukey halfspace depth"]#, "Distance to capsule", "Distance to CMJ"]
    dpf = [l2_depth, spatial_depth, tukey_depth] #, boundary_depth, boundary_depth]

    for (jd, depthfun) in enumerate(dpf)

        depths = build_depths([0.5], depthfun, annots, glom_types)

        clin[!, :Race] = [ismissing(x) ? missing : Int(x) for x in clin[:, :Race]]
        for x in unique(clin[:, :Race])
            clin[:, Symbol("Race$(x)")] = (clin[:, :Race] .== x)
        end
        depths = leftjoin(depths, clin, on = :Scanner_ID)

        if jd == 1
            CSV.write("depths_data.csv", depths)
        end

        ifig = plot_depths(depths, dn[jd], ifig)

        drslt = depth_analysis(depths, gt)
        crslt = clinical_analysis(depths, gt)
        write(out, @sprintf("=== %s ===\n\n", dn[jd]))
        write(out, "Differences in median depth based on glomerulus type:\n")
        write(out, dfclean(drslt))
        write(out, "\n\nAssociations between depth and clinical/morphometric variables:\n")
        write(out, dfclean(crslt))
        write(out, "\n\n")
    end

    close(out)
    return ifig
end

annotsx = glom_centroids(annots)

ifig = main(annotsx)

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=depth_plots.pdf $f`
run(c)
