using DataFrames, LinearAlgebra, Printf, StaticArrays, Statistics, PyPlot

include("defs.jl")
include("annot_utils.jl")
include("clinical_utils.jl")
include("depth_utils.jl")

outlog = open("depth_log.txt", "w")

plabel = Dict(:FGGS=>"GSG", :FSGS=>"SSG", :Imploding=>"impl.", :Ischemic=>"isch.", :Normal=>"normal",
              :All_glomeruli=>"all")
qlabel = Dict(:FGGS=>"GSG", :FSGS=>"SSG", :Imploding=>"imploding", :Ischemic=>"ischemic", :Normal=>"normal",
              :All_glomeruli=>"all")

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
        gtp = get(plabel, Symbol(gt), gt)
        gtq = get(qlabel, Symbol(gt), gt)
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
        PyPlot.ylabel(@sprintf("%s glomeruli median depth", gtq), size = 15)
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
        PyPlot.xlabel(@sprintf("(Normal + %s)/2", gtp), size = 15)
        PyPlot.ylabel(@sprintf("%s - normal median depth", gtq), size = 15)
        PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
        ifig += 1

        PyPlot.clf()
        PyPlot.title(title)
        PyPlot.hist(xx[:, 2] - xx[:, 1], ec = "black", fc = "none")
        PyPlot.ylabel("Frequency", size = 15)
        PyPlot.xlabel(@sprintf("%s depth - normal depth", gtq), size = 15)
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
            @sprintf("(%s depth - normal depth) / |%s depth + normal depth|", gtq, gtq),
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
        if gt == "Empty BC"
            continue
        end
        push!(lab, get(qlabel, Symbol(gt), gt))
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

function single_summary(out, na, dd, da)

    if !haskey(qlabel, Symbol(na))
        return
    end
    nax = qlabel[Symbol(na)]
    dp = dd ./ da
    dd = collect(skipmissing(dd))
    dp = 100 * collect(skipmissing(dp))
    @assert length(dp) == length(dd)
    write(out, @sprintf("%8.2f mean number of %s glomeruli\n", mean(dd), nax))
    write(out, @sprintf("%8.2f SD of number of %s glomeruli\n", std(dd), nax))
    write(out, @sprintf("%8.2f median number of %s glomeruli\n", median(dd), nax))
    write(out, @sprintf("%8.2f 25th percentile of number of %s glomeruli\n", quantile(dd, 0.25), nax))
    write(out, @sprintf("%8.2f 75th percentile of number of %s glomeruli\n", quantile(dd, 0.75), nax))

    if na != "All_glomeruli"
        write(out, @sprintf("%8.2f mean percentage of %s glomeruli\n", mean(dp), nax))
        write(out, @sprintf("%8.2f SD of percentage of %s glomeruli\n", std(dp), nax))
        write(out, @sprintf("%8.2f median percentage of %s glomeruli\n", median(dp), nax))
        write(out, @sprintf("%8.2f 25th percentile of percentage of %s glomeruli\n", quantile(dp, 0.25), nax))
        write(out, @sprintf("%8.2f 75th percentile of percentage of %s glomeruli\n\n", quantile(dp, 0.75), nax))
    else
        write(out, "\n")
    end
end

function depth_summary(depths)

    println(names(depths))

    out = open("depth_summary.txt", "w")

    println(names(depths))
    n0 = length(unique(depths[:, :Scanner_ID]))
    n1 = length(unique(depths[:, :TPC_ID]))
    n2 = length(unique(depths[:, :Precise_ID]))
    n3 = length(unique(depths[:, :ID]))
    write(out, @sprintf("%d scanner ID's\n", n0))
    write(out, @sprintf("%d TPC ID's\n", n1))
    write(out, @sprintf("%d Precise ID's\n", n2))
    write(out, @sprintf("%d samples\n\n", n3))

    da = depths[:, :All_glomeruli_n]
    for g in glom_types
        println(g)
        a = Symbol(@sprintf("%s_n", g))
        if string(a) in names(depths)
            single_summary(out, g, depths[:, a], da)
        end
    end

    close(out)

end

function main(annots)
    out = open("depth.txt", "w")

    ifig = 0
    pp = 0.5
    gt0 = [@sprintf("%s_q%02.0f", g, 100 * pp) for g in glom_types]
    gt = [@sprintf("%s_q%02.0f", get(plabel, Symbol(g), g), 100 * pp) for g in glom_types]

    rkeys = ["Normal", "Normal"]#, "Capsule", "CMJ"]

    # DEBUG
    dn = ["L2 depth"]#, "Spatial depth", "Tukey halfspace depth"]#, "Distance to capsule", "Distance to CMJ"]
    dpf = [l2_depth]#, spatial_depth, tukey_depth] #, boundary_depth, boundary_depth]

    for (jd, depthfun) in enumerate(dpf)

        write(outlog, @sprintf("%s\n", dn[jd]))

        depths = build_depths([0.5], depthfun, annots, glom_types)

        write(outlog, @sprintf("%d samples\n", size(depths, 1)))
        write(outlog, @sprintf("%d distinct subjects\n", length(unique(depths[:, :Scanner_ID]))))

        clin[!, :Race] = [ismissing(x) ? missing : Int(x) for x in clin[:, :Race]]
        for x in unique(clin[:, :Race])
            clin[:, Symbol("Race$(x)")] = (clin[:, :Race] .== x)
        end
        depths = leftjoin(depths, clin, on = :Scanner_ID)

        # Save all the depth data
        ff = open(@sprintf("depths_data_%s.csv", lowercase(first(split(dn[jd])))), "w")
        CSV.write(ff, depths)

        ifig = plot_depths(depths, dn[jd], ifig)

        n1 = length(unique(depths[:, :Scanner_ID]))
        n2 = length(unique(depths[:, :Precise_ID]))
        n3 = length(unique(depths[:, :TPC_ID]))
        n4 = length(unique(depths[:, :ID]))

        drslt = depth_analysis(depths, gt0)
        if jd == 1
            depth_summary(depths)
        end
        crslt = clinical_analysis(depths, gt0)
        println("crslt=", crslt)

        write(out, @sprintf("=== %s ===\n\n", dn[jd]))
        write(out, @sprintf("%d distinct samples in non-clinical depth analysis\n", n4))
        write(out, @sprintf("%d distinct scanner id's in non-clinical depth analysis\n", n1))
        write(out, @sprintf("%d distinct Precise id's in non-clinical depth analysis\n", n2))
        write(out, @sprintf("%d distinct TPC id's in non-clinical depth analysis\n\n", n3))
        write(out, "Differences in median depth based on glomerulus type:\n")
        write(out, dfclean(drslt))
        write(out, "\n\nAssociations between depth and clinical/morphometric variables:\n")
        write(out, dfclean(crslt))
        write(out, "\n\n")

        write(outlog, "\n\n")
    end

    close(out)
    return ifig
end

annotsx = glom_centroids(annots)
write(outlog, @sprintf("%d samples\n", length(annotsx)))

ifig = main(annotsx)

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=depth_plots.pdf $f`
run(c)

close(outlog)
