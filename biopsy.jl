using PyPlot, Statistics, Printf, DataFrames, LinearAlgebra, GeometricalPredicates, CSV

# TODO: add counts for specific subtypes
# rename columns to more informative names

# 0.25 microns/pixel = 4 pixels/micron
# length = 1.2 cm * 10000 micron / cm * 4 pixels / micron = 48000 pixels
# radius = 0.6 mm * 1000 micron / mm * 4 pixels / micron = 2400 pixels

include("defs.jl")
include("annot_utils.jl")

rm("plots", force = true, recursive = true)
mkdir("plots")

# Use a running PCA to identify a tangent vector to the capsule.
function smooth(x::Matrix{Float64}, j::Int, xw::Matrix{Float64})

    # Length of capsule path
    n = size(x, 2)

    # Maximum number of points used for linear approximation
    m = size(xw, 2)

    # The ideal window is j +/- s
    s = div(size(xw, 2), 2)
    @assert m == 2 * s + 1

    # The window runs from i1 to i2, in case of truncation
    i1 = max(1, j - s)
    i2 = min(n, j + s)

    ii = 1
    xw .= 0
    for i = i1:i2
        xw[:, ii] = x[:, i] - x[:, j]
        ii += 1
    end

    xx = Symmetric(xw * xw')
    a, b = eigen(xx)
    jj = argmax(a)
    return b[:, jj]
end

# Get the four corners of a quadrilateral polygon that represents the biopsy needle
function getpoly(zz, zt, zd)
    xy = zeros(4, 2)
    xy[1, :] = zz + 2400 * zt + 100 * zd
    xy[2, :] = zz + 2400 * zt + 48000 * zd
    xy[3, :] = zz - 2400 * zt + 48000 * zd
    xy[4, :] = zz - 2400 * zt + 100 * zd
    u = xy[2, :] - xy[1, :]
    v = xy[3, :] - xy[1, :]
    f = dot(u, v) / dot(u, u)
    xy[2, :] = xy[1, :] + f * u
    return xy
end

# Count the number of glomeruli that are inside the needle, and the total number
# of gloms in the section.
function capture(needle, gloms)

    po = Polygon(
        Point(needle[1, 1], needle[1, 2]),
        Point(needle[2, 1], needle[2, 2]),
        Point(needle[3, 1], needle[3, 2]),
        Point(needle[4, 1], needle[4, 2]),
    )

    n = 0
    for g in gloms
        n += inpolygon(po, Point(g[1], g[2]))
    end

    return n, length(gloms)
end

function make_plot(id, neph, ixp)

    counts = []

    # Capsule boundary
    capsule = neph["Capsule"]

    # Get centroid of all glomeruli
    ctr = sum(neph["All_glomeruli"]) / length(neph["All_glomeruli"])

    PyPlot.clf()
    ax = PyPlot.axes()
    ax.axis("equal")
    PyPlot.axis("off")

    PyPlot.title(id)

    xw = zeros(2, 101)

    # Loop over pieces of capsule boundary
    for (k, tp) in enumerate(neph["Capsule"])

        n = size(tp, 2)
        for i = 10:n-10

            # The point where the needle enters
            zz = tp[:, i]

            # Tangent vector to the capsule
            zt = smooth(tp, i, xw)

            dd = [norm(zz - x) for x in neph["All_glomeruli"]]
            ii = findall(dd .< quantile(dd, 0.2))
            ctr = sum(neph["All_glomeruli"][ii]) / length(ii)

            # Direction vector toward center
            dc = ctr - zz
            dc = dc ./ norm(dc)

            # Normal vector to the capsule
            zn = dc - dot(dc, zt) * zt
            zn = zn ./ norm(zn)
            if dot(zn, dc) < 0
                zn .= -zn
            end
            zn = zn ./ norm(zn)

            # Let the needle enter at a random angle centered on orthogonal
            aa = pi / 5
            a = 2 * aa * rand() - aa
            zd = cos(a) * zn + sin(a) * zt

            # Allow the needle to enter occasionally
            if i % 50 == 1
                xy = getpoly(zz, zt, zd)
                n_glom, t_glom = capture(xy, neph["All_glomeruli"])

                c = [n_glom, t_glom]
                for atp in atypical_glom_types
                    n_atp_glom, t_atp_glom = 0, 0
                    if haskey(neph, atp)
                        n_atp_glom, t_atp_glom = capture(xy, neph[atp])
                    end
                    push!(c, n_atp_glom, t_atp_glom)
                end
                push!(counts, c)

                # Draw the needle
                if i % 500 == 1
                    pa = PyPlot.matplotlib.patches.Polygon(
                        xy,
                        fill = true,
                        edgecolor = "grey",
                        facecolor = "lightgrey",
                    )
                    ax.add_patch(pa)
                end
            end
        end

        # Plot a fragment of the capsule boundary
        PyPlot.plot(tp[1, :], tp[2, :], "-", color = "yellow")
    end

    # Plot all glomeruli
    for g in neph["All_glomeruli"]
        PyPlot.plot(g[1], g[2], "o", color = "grey", mfc = "none")
    end

    # Plot the atypical glomeruli
    for atp in atypical_glom_types
        if haskey(neph, atp)
            x = [u[1] for u in neph[atp]]
            y = [u[2] for u in neph[atp]]
            PyPlot.plot(x, y, "o", label = atp, mfc = "none")
        end
    end

    ha, lb = PyPlot.gca().get_legend_handles_labels()
    leg = PyPlot.figlegend(ha, lb)
    leg.draw_frame(false)
    PyPlot.savefig(@sprintf("plots/%03d.pdf", ixp))
    return counts, ixp + 1
end

function do_all(annots)

    counts, src = [], []
    ixp = 0
    for (id, a) in annots
        println(id)

        # Skip nephrectomies with no capsule
        if !haskey(a, "Capsule") || length(a["Capsule"]) == 0
            continue
        end

        counts1, ixp = make_plot(id, a, ixp)
        push!(counts, counts1...)

        for _ in eachindex(counts1)
            push!(src, id)
        end
    end

    cnt = hcat(counts...)
    na = String["All_captured", "All_total"]
    for u in atypical_glom_types
        push!(na, "$(u)_captured")
        push!(na, "$(u)_total")
    end
    cnt = DataFrame(cnt', na)
    cnt[:, :ID] = parse.(Int, src)
    cnt = sort(cnt, :ID)

    return cnt, ixp
end

annotsx = glom_centroids(annots)
annotsx = major_components(annotsx)
cnt, ixp = do_all(annotsx)

CSV.write("biopsy_counts.csv", cnt)

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ixp-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=biopsy.pdf $f`
run(c)

ifig = 0
rm("plots", force = true, recursive = true)
mkdir("plots")

function count_scatterplots(ifig)

    # Scatterplot the number of glomeruli detected by the biopsy against the total
    # number of glomeruli in the nephrectomy sample
    PyPlot.clf()
    PyPlot.grid(true)
    PyPlot.plot(cnt[:, :All_total], cnt[:, :All_captured], "o", mfc = "none", alpha = 0.8)
    PyPlot.ylabel("Biopsied glomeruli", size = 15)
    PyPlot.xlabel("Total glomeruli", size = 15)
    PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
    ifig += 1

    # Same as above, except using only atypical glmoeruli
    for gt in atypical_glom_types
        PyPlot.clf()
        PyPlot.grid(true)
        PyPlot.plot(
            cnt[:, "$(gt)_total"],
            cnt[:, "$(gt)_captured"],
            "o",
            mfc = "none",
            alpha = 0.8,
        )
        PyPlot.ylabel("Biopsied $(gt) glomeruli", size = 15)
        PyPlot.xlabel("Total $(gt) glomeruli", size = 15)
        PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
        ifig += 1
    end

    return ifig
end

ifig = count_scatterplots(ifig)

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=biopsy_scatterplots.pdf $f`
run(c)
