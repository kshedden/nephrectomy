using PyPlot, Statistics, Printf, DataFrames, LinearAlgebra, CSV
using UnicodePlots
import GeometryBasics
import PolygonOps

# 0.25 microns/pixel = 4 pixels/micron
# length = 1.2 cm * 10000 micron / cm * 4 pixels / micron = 48000 pixels
# radius = 0.6 mm * 1000 micron / mm * 4 pixels / micron = 2400 pixels
# offset = 0.25 cm * 10000 micron / cm * 4 pixels / micron = 10000 pixels

include("defs.jl")
include("annot_utils.jl")

rm("plots", force = true, recursive = true)
mkdir("plots")

# Use a running PCA to identify a tangent vector to the capsule.
function tangent(x::Matrix{Float64}, j::Int, xw::Matrix{Float64})

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
    xy = zeros(2, 4)
    xy[:, 1] = zz + 2400 * zt + 100 * zd   # ll
    xy[:, 2] = zz + 2400 * zt + 48000 * zd # ul
    xy[:, 3] = zz - 2400 * zt + 48000 * zd # ur
    xy[:, 4] = zz - 2400 * zt + 100 * zd   # lr
    u = xy[:, 2] - xy[:, 1]
    v = xy[:, 3] - xy[:, 1]
    f = dot(u, v) / dot(u, u)
    xy[:, 2] = xy[:, 1] + f * u

    return xy
end

# Count the number of glomeruli that are inside the needle, and the total number
# of gloms in the section.
function capture(needle, gloms)

    po = [
        SVector(needle[1, 1], needle[2, 1]),
        SVector(needle[1, 2], needle[2, 2]),
        SVector(needle[1, 3], needle[2, 3]),
        SVector(needle[1, 4], needle[2, 4]),
        SVector(needle[1, 1], needle[2, 1]),
    ]

    n = 0
    for g in gloms
        n += PolygonOps.inpolygon(SVector(g[1], g[2]), po)
    end

    return n, length(gloms)
end


function tissue_frac(tis, xy)

    Gpoint = GeometryBasics.Point
    Gline = GeometryBasics.Line

    @assert length(tis) == 1
    tix = first(tis)

    # Two long sides of the needle
    f1 = Gline(Gpoint(xy[1, 1], xy[2, 1]), Gpoint(xy[1, 2], xy[2, 2]))
    f3 = Gline(Gpoint(xy[1, 3], xy[2, 3]), Gpoint(xy[1, 4], xy[2, 4]))

    frac = 1.0
    for i = 1:size(tix, 2)-1

        # One line segement of the tissue
        lt = Gline(Gpoint(tix[1, i], tix[2, i]), Gpoint(tix[1, i+1], tix[2, i+1]))

        i1, p1 = GeometryBasics.intersects(f1, lt)
        i3, p3 = GeometryBasics.intersects(f3, lt)
        if i1 || i3
            pp = i1 ? p1 : p3
            # Distance to inner left
            h1 = norm(pp - Gpoint(xy[1, 1], xy[2, 1]))
            # Distance to upper left
            h2 = norm(pp - Gpoint(xy[1, 2], xy[2, 2]))
            frac = min(frac, h1 / (h1 + h2))
        end
    end
    return frac
end

function make_plot(id, neph, offset, ixp)

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

    # Plot the tissue boundary
    tis = get(neph, "Tissue", [])
    for ti in tis
        PyPlot.plot(ti[1, :], ti[2, :], "-", color = "black")
    end

    # Loop over pieces of capsule boundary
    for (k, caps) in enumerate(neph["Capsule"])

        # Number of points on the capsule boundary
        n = size(caps, 2)
        for i = 10:n-10

            # The point where the needle enters
            zz = caps[:, i]

            # Tangent vector to the capsule
            zt = tangent(caps, i, xw)

            # Get the local glomeruli to the needle entry point
            # This is used to determine which side of the capsule
            # is the kidney
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

            # Push the needle inward so that it starts below the surface.
            zz1 = zz + offset * zd

            # Allow the needle to enter occasionally
            if i % 10 == 1
                xy0 = getpoly(zz, zt, zd)
                xy = getpoly(zz1, zt, zd)
                n_glom, t_glom = capture(xy, neph["All_glomeruli"])
                fr = tissue_frac(tis, xy0)
                fr0 = tissue_frac(tis, xy)
                if min(fr, fr0) < 0.75
                    continue
                end

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
                if i % 200 == 1
                    pa = PyPlot.matplotlib.patches.Polygon(
                        xy',
                        fill = true,
                        edgecolor = "grey",
                        facecolor = "lightgrey",
                    )
                    ax.add_patch(pa)
                end
            end
        end

        # Plot a fragment of the capsule boundary
        PyPlot.plot(caps[1, :], caps[2, :], "-", color = "yellow")
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

function do_all(annots, offset)

    idx = [k for k in keys(annots)]
    sort!(idx)

    counts, src = [], []
    ixp = 0
    for id in idx
        println(id)
        a = annots[id]

        # Skip nephrectomies with no capsule
        if !haskey(a, "Capsule") || length(a["Capsule"]) == 0
            continue
        end

        counts1, ixp = make_plot(id, a, offset, ixp)
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
    cnt[:, :ID] = src
    cnt = sort(cnt, :ID)

    return cnt, ixp
end

annotsx = glom_centroids(annots)
for offset in [0, 10000]

    rm("plots", force = true, recursive = true)
    mkdir("plots")

    cnt, ixp = do_all(annotsx, offset)

    s = offset > 0 ? "offset" : "no_offset"
    CSV.write(@sprintf("biopsy_counts_%s.csv", s), cnt)

    f = [@sprintf("plots/%03d.pdf", j) for j = 0:ixp-1]
    c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=biopsy_$s.pdf $f`
    run(c)
end
