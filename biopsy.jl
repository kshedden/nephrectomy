using PyPlot, Statistics, Printf, DataFrames, LinearAlgebra, GeometricalPredicates, CSV

# 0.25 microns/pixel = 4 pixels/micron
# length = 1.2 cm * 10000 micron / cm * 4 pixels / micron = 48000 pixels  
# radius = 0.6 mm * 1000 micron / mm * 4 pixels / micron = 2400 pixels

include("defs.jl")

rm("plots", force = true, recursive = true)
mkdir("plots")

# Use a running PCA to identify the main axis of the capsule.
function smooth(x::Array{Float64,2}, j::Int, xw::Array{Float64,2})

	# Length of capsule path
	n = size(x, 2)  

	# Maximum number of points used for lnear approximation
	m = size(xw, 2)

	# The ideal window is j +/- s
	s = div(size(xw, 2), 2)
	@assert m == 2*s + 1

	# The window runs from i1 to i2, in case of truncation
	i1 = max(1, j-s)
	i2 = min(n, j+s)

	ii = 1
	xw .= 0
	for i in i1:i2
		xw[1, ii] = x[1, i] - x[1, j] 
		xw[2, ii] = x[2, i] - x[2, j] 
		ii += 1
	end

	xx = Symmetric(xw * xw')
	a, b = eigen(xx)
	jj = argmax(a)
	return b[:, jj]
end

# Get the four corners of a quadrilateral polygon that represents the biopsy needle
function getpoly(zz::Vector{Float64}, zt::Vector{Float64}, zd::Vector{Float64}, xy::Array{Float64,2})
	xy[1,:] = zz + 2400*zt + 100*zd
	xy[2,:] = zz + 2400*zt + 48000*zd
	xy[3,:] = zz - 2400*zt + 48000*zd
	xy[4,:] = zz - 2400*zt + 100*zd
	u = xy[2, :] - xy[1, :]
	v = xy[3, :] - xy[1, :]
	f = dot(u, v) / dot(u, u)
	xy[2, :] = xy[1, :] + f*u
	return xy
end

# Count the number of glomeruli that are inside the needle, and the total number
# of gloms in the section.
function capture(needle, gloms)

	po = Polygon(Point(needle[1,1], needle[1,2]), Point(needle[2,1], needle[2,2]), 
	             Point(needle[3,1], needle[3,2]), Point(needle[4,1], needle[4,2]))

	n = 0
	for p in eachcol(gloms)
		n += inpolygon(po, Point(p[1], p[2]))
	end

	return tuple(n, size(gloms, 2))

end

function make_plot(fn, a, ixp)

	counts = []

	# Capsule boundary
	capsule = a["Capsule"]

	# Get centroid of all glomeruli
	ctr = [0.0, 0.0]
	agl = a["All Glomeruli"]
	gloms = zeros(2, length(agl))
	for (i,g) in enumerate(agl)
		m = mean(g, dims=2)
		ctr += m
		gloms[:, i] = m
	end
	ctr = ctr ./ length(agl)

	# Get centroids of atypical gloms
	xgl = a["Atypical"]
	atp_gloms = zeros(2, length(xgl))
	for (i,g) in enumerate(xgl)
		m = mean(g, dims=2)
		atp_gloms[:, i] = m
	end

	PyPlot.clf()
	ax = PyPlot.axes()
	ax.axis("equal")
    PyPlot.axis("off")

    fni = replace(fn, ".xml" => "")
    fni = parse(Int, fni)
	PyPlot.title(fni)

	xw = zeros(2, 101)

	for (k,tp) in enumerate(a["Capsule"])

		n = size(tp, 2)
		for i in 10:n-10

			# The point where the needle enters
			zz = tp[:, i]

			# Tangent vector to the capsule
			zt = smooth(tp, i, xw)

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
			a = 2*aa*rand() - aa
			zd = (cos(a)*zn + sin(a)*zt)[:,1]

			# Allow the needle to enter occasionally
			if i%50 == 1
				xy = zeros(4, 2)
				getpoly(zz, zt, zd, xy)
				n_glom, t_glom = capture(xy, gloms)
				n_atp_glom, t_atp_glom = capture(xy, atp_gloms)

				# Require at least 20 gloms in the biopsy.
				if n_glom >= 10
					push!(counts, [n_glom, t_glom, n_atp_glom, t_atp_glom])
				end

				if i % 500 == 1			
					pa = PyPlot.matplotlib.patches.Polygon(xy, fill=true, edgecolor="grey", facecolor="lightgrey")
					ax.add_patch(pa)
				end
			end

		end

		# Plot a fragment of th capsule boundary
		PyPlot.plot(tp[1, :], tp[2, :], "-", color="yellow")

	end

	# Plot all glomeruli
	for j in 1:size(gloms, 2)
		PyPlot.plot(gloms[1, j], gloms[2, j], "o", color="grey", mfc="none")
	end

	# Plot the atypical glomeruli
	for j in 1:size(atp_gloms, 2)
		PyPlot.plot(atp_gloms[1, j], atp_gloms[2, j], "o", color="red", mfc="none")
	end

	PyPlot.savefig(@sprintf("plots/%03d.pdf", ixp))
	return tuple(counts, ixp + 1)
	
end

function do_all()

	counts, src = [], []
    ixp = 0
    for fn in fi

        println(fn)
        a = read_annot(fn)

		# Skip nephrectomies with no capsule
	    if !haskey(a, "Capsule") || length(a["Capsule"]) == 0
		    continue
	    end

	    a = condense(a)
	    counts1, ixp = make_plot(fn, a, ixp)
		push!(counts, counts1...)

		ti = replace(fn, ".xml"=>"")
		ti = parse(Int, ti)
		for _ in 1:length(counts1)
			push!(src, ti)
		end
    end

	cnt = hcat(counts...)
	cnt = DataFrame(:n_glom=>cnt[1,:], :t_glom=>cnt[2,:], :n_xglom=>cnt[3,:], :t_xglom=>cnt[4,:],
	                :src=>src)

	cnt[:, :src] = Array{Int64}(cnt[:, :src])
	cnt = sort(cnt, :src)

    return tuple(cnt, ixp)
end

cnt, ixp = do_all()

CSV.write("counts.csv", cnt)

# Scatterplot the number of glomeruli detected by the biopsy against the total
# number of glomeruli in the nephrectomy sample
PyPlot.clf()
PyPlot.grid(true)
PyPlot.plot(cnt[:, :t_glom], cnt[:, :n_glom], "o", mfc="none", alpha=0.8)
PyPlot.ylabel("Biopsied glomeruli", size=15)
PyPlot.xlabel("Total glomeruli", size=15)
PyPlot.savefig("all_gloms.pdf")

# Same as above, except using only atypical glmoeruli
PyPlot.clf()
PyPlot.grid(true)
PyPlot.plot(cnt[:, :n_xglom], cnt[:, :t_xglom], "o", mfc="none", alpha=0.8)
PyPlot.ylabel("Biopsied atypical glomeruli", size=15)
PyPlot.xlabel("Total atypical glomeruli", size=15)
PyPlot.savefig("atypical_gloms.pdf")


f = [@sprintf("plots/%03d.pdf", j) for j = 0:ixp-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=biopsy.pdf $f`
run(c)
