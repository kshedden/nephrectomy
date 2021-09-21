using Statistics, PyPlot, Printf

rm("plots", force = true, recursive = true)
mkdir("plots")

include("defs.jl")

# Plot as paths
paf = ["Tissue", "Capsule", "CMJ"]

colors = Dict{String,String}(
    "All Glomeruli" => "grey",
    "FGGS" => "blue",
    "FSGS" => "magenta",
    "BSPC" => "green",
    "Ischemic" => "cyan",
    "Normal" => "black",
    "Imploding" => "red",
    "Capsule" => "purple",
    "CMJ" => "orange",
    "Tissue" => "yellow",
)


function hough(
    x::Array{Float64,1},
    y::Array{Float64,1},
    rmax::Float64,
    rsteps::Int,
    asteps::Int,
)

    hm = zeros(asteps, rsteps)
    hd = Array{Array{Int},2}(undef, asteps, rsteps)

    for j = 1:asteps
        for r = 1:rsteps
            hd[j, r] = Array{Int,1}()
        end
    end

    for j = 1:asteps
        for i in eachindex(x)
            a = 2 * pi * (j - 1) / asteps
            r = x[i] * cos(a) + y[i] * sin(a)
            if r < 0
                r = -r
                j = mod(j - 1 + div(asteps, 2), asteps) + 1
            end
            r = clamp(r / rmax, 0, 1)
            rr = Int(ceil(r * rsteps))
            hm[j, rr] = hm[j, rr] + 1

            push!(hd[j, rr], i)
        end
    end

    return hm, hd

end

for (ixp, fn) in enumerate(fi)

    println(fn)
    a = read_annot(fn)

    z = a["All Glomeruli"]

    x = [mean(u[1, :]) for u in z]
    y = [mean(u[2, :]) for u in z]

    hm, hd = hough(x, y, 200000.0, 200, 200)

    mm = findmax(hm)
    ii = hd[mm[2][1], mm[2][2]]

    PyPlot.clf()
    PyPlot.title(replace(fn, ".xml" => ""))

    PyPlot.axis("off")

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

    end

    PyPlot.plot(x, y, "o", color = "grey", mfc = "none", alpha = 0.5)
    PyPlot.plot(x[ii], y[ii], "o", mfc = "none", color = "red")

    PyPlot.savefig(@sprintf("plots/%03d.pdf", ixp))
    ixp += 1

end

f = [@sprintf("plots/%03d.pdf", j) for j = 1:length(fi)]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=lines.pdf $f`
run(c)
