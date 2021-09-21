using PyPlot, Statistics, Printf, DataFrames, IterTools

include("defs.jl")
include("pair_corr_utils.jl")

rm("plots", force = true, recursive = true)
mkdir("plots")

ptq = ["All Glomeruli", "FGGS", "Ischemic", "FSGS", "Imploding"]

di = Dict{Tuple{String,String},Array{Float64,1}}()

# Count the number of samples with at least 1/3/5 glomeruli of a given type.
smry = Dict{String,Array{Int,1}}()

# Loop over all samples
for fn in fi

    println(fn)
    a = read_annot(fn)

    for q in ptq
        if !haskey(smry, q)
            smry[q] = [0, 0, 0]
        end
        if haskey(a, q)
            m = length(a[q])
            smry[q] = smry[q] + [Int(m >= 1), Int(m >= 3), Int(m >= 5)]
        end
    end

    dit = pair_corr(a, use = ptq)

    # Normalize the distances to the all vs. all values.
    vr = dit[tuple("All Glomeruli", "All Glomeruli")]
    for (k, v) in dit
        if !haskey(di, k)
            di[k] = []
        end
        if length(v) > 0
            u = v ./ vr
            push!(di[k], u...)
        end
    end

    for (k, v) in di
        di[k] = sort(v)
    end

end

# If mode is "same", plot results for comparing within a group, if mode is
# "diff", plot results for comparing between groups, if mode is "all", plot
# all results.
function plot_cumulative(ixp, mode; xmax = 0, ymax = 16)

    PyPlot.clf()
    PyPlot.figure(figsize = (8, 5))
    PyPlot.axes([0.1, 0.1, 0.64, 0.8])
    PyPlot.grid(true)

    for (k, v) in di

        if (k[1] == "All Glomeruli") || (k[2] == "All Glomeruli")
            continue
        end

        if mode == "same"
            if k[1] != k[2]
                continue
            end
        elseif mode == "diff"
            if k[1] == k[2]
                continue
            end
        end

        n = length(v)
        pp = range(1 / n, 1 - 1 / n, length = n)
        lab = @sprintf("%s/%s", k[1], k[2])
        PyPlot.plot(pp, log.(v) / log(2), alpha = 0.5, label = lab)
    end

    PyPlot.xlabel("Probability", size = 15)
    PyPlot.ylabel("Log relative distance", size = 15)

    ha, lb = PyPlot.gca().get_legend_handles_labels()
    leg = PyPlot.figlegend(ha, lb, "center right")
    leg.draw_frame(false)

    PyPlot.savefig(@sprintf("plots/%03d.pdf", ixp))
    return ixp + 1

end

function genplots()
    ixp = 0
    for mode in ["same", "diff", "all"]
        ixp = plot_cumulative(ixp, mode, xmax = 0, ymax = 16)
        ixp = plot_cumulative(ixp, mode, xmax = -6, ymax = 11.5)
    end
    return ixp
end

ixp = genplots()

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ixp-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=pair_corr.pdf $f`
run(c)

ss = DataFrame(smry)
ss[!, "Samples"] = [">=1", ">=3", ">=5"]
open("summary.txt", "w") do io
    write(io, string(ss))
end
