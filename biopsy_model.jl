using CSV, DataFrames, GEE, StatsModels, GLM, Printf, Statistics, PyPlot
using StatsBase

# Correct for glom density, specimen size, and extent of glom clustering

rm("plots", recursive = true, force = true)
mkdir("plots")

out = open("biopsy_model_gee.txt", "w")

function load_counts()
    cnt = []
    for offset in [false, true]
        fname = offset ? "offset" : "no_offset"
        fname = "biopsy_counts_$(fname).csv"
        cnx = open(fname) do io
            CSV.read(io, DataFrame)
        end
        cnx[!, :IDx] = string.(cnx[:, :ID])
        cnx[!, :All_total] = Float64.(cnx[:, :All_total])
        cnx[:, :offset] .= offset ? 1 : 0
        push!(cnt, cnx)
    end
    cnt = vcat(cnt...)
    cnt = sort(cnt, :IDx)
    return cnt
end

function fitmodels_getnames(gt)
    pr = Symbol("$(gt)_prop")
    pr1 = Symbol("$(gt)_prop_f1")
    pr2 = Symbol("$(gt)_prop_f2")
    pr3 = Symbol("$(gt)_prop_f3")
    nc = Symbol("$(gt)_captured")
    return nc, pr, pr1, pr2, pr3
end

function fitmodels_setup(gt, cnt)

    nc, pr, pr1, pr2, pr3 = fitmodels_getnames(gt)

    # Proportion of all gloms of the given atypical type.
    pr = Symbol("$(gt)_prop")
    tl = Symbol("$(gt)_total")
    cnt[:, pr] = cnt[:, tl] ./ cnt[:, :All_total]

    q1 = quantile(cnt[:, pr], 1 / 3)
    f1 = x -> x
    cnt[:, pr1] = f1.(cnt[:, pr])

    q2 = quantile(cnt[:, pr], 2 / 3)
    f2 = x -> sqrt(x)
    cnt[:, pr2] = f2.(cnt[:, pr])

    f3 = x -> log(0.001 + x)
    cnt[:, pr3] = f3.(cnt[:, pr])

    return cnt, f1, f2, f3
end

function fitmodels_logistic(gt, trx, cnt, fpx, ifig)

    fp1, fp2, fp3 = fpx
    nc, pr, pr1, pr2, pr3 = fitmodels_getnames(gt)

    xmax = Dict("FGGS" => 0.5, "FSGS" => 0.1, "Imploding" => 0.2, "Ischemic" => 0.3)

    ypx = []
    for tr in trx
        # Logistic regression for having >=tr gloms of a particular type.
        tnc = Symbol("T_$(gt)_captured")
        cnt[!, tnc] = [x >= tr ? 1 : 0 for x in cnt[:, nc]]
        if sum(cnt[:, tnc]) < 200
            push!(ypx, [missing, missing])
            continue
        end

        # Model without offset
        f = term(tnc) ~ sum(term.([pr1, pr2]))
        m0 = gee(f, cnt, cnt[:, :IDx], LogitLink(), BinomialVar(), IndependenceCor())
        write(out, @sprintf(">= %.0f %s glomeruli\n\n", tr, gt))
        write(out, string(m0))
        write(out, "\n\n")

        # Model with offset
        ff = Any[term(:offset)]
        for p in [pr1, pr2]
            push!(ff, term(p))
            push!(ff, term(p) & term(:offset))
        end
        f = term(tnc) ~ sum(ff)
        m1x = gee(
            f,
            cnt,
            cnt[:, :IDx],
            LogitLink(),
            BinomialVar(),
            IndependenceCor();
            dofit = false,
        )
        m1 = gee(f, cnt, cnt[:, :IDx], LogitLink(), BinomialVar(), IndependenceCor())
        write(out, @sprintf(">= %.0f %s glomeruli", tr, gt))
        write(out, string(m1))
        write(out, "\n\n")

        st = scoretest(m1x.model, m0.model)
        write(out, @sprintf("Score test p=%.4f\n\n\n", st.Pvalue))

        # Get prediction from GEE model for count.
        dp = copy(cnt[1:200, :])
        dp[1:100, pr] = range(0, xmax[gt], 100)
        dp[1:100, pr1] = fp1.(dp[1:100, pr])
        dp[1:100, pr2] = fp2.(dp[1:100, pr])
        dp[1:100, pr3] = fp3.(dp[1:100, pr])
        dp[1:100, :offset] .= 0
        dp[101:200, :] = dp[1:100, :]
        dp[101:200, :offset] .= 1
        yp = predict(m1, dp)
        push!(ypx, [dp, yp])
    end

    color = ["blue", "purple", "green", "orange", "yellow"]
    tab10 = PyPlot.cm.get_cmap("tab10")

    # Plot fitted logistic regression
    for dfp in [false, true]
        PyPlot.clf()
        PyPlot.axes([0.16, 0.1, 0.72, 0.8])
        PyPlot.grid(true)
        for i in eachindex(ypx)
            dp, yp = ypx[i]
            if ismissing(dp)
                continue
            end
            c = tab10((i - 1) / length(ypx))

            i0 = findall(dp[:, :offset] .== 0)
            i1 = findall(dp[:, :offset] .== 1)

            if dfp
                yd = yp[i0] - yp[i1]
                PyPlot.plot(
                    dp[i0, pr],
                    yd,
                    "-",
                    color = c,
                    label = @sprintf("%.0f", trx[i]),
                )
            else
                PyPlot.plot(
                    dp[i0, pr],
                    yp[i0],
                    "-",
                    color = c,
                    label = @sprintf("%.0f", trx[i]),
                )
                PyPlot.plot(dp[i1, pr], yp[i1], "--", color = c)
            end
        end
        ha, lb = PyPlot.gca().get_legend_handles_labels()
        leg = PyPlot.figlegend(ha, lb, "center right")
        leg.draw_frame(false)
        leg.set_title("T")
        PyPlot.xlim(0, xmax[gt])
        if dfp
            PyPlot.ylim(-0.3, 0.6)
            PyPlot.ylabel(
                "Difference in probabilities of T or more\n$(gt) glomeruli in biopsy",
                size = 14,
            )
        else
            PyPlot.ylim(0, 1)
            PyPlot.ylabel("Probability of T or more\n$(gt) glomeruli in biopsy", size = 14)
        end
        PyPlot.xlabel("Proportion of $(gt) glomeruli in nephrectomy", size = 14)
        PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
        ifig += 1
    end

    return ifig
end

function fitmodels_helper(gt, cnt, ifig)
    cnt, fp1, fp2, fp3 = fitmodels_setup(gt, cnt)
    ifig = fitmodels_logistic(gt, [1.0, 2, 3, 4, 5], cnt, [fp1, fp2, fp3], ifig)
    return ifig
end

function fitmodels(ifig)
    cnt = load_counts()
    for gt in ["FGGS", "FSGS", "Imploding", "Ischemic"]
        ifig = fitmodels_helper(gt, cnt, ifig)
    end
    return ifig
end

ifig = 0
ifig = fitmodels(ifig)

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=biopsy_model.pdf $f`
run(c)

close(out)
