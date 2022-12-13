using CSV, DataFrames, GEE, StatsModels, GLM, Printf, Statistics, PyPlot
using StatsBase

rm("plots", recursive = true, force = true)
mkdir("plots")

xmax = Dict("FGGS"=>0.5, "FSGS"=>0.1, "Imploding"=>0.2, "Ischemic"=>0.3, "Normal"=>1.0)
color = ["blue", "purple", "green", "orange", "yellow"]
tab10 = PyPlot.cm.get_cmap("tab10")

plabel = Dict(:FGGS=>"GSG", :FSGS=>"SSG", :Imploding=>"impl.", :Ischemic=>"isch.", :Normal=>"normal")
qlabel = Dict(:FGGS=>"GSG", :FSGS=>"SSG", :Imploding=>"imploding", :Ischemic=>"ischemic", :Normal=>"normal")

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

    cnt1 = filter(r->r.offset > 0, cnt)
    id = unique(cnt1[:, :ID])
    cnt = filter(r->r.ID in id, cnt)

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

    f1 = x -> x
    cnt[:, pr1] = f1.(cnt[:, pr])

    f2 = x -> x^2
    cnt[:, pr2] = f2.(cnt[:, pr])

    f3 = x -> sqrt(x)
    cnt[:, pr3] = f3.(cnt[:, pr])

    return cnt, f1, f2, f3
end

function local_logistic(y, x, xp; bw=0.1)

    yp = zeros(length(xp))

    n = length(x)
    xm = ones(n, 2)
    xm[:, 2] = x

    for (i, xx) in enumerate(xp)
        w = (x .- xx) ./ bw
        w = exp.(-w.^2/2)
        w ./= sum(w)

        m = glm(xm, y, Binomial(); wts=w)
        yp[i] = predict(m, [1 xx;])[1]
    end

    return yp
end

function fitmodels_local(gt, trx, cnt, ifig)

    nc, pr, pr1, pr2, pr3 = fitmodels_getnames(gt)
    xp = range(0, xmax[gt], 50)
    bw = 0.2

    ypx = []
    for tr in trx
        # Local regression for having >=tr gloms of a particular type.
        tnc = Symbol("T_$(gt)_captured")
        cnt[!, tnc] = [x >= tr ? 1 : 0 for x in cnt[:, nc]]

        cnt1 = filter(r->r.offset > 0, cnt)
        cnt0 = filter(r->r.offset == 0, cnt)

        yh = []
        xx = []
        offset = []
        for (jj, cnx) in enumerate([cnt0, cnt1])
            y = cnx[:, tnc]
            x = cnx[:, pr]
            yhat = local_logistic(y, x, xp; bw=bw)
            push!(yh, yhat)
            push!(xx, xp)
            push!(offset, (jj-1)*ones(length(xp)))
        end
        dp = DataFrame(:offset=>vcat(offset...), pr=>vcat(xx...))
        push!(ypx, [dp, vcat(yh...)])
    end

    ifig = plot_fitted(gt, ypx, trx, pr, true, ifig)

    return ifig
end

function plot_fitted(gt, ypx, trx, pr, loc, ifig)

    # Plot fitted logistic regression
    for dfp in [1, 2, 3]
        PyPlot.clf()
        PyPlot.figure(figsize=(7, 5))
        PyPlot.axes([0.16, 0.13, 0.72, 0.8])
        PyPlot.grid(true)
        for i in eachindex(ypx)
            dp, yp = ypx[i]
            if ismissing(dp)
                continue
            end
            c = tab10((i-1)/length(ypx))

            i0 = findall(dp[:, :offset] .== 0)
            i1 = findall(dp[:, :offset] .== 1)

            if dfp == 1
                yd = yp[i0] - yp[i1]
                PyPlot.plot(dp[i0, pr], yd, "-", color=c, label=@sprintf("%.0f", trx[i]))
            elseif dfp == 2
                yr = yp[i0] ./ yp[i1]
                ii = findall(yr .< 5)
                PyPlot.plot(dp[i0[ii], pr], yr[ii], "-", color=c, label=@sprintf("%.0f", trx[i]))
            elseif dfp == 3
                PyPlot.plot(dp[i0, pr], yp[i0], "-", color=c, label=@sprintf("%.0f", trx[i]))
                PyPlot.plot(dp[i1, pr], yp[i1], "--", color=c)
            else
                error("")
            end
        end
        ha, lb = PyPlot.gca().get_legend_handles_labels()
        leg = PyPlot.figlegend(ha, lb, "center right")
        leg.draw_frame(false)
        leg.set_title("T")
        if loc
            PyPlot.title("Local logistic regression")
        else
            PyPlot.title("Logistic regression")
        end

        if gt == "FGGS"
            PyPlot.xlim(0, 0.3)
        else
            PyPlot.xlim(0, 0.1)
        end

        if dfp == 1
            PyPlot.ylim(-0.3, 0.8)
            PyPlot.ylabel("Difference in probabilities of T or more\n$(qlabel[Symbol(gt)]) glomeruli in biopsy", size=14)
        elseif dfp == 2
            PyPlot.ylim(0, 5)
            PyPlot.ylabel("Ratio of probabilities of T or more\n$(qlabel[Symbol(gt)]) glomeruli in biopsy", size=14)
        elseif dfp == 3
            PyPlot.ylim(0, 1)
            PyPlot.ylabel("Probability of T or more\n$(qlabel[Symbol(gt)]) glomeruli in biopsy", size=14)
        else
            error("")
        end
        PyPlot.xlabel("Proportion of $(qlabel[Symbol(gt)]) glomeruli in nephrectomy", size=14)
        PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
        ifig += 1
    end

    return ifig
end

function fitmodels_logistic(gt, trx, cnt, fpx, ifig)

    fp1, fp2, fp3 = fpx
    nc, pr, pr1, pr2, pr3 = fitmodels_getnames(gt)

    funcs = [pr1, pr2, pr3]

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
        f = term(tnc) ~ sum(term.(funcs))
        m0 = gee(f, cnt, cnt[:, :IDx], LogitLink(), BinomialVar(), IndependenceCor())
        write(out, @sprintf(">= %.0f %s glomeruli\n\n", tr, gt))
        write(out, @sprintf(">= %.0f %s glomeruli\n\n", tr, gt))
        write(out, @sprintf("%d nephrectomies\n", length(unique(cnt[:, :IDx]))))
        write(out, @sprintf("%d observations\n", size(cnt[:, :IDx])[1]))
        write(out, string(m0))
        write(out, "\n\n")

        # Model with offset
        ff = Any[term(:offset)]
        for p in funcs
            push!(ff, term(p))
            push!(ff, term(p) & term(:offset))
        end
        f = term(tnc) ~ sum(ff)
        m1x = gee(f, cnt, cnt[:, :IDx], LogitLink(), BinomialVar(), IndependenceCor(); dofit=false)
        m1 = gee(f, cnt, cnt[:, :IDx], LogitLink(), BinomialVar(), IndependenceCor())
        write(out, @sprintf(">= %.0f %s glomeruli\n\n", tr, gt))
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

    ifig = plot_fitted(gt, ypx, trx, pr, false, ifig)

    return ifig
end

function fitmodels_helper(gt, cnt, ifig)
    cnt, fp1, fp2, fp3 = fitmodels_setup(gt, cnt)
    ifig = fitmodels_logistic(gt, [1., 2, 3, 4, 5], cnt, [fp1, fp2, fp3], ifig)
    ifig = fitmodels_local(gt, [1., 2, 3, 4, 5], cnt, ifig)
    return ifig
end

function fitmodels(cnt, ifig)
    for gt in ["FGGS", "FSGS", "Imploding", "Ischemic", "Normal"]
        ifig = fitmodels_helper(gt, cnt, ifig)
    end
    return ifig
end

function specimen_marginals(cnt, ifig)

    v = [:FGGS, :FSGS, :Imploding, :Ischemic, :Normal]
    vm = [Symbol("$(x)_total") for x in v]
    cc = combine(groupby(cnt, :ID), first)[:, vm]
    for x in names(cc)
        if occursin("total", string(x))
            cc = rename(cc, x=>replace(string(x), "_total"=>""))
        end
    end
    cc[:, :Atypical] = cc[:, :FGGS] + cc[:, :FSGS] + cc[:, :Imploding] + cc[:, :Ischemic]
    la0 = [:FGGS, :FSGS, :Imploding, :Ischemic, :Atypical, :Normal]
    cc = cc[:, la0]
    la = replace(la0, :FGGS=>:GSG, :FSGS=>:SSG)

    for pct in [false, true]
        if pct
            cc[:, :total] = cc[:, :Atypical] + cc[:, :Normal]
            for x in [:FGGS, :FSGS, :Ischemic, :Imploding, :Atypical, :Normal]
                cc[!, x] = 100 * cc[:, x] ./ cc[:, :total]
            end
        end
        for ff in [false, true]
            PyPlot.clf()
            PyPlot.figure(figsize=(6, 4))
            PyPlot.axes([0.12, 0.1, 0.8, 0.8])
            if ff
                cc1 = cc[:, la0]
                PyPlot.boxplot([x for x in eachcol(cc1)], labels=la)
            else
                # Don't plot normal and atypical
                cc1 = cc[:, la0[1:end-2]]
                PyPlot.boxplot([x for x in eachcol(cc1)], labels=la[1:end-2])
            end
            if pct
                PyPlot.ylabel("Percent of glomeruli per nephrectomy")
            else
                PyPlot.ylabel("Number of glomeruli per nephrectomy")
            end
            PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
            ifig += 1
        end
    end

    return ifig
end


function biopsy_marginals(cnt, ifig)

    cnt1 = filter(r->r.offset == 1, cnt)
    cnt0 = filter(r->r.offset == 0, cnt)

    for stat in [median, mean, minimum, maximum]
        v = [:FGGS, :FSGS, :Imploding, :Ischemic, :Normal]
        cc = []
        for offset in [0, 1]
            cnt = offset == 1 ? cnt1 : cnt0
            pp = [Symbol("$(x)_captured")=>stat for x in v]
            c = combine(groupby(cnt, :ID), pp...)
            vm = [Symbol("$(x)_captured_$(stat)") for x in v]
            push!(cc, c[:, vm])
        end

        # Interleave
        cx, la = [], []
        for x in v
            u = Symbol("$(x)_captured_$(stat)")
            push!(cx, cc[1][:, u])
            push!(la, @sprintf("%s", string(x)))
            push!(cx, cc[2][:, u])
            push!(la, @sprintf("%s\ndeep", string(x)))
        end
        for i in eachindex(la)
            la[i] = replace(la[i], "FGGS"=>"GSG")
            la[i] = replace(la[i], "FSGS"=>"SSG")
            la[i] = replace(la[i], "Imploding"=>"Impl.")
            la[i] = replace(la[i], "Ischemic"=>"Isch.")
        end

        PyPlot.clf()
        PyPlot.figure(figsize=(7, 5))
        PyPlot.axes([0.12, 0.15, 0.8, 0.7])
        PyPlot.boxplot(cx, labels=la)
        PyPlot.ylabel(@sprintf("%s number of captured\nglomeruli per biopsy", titlecase(string(stat))), size=12)

        #for x in PyPlot.gca().get_xticklabels()
        #	x.set_rotation(-90)
        #end

        PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
        ifig += 1
    end

    return ifig
end

function biopsy_scatterplots(cnt, ifig)

    cnt1 = filter(r->r.offset == 1, cnt)
    cnt0 = filter(r->r.offset == 0, cnt)

    for stat in [median, mean, minimum, maximum]
        v = [:FGGS, :FSGS, :Imploding, :Ischemic, :Normal]
        cc = []
        for offset in [0, 1]
            cnt = offset == 1 ? cnt1 : cnt0
            pp = [Symbol("$(x)_captured")=>stat for x in v]
            c = combine(groupby(cnt, :ID), pp...)
            vm = [Symbol("$(x)_captured_$(stat)") for x in v]
            push!(cc, c[:, vm])
        end

        for x in v
            PyPlot.clf()
            PyPlot.figure(figsize=(6, 4))
            PyPlot.axes([0.16, 0.18, 0.72, 0.7])
            PyPlot.grid(true)
            s = Symbol("$(x)_captured_$(stat)")
            PyPlot.plot(cc[1][:, s], cc[2][:, s], "o", mfc="none")
            PyPlot.axline((0, 0), slope=1, color="black")
            PyPlot.xlabel(@sprintf("%s number of captured %s\nglomeruli per biopsy",
                          titlecase(string(stat)), qlabel[x]), size=12)
            PyPlot.ylabel(@sprintf("%s number of captured %s\nglomeruli per deep biopsy",
                          titlecase(string(stat)), qlabel[x]), size=12)
            PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
            ifig += 1
        end
    end

    return ifig
end

cnt = load_counts()

ifig = 0
ifig = specimen_marginals(cnt, ifig)
ifig = biopsy_marginals(cnt, ifig)
ifig = biopsy_scatterplots(cnt, ifig)
ifig = fitmodels(cnt, ifig)

f = [@sprintf("plots/%03d.pdf", j) for j = 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dAutoRotatePages=/None -sOutputFile=biopsy_model.pdf $f`
run(c)

close(out)
