using GZip, CSV, DataFrames, IterTools, LinearAlgebra, UnicodePlots, Printf, PyPlot, Distributions

rm("plots", force=true, recursive=true)
mkdir("plots")

include("defs.jl")
include("pair_corr_utils.jl")
include("clinical_utils.jl")

# Pairwise correlation quantiles
idpcq, pcq = get_normalized_paircorr()

function analyze(vname, ifig, out)

    y, x = get_response(vname, idpcq, pcq)

	for j in 1:size(x,2)
	    x[:, j] = x[:, j] .- mean(x[:, j])
	end
	y = (y .- mean(y)) ./ std(y)

	u, s, v = svd(x)
	sd = std(u[:, 1])

	p = size(x, 2)
	F = zeros(p-2, p)
	for i in 1:p-2
	    F[i, i:i+2] = [1, -2, 1]
	end
	G = F' * F

	println(vname)
	bl = []
	for la in [0.1, 1.0, 10., 100.]
        pm = x'*x + la*G
	    b = pm \ (x' * y)
	    push!(bl, b)
	    yh = x * b
	    sig2 = mean((yh - y).^2)
	    cm = sig2 * (pm \ (x' * x) / pm)
	    cs = b' * (cm \ b)
	    pv = 1 - cdf(Chisq(20), cs)
    	r = cor(yh, y)
    	write(out, @sprintf("%s,%.2f,%.3f,%.3f\n", vname, la, pv, r))
	end

    # Plot the PC loading vector
	PyPlot.clf()
	PyPlot.axes([0.13, 0.12, 0.75, 0.8])
	PyPlot.grid(true)
	pr = collect(range(0, 1, length=size(v, 1)))
	pr = pr ./ maximum(pr)
	cm = PyPlot.cm.get_cmap("jet")
	for (j,b) in enumerate(bl)
	    col = cm(j/(length(bl)+1))
    	PyPlot.plot(pr, b, "-", color=col)
    end
    PyPlot.title(vname)
	PyPlot.xlabel("Probability point", size=15)
	PyPlot.ylabel("Coefficient", size=15)
	PyPlot.savefig(@sprintf("plots/%03d.pdf", ifig))
	ifig += 1

	return ifig

end

avn = [:Female, :NonWhite, :age, :htn_code_pre, :dm_code_pre, :bmi, :htn_sbp, :htn,
       :dm_lab_pre, :dm, :bl_cr_prenx, :bl_egfr_prenx, :bl_ckd_prenx, :kdigo,
       :"Hyper-Bin", :"T2D-Bin", :"eGFR Min4 Mean", :"I/L Ratio", :">100%", :">50%",
       :NormalGlomPercent, :AbnormalPercent, :"AKI-Bin", :GGS, :Imploding, :GlomVol,
       :"% CA Glomerular", :"Nephron Density (All)", :"Nephron Density (Normal)",
       :"GV/Podocyte", :PodoDensity, :MesIndex, :PodoPerGlom, :MesVol, :PodoVol,
       :FIA, :SGS, :Ischemic, :vGlepp]

function main()
    ifig = 0
    out = open("clinical_reg_results.csv", "w")
    write(out, "Variable,Lambda,P-value,Correlation\n")
    for av in avn
	    ifig = analyze(av, ifig, out)
    end
    close(out)
    return ifig
end

ifig = main()

f = [@sprintf("plots/%03d.pdf", j) for j in 0:ifig-1]
c = `gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=clinical_reg.pdf $f`
run(c)
