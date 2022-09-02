using DataFrames, LinearAlgebra, Printf, StaticArrays, Statistics

include("defs.jl")
include("annot_utils.jl")
include("clinical_utils.jl")
include("depth_utils.jl")

function dfclean(df)
    x = split(string(df), "\n")
    x = vcat(x[2], x[4:end])
    x = join(x, "\n")
    return x
end

function main(annots)
    out = open("depth.txt", "w")

    pp = 0.5
    gt = [@sprintf("%s_%.2f", g, pp) for g in glom_types]

    rkeys = ["Normal"]#, "Normal", "Capsule", "CMJ"]

    dn = ["L2 depth"]#, "Tukey halfspace depth"]#, "Distance to capsule", "Distance to CMJ"]
    dpf = [l2_depth]#, tukey_depth] #, boundary_depth, boundary_depth]

    for (jd, depthfun) in enumerate(dpf)

        depths = build_depths([0.5], depthfun, annots, glom_types)

        clin[!, :Race] = [ismissing(x) ? missing : Int(x) for x in clin[:, :Race]]
        for x in unique(clin[:, :Race])
            clin[:, Symbol("Race$(x)")] = clin[:, :Race] .== x
        end
        depths = leftjoin(depths, clin, on = :Scanner_ID)

        drslt = depth_analysis(depths, gt)
        crslt = clinical_analysis(depths, gt)
        write(out, @sprintf("=== %s ===\n\n", dn[jd]))
        write(out, "Differences in median depth based on glomerulus type:\n")
        write(out, dfclean(drslt))
        write(out, "\n\nAssociations between depth and clinical variables:\n")
        write(out, dfclean(crslt))
        write(out, "\n\n")
    end

    close(out)
end

annotsx = glom_centroids(annots)
annotsx = major_components(annotsx)

main(annotsx)
