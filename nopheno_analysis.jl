using Printf, Statistics, IterTools, LinearAlgebra

include("defs.jl")
include("annot_utils.jl")
include("pair_corr_utils.jl")

annotsx = glom_centroids(annots)
annotsx = Dict(k => condense(v) for (k, v) in annotsx)

# Only consider samples with a minimum number of total glomeruli
annotsx = Dict(k => v for (k, v) in annotsx if length(v["All_glomeruli"]) > 100)

# Pairwise correlation quantiles
# pcq are the log ratios of atypical/atypical distances versus typical/typical distances
# pcqn are the atypical/atypical distances
# pcqd are the typical/typical distances
scid, pcq, pcqn, pcqd = get_normalized_paircorr(annotsx)

n = size(pcq, 1)
z = sqrt(n) * mean(pcq, dims = 1) ./ std(pcq, dims = 1)
