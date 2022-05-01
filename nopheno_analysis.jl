using Printf, Base.Iterators

include("defs.jl")
include("annot_utils.jl")
include("pair_corr_utils.jl")

# Pairwise correlation quantiles
# pcq are the log ratios of atypical/atypical distances versus typical/typical distances
# pcqn are the atypical/atypical distances
# pcqd are the typical/typical distances
idpcq, pcq, pcqn, pcqd = get_normalized_paircorr(annots, nothing)
