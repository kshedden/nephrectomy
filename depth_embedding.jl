include("depth_utils.jl")

# Set up the glomeruli
annotsx = glom_centroids(annots)
annotsx = major_components(annotsx)
annotsc = Dict(k => condense(v) for (k, v) in annotsx)

# Calculate depth quantiles for these categories of glomeruli,
# always relative to normal references.
gt = ["Atypical", "Normal"]

# Calculate depth quantiles at a sequence of probability points
pp = range(0.1, 0.9, length=9)
dd = build_depths(pp, l2_depth, annotsc, gt)
dd = dd[completecases(dd), :]

# Depths of atypical glomeruli
na = [@sprintf("Atypical_%.2f", p) for p in pp]
da = dd[:, na]

# Depths of normal glomeruli
na = [@sprintf("Normal_%.2f", p) for p in pp]
dn = dd[:, na]

# Log ratio of atypical depths to normal depths
anr = log.(Matrix(da) ./ Matrix(dn))

id = dd[:, :Scanner_ID]
