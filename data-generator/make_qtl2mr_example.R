# Deterministic synthetic teaching fixture for qtl2mr.
#
# This script creates data/qtl2mr_example.rda. It contains no animal,
# application, or biological data. The fixed seed and generating mechanism are
# deliberately recorded so the fixture can be inspected or regenerated.

set.seed(301)
n <- 96L
founders <- LETTERS[1:8]
markers <- paste0("m", 1:4)
ids <- paste0("id", seq_len(n))

allele_array <- array(
  NA_real_,
  dim = c(n, length(founders), length(markers)),
  dimnames = list(ids, founders, markers)
)
for (j in seq_along(markers)) {
  raw_probabilities <- matrix(stats::rexp(n * length(founders)), nrow = n)
  allele_array[, , j] <- raw_probabilities / rowSums(raw_probabilities)
}
allele_probs <- list("1" = allele_array)
class(allele_probs) <- c("calc_genoprob", "list")

helmert <- stats::contr.helmert(length(founders))
helmert <- apply(helmert, 2L, function(x) x / sqrt(sum(x^2)))
H <- cbind(mean = rep(1, length(founders)) / sqrt(length(founders)), helmert)
rownames(H) <- founders
colnames(H) <- c("mean", paste0("FC", 1:7))
z <- lapply(seq_along(markers), function(j) allele_array[, , j] %*% H)

true_theta <- 0.75
x <- as.numeric(
  1.4 * z[[1L]][, "FC1"] + 1.1 * z[[2L]][, "FC2"] +
    0.9 * z[[3L]][, "FC3"] + stats::rnorm(n, sd = 0.25)
)
y <- as.numeric(true_theta * x + 0.2 * z[[4L]][, "FC4"] + stats::rnorm(n, sd = 0.25))

qtl2mr_example <- list(
  allele_probs = allele_probs,
  pheno = data.frame(x = x, y = y, row.names = ids),
  covar = data.frame(
    sex = rep(0:1, length.out = n),
    age = as.numeric(scale(seq_len(n))),
    row.names = ids
  ),
  kinship = NULL,
  marker_map = data.frame(
    marker = markers,
    chr = "1",
    pos = seq_along(markers),
    exposure_lod = c(8.0, 7.2, 4.2, 3.5),
    stringsAsFactors = FALSE
  ),
  selected_loci = c("m1", "m2"),
  exposure = "x",
  outcome = "y",
  true_theta = true_theta,
  metadata = list(
    purpose = "synthetic teaching and integration fixture only",
    seed = 301L,
    founders = founders,
    n_samples = n,
    kinship_note = "NULL only for the small self-contained teaching fixture; use harmonized LOCO kinship for real DO analyses.",
    exposure_lod_note = "Illustrative marker labels only; these values are not qtl2-derived calibration."
  )
)

save(qtl2mr_example, file = "data/qtl2mr_example.rda", compress = "xz")
