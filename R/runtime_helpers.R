qtl2mr_default_output_root <- function() {
  opt <- getOption("qtl2mr.output_root")
  if (!is.null(opt) && nzchar(opt)) return(normalizePath(opt, mustWork = FALSE))
  env <- Sys.getenv("QTL2MR_OUTPUT_ROOT", unset = "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  tempdir()
}

# Peaks of one chromosome's LOD profile, given in map order; returns indices into
# lod. Non-finite LOD values are skipped: they neither start a peak nor count as
# a dip. "peakdrop" follows qtl2::find_peaks() (plain, no intervals) except that
# the threshold is LOD >= lod_thresh; "lod_island" takes the maximum of each
# contiguous above-threshold run, the rule used before 0.7.0.
qtl2mr_get_tophits <- function(lod, lod_thresh = 4, rule = c("peakdrop", "lod_island"), peakdrop = 1.5) {
  lod_thresh <- qtl2mr_validate_lod_thresh(lod_thresh)
  rule <- match.arg(rule)
  lod <- as.numeric(lod)
  finite <- which(is.finite(lod))
  peaks <- if (identical(rule, "peakdrop")) {
    peakdrop <- qtl2mr_validate_peakdrop(peakdrop)
    qtl2mr_peakdrop_peaks(lod[finite], lod_thresh, peakdrop)
  } else {
    qtl2mr_lod_island_peaks(lod[finite], lod_thresh)
  }
  finite[peaks]
}

# A new peak starts only after the LOD has dropped by more than peakdrop below
# the current peak and risen again by more than peakdrop; otherwise a higher LOD
# moves the current peak. Ties keep the first marker in map order.
qtl2mr_peakdrop_peaks <- function(lod, lod_thresh, peakdrop) {
  peaks <- integer()
  last_peak <- NA_real_
  min_since <- Inf
  for (i in seq_along(lod)) {
    min_since <- min(min_since, lod[[i]])
    if (lod[[i]] < lod_thresh) next
    if (!length(peaks)) {
      peaks <- i
    } else if (last_peak > min_since + peakdrop && lod[[i]] > min_since + peakdrop) {
      peaks <- c(peaks, i)
    } else if (lod[[i]] > last_peak) {
      peaks[[length(peaks)]] <- i
    } else {
      next
    }
    last_peak <- lod[[i]]
    min_since <- lod[[i]]
  }
  peaks
}

# Peaks are returned in decreasing LOD, ties in map order.
qtl2mr_lod_island_peaks <- function(lod, lod_thresh) {
  lod[lod < lod_thresh] <- 0
  peaks <- integer()
  while (any(lod >= lod_thresh)) {
    idx <- which.max(lod)
    left_hit <- which(lod[idx:1L] < lod_thresh)[1L]
    right_hit <- which(lod[idx:length(lod)] < lod_thresh)[1L]
    left <- if (is.na(left_hit)) 1L else idx - left_hit + 2L
    right <- if (is.na(right_hit)) length(lod) else idx + right_hit - 2L
    peaks <- c(peaks, idx)
    lod[seq.int(max(left, 1L), min(right, length(lod)))] <- 0
  }
  peaks
}

# Greedy minimum spacing: visit peaks in decreasing LOD (ties in map order) and
# drop any within min_spacing of a peak already kept.
qtl2mr_space_peaks <- function(peaks, lod, pos, min_spacing) {
  kept <- integer()
  for (p in peaks[order(-lod[peaks], peaks)]) {
    if (!any(abs(pos[[p]] - pos[kept]) <= min_spacing)) kept <- c(kept, p)
  }
  kept
}

qtl2mr_draw_explicit_device <- function(path, expr,
                                        width = 1800,
                                        height = 1300,
                                        res = 180,
                                        pdf_width = 8.5,
                                        pdf_height = 6.5) {
  ext <- tolower(tools::file_ext(path))
  if (identical(ext, "png")) {
    grDevices::png(path, width = width, height = height, res = res)
  } else if (identical(ext, "pdf")) {
    grDevices::pdf(path, width = pdf_width, height = pdf_height)
  } else {
    stop("Unsupported figure extension: ", path, call. = FALSE)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}
