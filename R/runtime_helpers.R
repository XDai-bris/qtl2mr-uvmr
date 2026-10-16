qtl2mr_default_output_root <- function() {
  opt <- getOption("qtl2mr.output_root")
  if (!is.null(opt) && nzchar(opt)) return(normalizePath(opt, mustWork = FALSE))
  env <- Sys.getenv("QTL2MR_OUTPUT_ROOT", unset = "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  tempdir()
}

qtl2mr_get_tophits <- function(lod, lod_thresh = 4) {
  lod_thresh <- qtl2mr_validate_lod_thresh(lod_thresh)
  lod <- as.numeric(lod)
  lod[!is.finite(lod) | lod < lod_thresh] <- 0
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
