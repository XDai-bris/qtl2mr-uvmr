# Session state for reusing exact intermediate results. The state lives in the
# package namespace, never in resources or results, so saveRDS() of user
# objects is unchanged and nothing here outlives the session.
qtl2mr_state <- new.env(parent = emptyenv())

qtl2mr_memo_list <- function(name) qtl2mr_state[[name]] %||% list()

qtl2mr_memo_push <- function(name, entry, max_entries) {
  memo <- c(list(entry), qtl2mr_memo_list(name))
  if (length(memo) > max_entries) memo <- memo[seq_len(max_entries)]
  assign(name, memo, envir = qtl2mr_state)
  invisible(entry)
}

# Kinship matrices that passed the matrix-only checks of validation. identical()
# returns at once for the same object, so a resources kinship is checked once.
qtl2mr_kinship_checked <- function(K) {
  for (k in qtl2mr_memo_list("kinship_checked")) if (identical(k, K)) return(TRUE)
  FALSE
}

qtl2mr_kinship_checked_add <- function(K) {
  qtl2mr_memo_push("kinship_checked", K, 64L)
}

