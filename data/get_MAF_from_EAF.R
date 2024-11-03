eaf2maf <- function(eaf = NULL) {
  if (any(is.infinite(eaf))) {
    warning("The 'eaf' vector contains infinite values. These will be turned into NAs.")
    is.na(eaf) <- is.infinite(eaf)
  }
  if (is.null(eaf)) {
    stop("No 'eaf' vector provided.")
  }
  if (is.character(eaf)) {
    stop("'eaf' must be a numeric vector.")
  }
  if (all(is.na(eaf))) {
    stop("All values in 'eaf' are NA.")
  }
  if (is.numeric(eaf)) {
    maf <- eaf
    ind <- which(eaf > 0.5)
    maf[ind] <- 1 - eaf[ind]
    return(maf)
  }
}