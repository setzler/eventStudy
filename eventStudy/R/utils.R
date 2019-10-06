#' @export
# modified max() function that returns 0 instead of -Inf in the case where all entries are NA and na.rm = TRUE
# otherwise, behaves just like max()
# used in ES_clean_data() when subsetting to valid treated / control units using, e.g., treated_subset_var and control_subset_var
max2 <- function(x, na.rm){
  y <- suppressWarnings(max(x, na.rm = na.rm))
  if(y==-Inf){
    y <- 0
  }
  return(y)
}

#' @export
# modified sum() function that returns NA instead of 0 in the case where all entries are NA and na.rm = TRUE
# otherwise, behaves just like sum()
# used in ES() at a variety of places to correctly assigned ref_event_time-specific covariates if they are present at the desired ref_event_time
sum2 <- function(x, na.rm){
  if (all(is.na(x))){
    y <- NA
  } else{
    y <- sum(x, na.rm = na.rm)
  }
  return(y)
}

#' @export
# function to draw a block bootstrap sample along provided 'unit_var'
# Thanks (without implicating) to Michael Graber for sharing this function
block_sample <- function(long_data, unit_var, cal_time_var){
  # vector of unique cluster IDs and total number of clusters
  setkeyv(long_data, c(unit_var, cal_time_var))
  ids <- unique(long_data[, get(unit_var)])
  N   <- length(ids)
  # sample ID's with replacement:
  # we sample ID's with replacement, then collapse the resulting draws to
  # a vector of unique IDs with a column representing the number of times the cluster
  # has been sampled. e.g. if unit_var = c("A","B") and count = c(1,3) then
  # cluster "A" has been sampled once, while cluster "B" has been sampled three times.
  clusters <- data.table(sample(ids, size = N, replace = TRUE))
  names(clusters) <- unit_var
  setkeyv(clusters, unit_var)
  clusters[,count := .N, by = get(unit_var)]
  clusters <- unique(clusters)
  # repeated merges to create block bootstrap sample (bs)
  max_iter <- max(clusters$count)
  bs       <- data.table()
  for (i in 1 : max_iter) {
    tmp <- merge(clusters[count >= i,], long_data, by = unit_var, all.x = TRUE,
                 all.y = FALSE, sort = FALSE)
    tmp[, eval(unit_var) := paste0(get(unit_var), "-", as.character(i))]  # create new unique ID: original ID "-" index of draw
    bs <- rbind(bs, tmp)
  }
  bs[,count := NULL]
  setkeyv(bs, c(unit_var, cal_time_var))
  gc()

  return(bs)
}

#' @export
# Adapted near-verbatim from 'msm' package
# g         # a formula or list of formulae (functions) giving the transformation g(x) in terms of x1, x2,  etc
# mean      # mean, MLE, or other consistent plug-in estimate of x
# cov       # covariance matrix of x
# ses=TRUE  # return standard errors, else return covariance matrix
delta_method <- function (g, mean, cov, ses = TRUE) {
  ## Var (G(x))  =  D_g * Var(X) * t(D_g)
  cov <- as.matrix(cov)
  n <- length(mean)
  if (!is.list(g)){
    g <- list(g)
  }
  if ((dim(cov)[1] != n) || (dim(cov)[2] != n)){
    stop(print(sprintf("'cov' is a %s X %s matrix, but should be a square %s X %s matrix.", dim(cov)[1], dim(cov)[2], n, n)))
  }
  syms <- paste("x", 1:n, sep = "")
  for (i in 1:n) assign(syms[i], mean[i])
  D_g <- t(sapply(g, function(form) {
    ## 1) differentiate each formula in the list
    ## 2) evaluate at the supplied estimated / plug-in value
    ## 3) take these elements and make Jacobian row-by-row
    as.numeric(attr(eval(deriv(form, syms)), "gradient"))
  }))
  new_cov <- D_g %*% cov %*% t(D_g)
  if (ses == TRUE) {
    result <- sqrt(diag(new_cov))
  } else{
    result <- new_cov
  }

  return(result)
}
