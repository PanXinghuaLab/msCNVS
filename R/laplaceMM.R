#' Fixed Mean Mixture Density Estimation with Background Absorption
#'
#' Fits a mixture model with fixed means to data, including a background component.
#'
#' @param x Numeric vector of data points.
#' @param means Numeric vector of fixed means for the peaks.
#' @param offmeans Optional numeric vector of additional (off-peak) means.
#' @param weights Optional numeric vector of weights for each data point.
#' @param dist_type Distribution type for peaks ("gaussian", "laplace").
#' @param background_type Background distribution type ("uniform", "cauchy").
#' @param max_iter Maximum number of EM iterations (default: 1000).
#' @param tol Convergence tolerance (default: 1e-6).
#' @param core_prob Central probability used to define core intervals (default: 0.95).
#' @return A list containing model parameters and results, including classification.
#'
#' @export
fixed_mean_mixture <- function(x, means, offmeans , weights = NULL,
                               dist_type = "gaussian",
                               background_type = "uniform",
                               max_iter = 1000, tol = 1e-6,
                               core_prob = 0.95) {
  # Input validation and initialization
  if (is.null(weights)) weights <- rep(1, length(x))
  if (length(weights) != length(x)) stop("weights length must match x length")
  if (!dist_type %in% c("gaussian", "laplace")) stop("Invalid dist_type")
  if (!background_type %in% c("uniform", "cauchy")) stop("Invalid background_type")
  n <- length(x)
  k_peaks     <- length(means)
  k_offpeaks  <- length(offmeans)
  k_total     <- k_peaks + 1 + ifelse(k_offpeaks > 0, k_offpeaks, 0)

  # Initialize parameters
  params <- list(
    weights_comp = rep(1/k_total, k_total),
    scale_params = rep(0.1, k_peaks),
    bg_params = if (background_type == "uniform") {
      list(min = min(x), max = max(x))
    } else {
      list(location = stats::median(x), scale = stats::IQR(x)/2)
    }
  )
  if(k_offpeaks){
    params$mix_params <-  rep(0.05,k_offpeaks)
  }

  # EM loop
  em_result <- run_em_algorithm(x, means, offmeans, weights, params,
                                dist_type, background_type,
                                max_iter, tol)


  processed_result <-process_em_results(x, means, offmeans, em_result,
                                        core_prob = core_prob)

  # Build final result
  result <- list(
    peak_params = data.frame(
      mean = means,
      scale = em_result$scale_params,
      weight = em_result$weights_comp[1:k_peaks]
    ),
    background_params = c(em_result$bg_params, weight = em_result$weights_comp[k_total]),
    gamma = em_result$gamma,
    classification = processed_result$classification,
    core_intervals = processed_result$core_intervals,
    log_lik = em_result$log_lik,
    converged = em_result$converged,
    iterations = em_result$iter,
    data = x,
    weights = weights,
    dist_type = dist_type,
    background_type = background_type
  )
  return(result)
}

# EM algorithm
run_em_algorithm <- function(x, means, offmeans, weights, params,
                             dist_type, background_type,
                             max_iter, tol) {
  n <- length(x)
  k_peaks     <- length(means)
  k_offpeaks  <- length(offmeans)
  k_total     <- k_peaks + 1 + ifelse(k_offpeaks > 0, k_offpeaks, 0)

  log_lik <- -Inf
  converged <- FALSE
  gamma <- matrix(1/k_total, nrow = n, ncol = k_total)
  for (iter in 1:max_iter) {
    # E-step: compute responsibilities
    log_dens <- calculate_log_densities(x, means, offmeans, params, dist_type, background_type)
    e_step <- perform_e_step(log_dens, params$weights_comp, weights)
    gamma <- e_step$gamma

    # Check convergence
    new_log_lik <- e_step$log_lik
    if (abs(new_log_lik - log_lik) < tol) {
      converged <- TRUE
      break
    }
    log_lik <- new_log_lik

    # M-step: update parameters
    params <- perform_m_step(x, means, offmeans, weights, gamma, params,
                             dist_type, background_type)
  }

  return(list(
    weights_comp = params$weights_comp,
    scale_params = params$scale_params,
    bg_params = params$bg_params,
    mix_params = params$mix_params,
    gamma = gamma,
    log_lik = log_lik,
    converged = converged,
    iter = iter
  ))
}

# Compute log densities for components
calculate_log_densities <- function(x, means,offmeans, params, dist_type, background_type) {

  n <- length(x)
  k_peaks     <- length(means)
  k_offpeaks  <- length(offmeans)
  k_total     <- k_peaks + 1 + ifelse(k_offpeaks > 0, k_offpeaks, 0)

  log_dens <- matrix(0, nrow = length(x), ncol = k_total)

  # Peak components
  for (j in 1:k_peaks) {
    if (dist_type == "gaussian") {
      log_dens[, j] <- stats::dnorm(x, means[j], params$scale_params[j], log = TRUE)
    } else {
      log_dens[, j] <- log(0.5/params$scale_params[j]) - abs(x - means[j]) / params$scale_params[j]
    }
  }

  # Background component
  if (background_type == "uniform") {
    log_dens[, k_total] <- stats::dunif(x, params$bg_params$min, params$bg_params$max, log = TRUE)
  } else {
    log_dens[, k_total ] <- stats::dcauchy(x, params$bg_params$location,
                                    params$bg_params$scale, log = TRUE)
  }
  if(k_offpeaks) {
    for(j in 1:k_offpeaks){
      log_dens[, k_peaks+j] <- log(0.5/params$mix_params[j]) - abs(x - offmeans[j]) / params$mix_params[j]
      #log_dens[, k_peaks+j] <- log(1 / (2 * params$mix_params$scale[j])) - abs(x - params$mix_params$mu[j]) / params$mix_params$scale[j]
    }
  }

  return(log_dens)
}

# E-step implementation
perform_e_step <- function(log_dens, weights_comp, weights) {
  if (length(weights_comp) != ncol(log_dens))
    stop("weights_comp length must equal number of columns in log_dens")

  weighted_log_dens <- sweep(log_dens, 2, log(weights_comp), "+")
  max_log_dens <- apply(weighted_log_dens, 1, max)
  centered <- sweep(weighted_log_dens, 1, max_log_dens, "-")
  gamma <- exp(centered)
  rs <- rowSums(gamma)
  gamma <- gamma / rs

  log_lik <- sum(weights * (max_log_dens + log(rs)))
  return(list(gamma = gamma, log_lik = log_lik))
}

# M-step implementation
perform_m_step <- function(x, means,offmeans, weights, gamma, params,
                           dist_type, background_type) {
  k_peaks     <- length(means)
  k_offpeaks  <- length(offmeans)
  k_total     <- k_peaks + 1 + ifelse(k_offpeaks > 0, k_offpeaks, 0)

  Nk <- colSums(weights * gamma)

  # Update mixture weights
  params$weights_comp <- Nk / sum(Nk)

  # Update peak scales
  for (j in 1:k_peaks) {
    if (dist_type == "gaussian") {
      params$scale_params[j] <- sqrt(sum(weights * gamma[, j] * (x - means[j])^2) / Nk[j])
    } else {
      params$scale_params[j] <- sum(weights * gamma[, j] * abs(x - means[j])) / Nk[j]
    }
  }

  # Update background params
  if(k_offpeaks>0){
    for(j in 1:k_offpeaks){
      #params$mix_params$mu[j] <- matrixStats::weightedMedian(x, w = weights * gamma[, j])
      params$mix_params[j] <- sum(weights * gamma[,k_peaks + j] * abs(x - offmeans[j])) / Nk[k_peaks+j]
    }
  }
  if (background_type == "cauchy") {
    bg_opt <- stats::optim(
      par = c(params$bg_params$location, log(params$bg_params$scale)),
      fn = function(par) {
        loc <- par[1]
        sc <- exp(par[2])
        -sum(weights * gamma[, k_total] * stats::dcauchy(x, loc, sc, log = TRUE))
      }
    )
    params$bg_params$location <- bg_opt$par[1]
    params$bg_params$scale <- exp(bg_opt$par[2])
  }

  return(params)
}

process_em_results <- function(x, means,offmeans, em_result,
                               core_prob ) {
  # Input checks
  if (length(x) != nrow(em_result$gamma)) {
    stop("Length of data does not match gamma matrix rows")
  }
  k_peaks     <- length(means)
  k_offpeaks  <- length(offmeans)
  k_total     <- k_peaks + 1 + ifelse(k_offpeaks > 0, k_offpeaks, 0)

  # 1. Classification ------------------------------------------------------
  gamma_peaks <- em_result$gamma
  classification <- apply(gamma_peaks, 1, which.max)

  # 2. Core intervals ------------------------------------------------------
  calc_core_interval <- function(points, prob) {
    if (length(points) > 0){
      stats::quantile(points, c((1-prob)/2, 1-(1-prob)/2))
    } else c(NA, NA)
  }
  id <-  unique(classification)
  id <- id[ id %in% 1:k_peaks]
  core_intervals <- t(sapply(id, function(j) {
    calc_core_interval(x[classification == j], core_prob)
  }))
  core_intervals <- as.data.frame(core_intervals)
  colnames(core_intervals) <- c("lower", "upper")
  row.names(core_intervals) <- means[id]

  # 3. Return --------------------------------------------------------------
  list(
    classification = classification,  # point class (1:k_peaks are peaks; k_peaks+1 is background)
    core_intervals = core_intervals  # core interval for each peak

  )
}

#' @rdname laplaceMM
#' @export
setMethod("laplaceMM", "CNVisionObj", function(object,core_prob = 0.95,dist_type = "gaussian",background_type = "uniform") {
  #df <- do.call(rbind,object@result$ScaledCNV)
  #df <- df[df$seg.mean>= 0.5 & df$seg.mean <= 6.5,]
  df <- object@result$df
  model <- fixed_mean_mixture(x =  df$seg.mean ,
                              weights = df$num.mark ,
                              means = object@result$kde_peaks$peak_x ,
                              offmeans = object@result$offpeaks$peak_x ,
                              dist_type = dist_type ,
                              background_type = background_type,
                              core_prob = core_prob)
  set_values_matrix <- function(x, intervals, out_val = 0) {
    res <-data.frame(seg.mean = x,lab = rep(out_val, length(x)))
    for (i in seq_len(nrow(intervals))) {
      res[x >= intervals[i, 1] & x <= intervals[i, 2],] <- round(as.numeric(row.names(intervals)))[i]

    }
    return(res)
  }
  AbsoluteCN <- lapply(object@result$ScaledCNV,function(df){

    df2 <- set_values_matrix(x = df$seg.mean,intervals =  model$core_intervals,out_val = 0)
    #df2$num.mark <- df$num.mark
    df <- df %>%
      dplyr::select(-dplyr::any_of(names(df2))) %>%
      dplyr::bind_cols(df2)
    return(df)
  })

  stats <- lapply(object@result$ScaledCNV,function(df){
    weighted_mean <- sum(df$seg.mean * df$num.mark) / sum(df$num.mark)
    weighted_var <- sum(df$num.mark * (df$seg.mean - weighted_mean)^2) / sum(df$num.mark)
    weighted_sd <- sqrt(weighted_var)
    weighted_cv <- weighted_sd / weighted_mean

    return(weighted_cv)
  })
  CV <- apply(object@cnvData$gc.ratio, 2, function(x) {
    stats::sd(x)/mean(x)
  })
  names(CV) <- object@config$cells

  Nrcd <- lapply(object@config$cells,function(cell){
    gc.ratio <- scale_transform(object@cnvData$gc.ratio,object, cell)
    x <- gc.ratio- rep(AbsoluteCN[[cell]]$seg.mean,AbsoluteCN[[cell]]$num.mark)
    x <- x[gc.ratio<=6.5]
    x <- sum(x^2)/object@config$bin_number/object@result$ploidy[[cell]]
    return(x)
  })
  names(Nrcd) <- object@config$cells

  object@result$AbsoluteCN <-AbsoluteCN
  object@result$model <- model
  object@result$stats <- stats
  object@result$CV <- CV
  object@result$Nrcd <- Nrcd
  return(object)
})
