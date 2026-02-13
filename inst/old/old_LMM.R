#' Fixed Mean Mixture Density Estimation with Background Absorption
#'
#' Fits a mixture model with fixed means to data, including a background component.
#'
#' @param x Numeric vector of data points.
#' @param means Numeric vector of fixed means for the peaks.
#' @param weights Optional numeric vector of weights for each data point.
#' @param dist_type Distribution type for peaks ("gaussian", "laplace").
#' @param background_type Background distribution type ("uniform", "cauchy").
#' @param max_iter Maximum number of EM iterations (default: 100).
#' @param tol Convergence tolerance (default: 1e-6).
#' @param gamma_threshold Gamma threshold for point filtering (default: 0.5).
#' @return A list containing model parameters and results, including classification.
#'
#' @export
fixed_mean_mixture <- function(x, means, weights = NULL,
                               dist_type = "gaussian",
                               background_type = "uniform",
                               max_iter = 1000, tol = 1e-6,
                               core_prob = 0.95) {
  # 输入验证和初始化
  if (is.null(weights)) weights <- rep(1, length(x))
  if (length(weights) != length(x)) stop("weights length must match x length")
  if (!dist_type %in% c("gaussian", "laplace")) stop("Invalid dist_type")
  if (!background_type %in% c("uniform", "cauchy")) stop("Invalid background_type")

  n <- length(x)
  k_peaks <- length(means)
  k_total <- k_peaks + 1  # 峰 + 背景

  # 初始化参数
  params <- list(
    weights_comp = rep(1/k_total, k_total),
    scale_params = rep(0.1, k_peaks),
    bg_params = if (background_type == "uniform") {
      list(min = min(x), max = max(x))
    } else {
      list(location = median(x), scale = IQR(x)/2)
    }
  )

  # EM算法主循环
  em_result <- run_em_algorithm(x, means, weights, params,
                                dist_type, background_type,
                                max_iter, tol)


  processed_result <-process_em_results(x, means, em_result,
                                        core_prob = core_prob)

  # 构建最终结果
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

# EM算法实现
run_em_algorithm <- function(x, means, weights, params,
                             dist_type, background_type,
                             max_iter, tol) {
  n <- length(x)
  k_peaks <- length(means)
  k_total <- k_peaks + 1

  log_lik <- -Inf
  converged <- FALSE
  gamma <- matrix(1/k_total, nrow = n, ncol = k_total)

  for (iter in 1:max_iter) {
    # E-step: 计算gamma
    log_dens <- calculate_log_densities(x, means, params, dist_type, background_type)
    e_step <- perform_e_step(log_dens, params$weights_comp, weights)
    gamma <- e_step$gamma

    # 检查收敛
    new_log_lik <- e_step$log_lik
    if (abs(new_log_lik - log_lik) < tol) {
      converged <- TRUE
      break
    }
    log_lik <- new_log_lik

    # M-step: 更新参数
    params <- perform_m_step(x, means, weights, gamma, params,
                             dist_type, background_type)
  }

  return(list(
    weights_comp = params$weights_comp,
    scale_params = params$scale_params,
    bg_params = params$bg_params,
    gamma = gamma,
    log_lik = log_lik,
    converged = converged,
    iter = iter
  ))
}

# 计算各成分的对数密度
calculate_log_densities <- function(x, means, params, dist_type, background_type) {
  k_peaks <- length(means)
  k_total <- k_peaks + 1
  log_dens <- matrix(0, nrow = length(x), ncol = k_total)

  # 峰成分
  for (j in 1:k_peaks) {
    if (dist_type == "gaussian") {
      log_dens[, j] <- dnorm(x, means[j], params$scale_params[j], log = TRUE)
    } else {
      log_dens[, j] <- log(0.5/params$scale_params[j]) - abs(x - means[j]) / params$scale_params[j]
    }
  }

  # 背景成分
  if (background_type == "uniform") {
    log_dens[, k_total] <- dunif(x, params$bg_params$min, params$bg_params$max, log = TRUE)
  } else {
    log_dens[, k_total] <- dcauchy(x, params$bg_params$location,
                                   params$bg_params$scale, log = TRUE)
  }

  return(log_dens)
}

# E-step实现
perform_e_step <- function(log_dens, weights_comp, weights) {
  weighted_log_dens <- sweep(log_dens, 2, log(weights_comp), "+")
  max_log_dens <- apply(weighted_log_dens, 1, max)
  weighted_log_dens <- sweep(weighted_log_dens, 1, max_log_dens, "-")
  gamma <- exp(weighted_log_dens)
  gamma <- gamma / rowSums(gamma)

  log_lik <- sum(weights * (max_log_dens + log(rowSums(exp(weighted_log_dens)))))

  return(list(gamma = gamma, log_lik = log_lik))
}

# M-step实现
perform_m_step <- function(x, means, weights, gamma, params,
                           dist_type, background_type) {
  k_peaks <- length(means)
  k_total <- k_peaks + 1
  Nk <- colSums(weights * gamma)

  # 更新混合权重
  params$weights_comp <- Nk / sum(Nk)

  # 更新峰参数
  for (j in 1:k_peaks) {
    if (dist_type == "gaussian") {
      params$scale_params[j] <- sqrt(sum(weights * gamma[, j] * (x - means[j])^2) / Nk[j])
    } else {
      params$scale_params[j] <- sum(weights * gamma[, j] * abs(x - means[j])) / Nk[j]
    }
  }

  # 更新背景参数
  if (background_type == "cauchy") {
    bg_opt <- optim(
      par = c(params$bg_params$location, log(params$bg_params$scale)),
      fn = function(par) {
        loc <- par[1]
        sc <- exp(par[2])
        -sum(weights * gamma[, k_total] * dcauchy(x, loc, sc, log = TRUE))
      }
    )
    params$bg_params$location <- bg_opt$par[1]
    params$bg_params$scale <- exp(bg_opt$par[2])
  }

  return(params)
}

process_em_results <- function(x, means, em_result,
                               core_prob ) {
  # 输入检查
  if (length(x) != nrow(em_result$gamma)) {
    stop("数据长度与gamma矩阵行数不匹配")
  }

  # 1. 分类结果 ------------------------------------------------------------
  gamma_peaks <- em_result$gamma
  max_gamma <- apply(gamma_peaks, 1, max)

  # 分配类别（最大值对应的峰，低于阈值归为背景）
  classification <- apply(gamma_peaks, 1, which.max)

  # 2. 核心区间 ------------------------------------------------------------
  calc_core_interval <- function(points, prob) {
    if (length(points)){
      quantile(points, c((1-prob)/2, 1-(1-prob)/2))
    } else {
      c(NA, NA)
    }
  }
  id <-  unique(classification)
  id <- id[ id != (length(means)+1)]
  core_intervals <- t(sapply(id, function(j) {
    calc_core_interval(x[classification == j], core_prob)
  }))
  core_intervals <- as.data.frame(core_intervals)
  colnames(core_intervals) <- c("lower", "upper")
  row.names(core_intervals) <- means[id]

  # 3. 返回结果 -----------------------------------------------------------
  list(
    classification = classification,  # 每个点的类别（1:k_peaks为峰，k_peaks+1为背景）
    core_intervals = core_intervals  # 各峰核心区间

  )
}

setMethod("laplaceMM", "CNVisionObj", function(object,core_prob = 0.95,dist_type = "gaussian",background_type = "uniform") {
  #df <- do.call(rbind,object@result$ScaledCNV)
  #df <- df[df$seg.mean>= 0.5 & df$seg.mean <= 6.5,]
  df <- object@result$df
  model <- fixed_mean_mixture(x =  df$seg.mean,
                              weights = df$num.mark,
                              means = object@result$kde_peaks$peak_x,
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
    df2$num.mark <- df$num.mark

    return(df2)
  })

  stats <- lapply(object@result$ScaledCNV,function(df){
    weighted_mean <- sum(df$seg.mean * df$num.mark) / sum(df$num.mark)
    weighted_var <- sum(df$num.mark * (df$seg.mean - weighted_mean)^2) / sum(df$num.mark)
    weighted_sd <- sqrt(weighted_var)
    weighted_cv <- weighted_sd / weighted_mean

    return(weighted_cv)
  })
  CV <- apply(object@cnvData$gc.ratio, 2, function(x) {
    sd(x)/mean(x)
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

