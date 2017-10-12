detect_anoms <- function(data, y, k = 0.05, alpha = 0.1, num_obs_per_period = NULL,
                         one_tail = FALSE, upper_tail = TRUE, verbose = FALSE) {


  if(is.null(num_obs_per_period)){
    stop("must supply period length for time series decomposition")
  }

  num_obs <- nrow(data)

  if(num_obs < num_obs_per_period * 2){
    stop("Anom detection needs at least 2 periods worth of data")
  }

  if (length(rle(is.na(c(NA,data[[2L]],NA)))$values)>3){
    stop("Data contains non-leading NAs. We suggest replacing NAs with interpolated values (see na.approx in Zoo package).")
  } else {
    data <- na.omit(data)
  }

  data_decomp <- stl(ts(data[[2L]], frequency = num_obs_per_period),
                     s.window = "periodic", robust = TRUE)

  judge <- simi_measure(data, y)

  if(judge){

    lag_info <- lag_caculation(data, y)
    length <- dim(subset(y, y[[1L]]>=data[[1L]][1] & y[[1L]]<=data[[1L]][nrow(data)]))[1]
    numbers <- dim(subset(y, y[[1L]]<data[[1L]][1]))[1] + 1
    index <- numbers - lag_info
    query <- reinterpolate(y[index:(length+index),][[2L]], dim(data)[1])
    reference <- data[[2L]]
    X2 <- cbind(query, 1)
    results <- solve.QP(t(X2) %*% X2, t(reference) %*% X2, cbind(c(min(y[[2L]]), 1), c(1, 0)), c(0, 0))
    baseline <- results$solution[2] + results$solution[1]*query
    data <- data.frame(timestamp = data[[1L]], count = (data[[2L]]-data_decomp$time.series[,"seasonal"]-baseline))
  } else {

    length <- dim(subset(y, y[[1L]]>=data[[1L]][1] & y[[1L]]<=data[[1L]][nrow(data)]))[1]
    numbers <- dim(subset(y, y[[1L]]<data[[1L]][1]))[1] + 1
    query <- reinterpolate(y[numbers:(length+numbers),][[2L]], dim(data)[1])
    reference <- data[[2L]]
    X2 <- cbind(query, 1)
    results <- solve.QP(t(X2) %*% X2, t(reference) %*% X2, cbind(c(min(y[[2L]]), 1), c(1, 0)), c(0, 0))
    smoothed_baseline <- lowess(data_decomp$time.series[,"trend"], f=0.2)
    data <- data.frame(timestamp = data[[1L]], count = (data[[2L]]-data_decomp$time.series[,"seasonal"]-smoothed_baseline$y))
  }


  max_outliers <- trunc(num_obs*k)

  func_ma <- match.fun(median)
  func_sigma <- match.fun(mad)

  n <- length(data[[2L]])
  R_idx <- as.Date(data[[1L]][1L:max_outliers])
  num_anoms <- 0L

  for (i in 1L:max_outliers){

    if(one_tail){
      if(upper_tail){
        ares <- data[[2L]] - func_ma(data[[2L]])
      } else {
        ares <- func_ma(data[[2L]]) - data[[2L]]
      }
    } else {
      ares = abs(data[[2L]] - func_ma(data[[2L]]))
    }

    data_sigma <- func_sigma(data[[2L]])
    if(data_sigma == 0)
      break

    ares <- ares/data_sigma
    R <- max(ares)

    ttttt <- which(ares == R)[1L]

    R_idx[i] <- data[[1L]][ttttt]

    data <- data[-which(data[[1L]] == R_idx[i]), ]

    if(one_tail){
      p <- 1 - alpha/(n-i+1)
    } else {
      p <- 1 - alpha/(2*(n-i+1))
    }

    t <- qt(p,(n-i-1L))
    lam <- t*(n-i) / sqrt((n-i-1+t**2)*(n-i+1))

    if(R > lam)
      num_anoms <- i
  }

  if(num_anoms > 0) {
    R_idx <- R_idx[1L:num_anoms]
  } else {
    R_idx = NULL
  }

  if(judge){
    return (list(R_idx = R_idx, coefficients = results$solution, lag=lag_info, similarity=judge))
  } else {
    return (list(R_idx = R_idx, coefficients = results$solution, lag=NULL, similarity=judge))
  }
}
