prediction_method <- function(data, y, verbose = FALSE, 
                              lag = NULL, message = NULL, pred_length = 16, pred_method=c('holt.winter',
                                                                                          'stlf.arima',
                                                                                          'stlf.ets',
                                                                                          'seasonal.arima',
                                                                                          'NN.NNAR','IB')) {
  
  pckg = c("lubridate", "dplyr", "quadprog", "forecast")
  
  temp <- unlist(lapply(pckg, require, character.only=T))
  
  if (!all(temp==1) ) {
    stop("This function relies on packages \"dtwclust\" \"ggplot2\" \"lubridate\" \"dplyr\" \"quadprog\" \"zoo\".
         Use ?install.packages if they are not yet installed. \n")
  }
  
  if (!all(temp==1) ) {
    stop("This function relies on packages \"dtwclust\" \"ggplot2\" \"lubridate\" \"dplyr\" \"quadprog\" \"zoo\".
         Use ?install.packages if they are not yet installed. \n")
  }
  
  if(dim(data)[1] > dim(y)[1]*4.34524){
    stop("Data quality issue: Install base has less data than Usage")
  }
  
  if(data[[1L]][1] < y[[1L]][1]){
    stop("Data quality issue: the starting date of Usage is earlier than IB")
  }
  
  if(data[[1L]][dim(data)[1]] > y[[1L]][dim(y)[1]]){
    stop("Data quality issue: the ending date of Usage is later than IB")
  }
  
  if(!is.data.frame(data)){
    stop("data must be a single data frame.")
  } else {
    if(ncol(data) != 2 || !is.numeric(data[[2]])){
      stop("data must be a 2 column data.frame, with the first column being a set of timestamp, and the second coloumn being numeric values.")
    }
    if (!(class(data[[1]])[1] == "Date")) {
      stop("the first column being a set of timestamp.")
    }
  }
  
  if(!is.data.frame(y)){
    stop("data must be a single data frame.")
  } else {
    if(ncol(y) != 2 || !is.numeric(y[[2]])){
      stop("data must be a 2 column data.frame, with the first column being a set of timestamp, and the second coloumn being numeric values.")
    }
    if (!(class(y[[1]])[1] == "Date")) {
      stop("the first column being a set of timestamp.")
    }
  }
  
  num_obs <- nrow(data)
  
  if(num_obs > 13 * 2){
    num_obs_per_period = 13
  } else if(num_obs > 4 * 2){
    num_obs_per_period = 4
  } else if(num_obs < 4 * 2){
    stop("prediction needs at least 2 periods worth of data")
  }
  
  if(length(pred_method) != 1) {
    stop("Pick only one of the following:
         c(\"IB\", \"holt.winter\", \"stlf.arima\", \"stlf.ets\", \"seasonal.arima\", \"tslm.basic\", \"NN.NNAR\")")
  }
  
  #method IB
  if(pred_method=="IB"){
    lag_info <- lag
    length <- dim(subset(y, y[[1L]]>=data[[1L]][1] & y[[1L]]<=data[[1L]][nrow(data)]))[1]
    numbers <- dim(subset(y, y[[1L]]<data[[1L]][1]))[1] + 1
    index <- numbers - lag_info
    query <- reinterpolate(y[index:(length+index),][[2L]], dim(data)[1])
    reference <- data[[2L]]
    X2 <- cbind(query, 1)
    results <- solve.QP(t(X2) %*% X2, t(reference) %*% X2, cbind(c(min(y[[2L]]), 1), c(1, 0)), c(0, 0))
    y$baseline <- results$solution[2]+ results$solution[1]*y[[2L]]
    timestamp_lag <- y[[1]]
    month(timestamp_lag) <- month(timestamp_lag) + lag_info
    y_transfer <- data.frame(timestamp=timestamp_lag, baseline=y$baseline)
    
    predicted_df <- data.frame(Date=seq(data[dim(data)[1], ][[1]]%m+%weeks(1), by = "week", length.out = pred_length), Magnitude=rep(0, pred_length))
    
    decop <- stl(ts(data[[2L]], frequency = num_obs_per_period), s.window = "periodic", robust = TRUE)
    seasonal_req <- decop$time.series[,"seasonal"][1:num_obs_per_period]
    req_len <- dim(data)[1]+pred_length
    predicted_df$seaonal <- rep(seasonal_req, len=req_len)[(dim(data)[1]+1):req_len]
    y_trim <- y_transfer[y_transfer[[1]] >= predicted_df[1,][[1]], ]
    
    interpo <- as.Date(tapply(as.character(predicted_df[[1]]), format(predicted_df[[1]], "%Y-%m"), min))
    interpo_df <- data.frame(Date=interpo, Magnitude=rep(0, length(interpo)))
    for(j in 1:length(interpo)){
      for(i in 1:dim(y_trim)[1]){
        if(format(y_trim[[1]], "%Y-%m")[i]==format(interpo[j], "%Y-%m")){
          interpo_df[[2]][j] <- y_trim[[2]][i]
        }
      }
    }
    predicted_df[predicted_df[[1]] %in% interpo_df[[1]], ][[2]] <- interpo_df[[2]]
    replaced <- predicted_df[(predicted_df[[2]] == 0), ]
    
    predicted_df[[2]] <- replace_outlier(predicted_df[, c(1,2)], replaced[, c(1,2)], replace_method="linear")[[2]]
    predicted_df[[2]] <- predicted_df[[2]] + predicted_df[[3]]
    predicted_df <- predicted_df[, c(1,2)]
  }
  
  else if(pred_method == "holt.winter"){
    predicted_df <- data.frame(Date=seq(data[dim(data)[1], ][[1]]%m+%weeks(1), by = "week", length.out = pred_length), Magnitude=rep(0, pred_length))
    s <- ts(data[[2L]], frequency = num_obs_per_period)
    model <- HoltWinters(s)
    fc <- forecast(model, h=pred_length)
    predicted_df[[2]] <- fc$mean
    y_transfer <- data.frame()
  }
  
  else if(pred_method == "stlf.arima"){
    predicted_df <- data.frame(Date=seq(data[dim(data)[1], ][[1]]%m+%weeks(1), by = "week", length.out = pred_length), Magnitude=rep(0, pred_length))
    s <- ts(data[[2L]], frequency = num_obs_per_period)
    fc <- stlf(s, h=pred_length, s.window = "periodic", method='arima', ic='bic')
    predicted_df[[2]] <- fc$mean
    y_transfer <- data.frame()
  }
  
  else if(pred_method == "stlf.ets"){
    predicted_df <- data.frame(Date=seq(data[dim(data)[1], ][[1]]%m+%weeks(1), by = "week", length.out = pred_length), Magnitude=rep(0, pred_length))
    s <- ts(data[[2L]], frequency = num_obs_per_period)
    fc <- stlf(s, h=pred_length, s.window = "periodic", method='ets', ic='bic', opt.crit='mae')
    predicted_df[[2]] <- fc$mean
    y_transfer <- data.frame()
  }
  
  else if(pred_method == "seasonal.arima"){
    predicted_df <- data.frame(Date=seq(data[dim(data)[1], ][[1]]%m+%weeks(1), by = "week", length.out = pred_length), Magnitude=rep(0, pred_length))
    s <- ts(data[[2L]], frequency = num_obs_per_period)
    model <- auto.arima(s, ic='bic', seasonal.test = 'ch')
    fc <- forecast(model, h=pred_length)
    predicted_df[[2]] <- fc$mean
    y_transfer <- data.frame()
  }
  
  else if(pred_method == "NN.NNAR"){
    predicted_df <- data.frame(Date=seq(data[dim(data)[1], ][[1]]%m+%weeks(1), by = "week", length.out = pred_length), Magnitude=rep(0, pred_length))
    s <- ts(data[[2L]], frequency = num_obs_per_period)
    model <- nnetar(s)
    fc <- forecast(model, h=pred_length)
    predicted_df[[2]] <- fc$mean
    y_transfer <- data.frame()
  }
  
  return(list(pred = predicted_df, IB_trasfer = y_transfer)) 
  
}