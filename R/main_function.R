main_function <- function(PN_usage, PN_IB, percentage=1, interval=26, range=16) {
  
  if(1 < percentage || percentage <= 0){
    stop("Warning: percentage is the percent of usage data, and is between 0 and 1")
  }
  
  if(26 < interval || interval <= 0){
    stop("Warning: interval is the prediction range, and is between 0 and 26")
  }
  
  tt <- nrow(PN_usage)
  l <- round(percentage*tt)
  if(l <= (2*interval + range + 1)){
    stop("Warning: please incease the percentage, Pn_usage is not enough for prediction")
  }
  
  if(!is.data.frame(PN_usage)){
    stop("data must be a single data frame.")
  } else {
    if(ncol(PN_usage) != 2 || !is.numeric(PN_usage[[2]])){
      stop("data must be a 2 column data.frame, with the first column being a set of timestamp, and the second coloumn being numeric values.")
    }
    if (!(class(PN_usage[[1]])[1] == "Date")) {
      stop("the first column being a set of timestamp.")
    }
  }
  if (any((names(PN_usage) == c("Date", "Magnitude")) == FALSE)) {
    colnames(PN_usage) <- c("Date", "Magnitude")
  }
  
  if(!is.data.frame(PN_IB)){
    stop("data must be a single data frame.")
  } else {
    if(ncol(PN_IB) != 2 || !is.numeric(PN_IB[[2]])){
      stop("data must be a 2 column data.frame, with the first column being a set of timestamp, and the second coloumn being numeric values.")
    }
    if (!(class(PN_IB[[1]])[1] == "Date")) {
      stop("the first column being a set of timestamp.")
    }
  }
  if (any((names(PN_IB) == c("Date", "Magnitude")) == FALSE)) {
    colnames(PN_IB) <- c("Date", "Magnitude")
  }
  

  #performance plot and results
  PN_ts <- ts(PN_usage[1:l, ][[2]])
  
  me_accuracy_list <- c()
  IB_accuracy_list <- c()
  HW_accuracy_list <- c()
  sa_accuracy_list <- c()
  se_accuracy_list <- c()
  arima_accuracy_list <- c()
  nn_accuracy_list <- c()
  accuracy_list <- c()
  lag <- c()
  lag_list <- c()
  com_forecasts <- c()
  
  for(i in 1:interval){
    actual <- PN_ts[(i+range):(i+interval+range-1)]
    me_pred <-rep(mean(PN_ts[i:(i+range-1)]), interval)
    me_accuracy <- mean(abs(pmin(actual, me_pred)/pmax(actual, me_pred)) * 100)
    me_accuracy_list <- c(me_accuracy_list, me_accuracy)
    
    res1 <- Anomaly_Detection(PN_usage[1:(i+range-1), ], PN_IB, alpha = 0.1)
    rep_PN_usage <- replace_outlier(PN_usage[1:(i+range-1), ], res1$anoms, replace_method="karman")
    lag <- c(lag, lag_caculation(rep_PN_usage, PN_IB))
    lag_list <- c(lag_list, round(mean(lag)))
    pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=round(mean(lag)), message =res1$message, pred_length=26, pred_method='IB')
    IB_pred <- pred1$pred[[2]][1:interval]
    IB_accuracy <- mean(abs(pmin(actual, IB_pred)/pmax(actual, IB_pred)) * 100)
    IB_accuracy_list <- c(IB_accuracy_list, IB_accuracy)
    
    pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=res1$lag, message =res1$message, pred_length=interval, pred_method='holt.winter')
    HW_pred <- pred1$pred[[2]]
    HW_accuracy <- mean(abs(pmin(actual, HW_pred)/pmax(actual, HW_pred)) * 100)
    HW_accuracy_list <- c(HW_accuracy_list, HW_accuracy)
    
    pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=res1$lag, message =res1$message, pred_length=interval, pred_method='stlf.arima')
    sa_pred <- pred1$pred[[2]]
    sa_accuracy <- mean(abs(pmin(actual, sa_pred)/pmax(actual, sa_pred)) * 100)
    sa_accuracy_list <- c(sa_accuracy_list, sa_accuracy)
    
    pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=res1$lag, message =res1$message, pred_length=interval, pred_method='stlf.ets')
    se_pred <- pred1$pred[[2]]
    se_accuracy <- mean(abs(pmin(actual, se_pred)/pmax(actual, se_pred)) * 100)
    se_accuracy_list <- c(se_accuracy_list, se_accuracy)
    
    pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=res1$lag, message =res1$message, pred_length=interval, pred_method='seasonal.arima')
    arima_pred <- pred1$pred[[2]]
    arima_accuracy <- mean(abs(pmin(actual, arima_pred)/pmax(actual, arima_pred)) * 100)
    arima_accuracy_list <- c(arima_accuracy_list, arima_accuracy)
    
    pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=res1$lag, message =res1$message, pred_length=interval, pred_method='NN.NNAR')
    nn_pred <- pred1$pred[[2]]
    nn_accuracy <- mean(abs(pmin(actual, nn_pred)/pmax(actual, nn_pred)) * 100)
    nn_accuracy_list <- c(nn_accuracy_list, nn_accuracy)
    
    indiv_forecasts <- cbind(IB_pred[1], HW_pred[1], sa_pred[1], se_pred[1], arima_pred[1], nn_pred[1])
    com_forecasts <- rbind(com_forecasts, indiv_forecasts)
    
  }
  
  ind_nam <- c("IB", "HW", "sa", "se", "arima", "nn")
  colnames(com_forecasts) <- ind_nam
  actualed <- PN_ts[(1+range):(interval+range)]
  
  for(i in (interval+1):(length(PN_ts)-range-interval)){
    
    actual <- PN_ts[(i+range):(i+interval+range-1)]
    me_pred <-rep(mean(PN_ts[i:(i+range-1)]), interval)
    me_accuracy <- mean(abs(pmin(actual, me_pred)/pmax(actual, me_pred)) * 100)
    me_accuracy_list <- c(me_accuracy_list, me_accuracy)
    
    res1 <- Anomaly_Detection(PN_usage[1:(i+range-1), ], PN_IB, alpha = 0.1)
    rep_PN_usage <- replace_outlier(PN_usage[1:(i+range-1), ], res1$anoms, replace_method="karman")
    lag <- c(lag, lag_caculation(rep_PN_usage, PN_IB))
    pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=round(mean(lag)), message =res1$message, pred_length=26, pred_method='IB')
    IB_pred <- pred1$pred[[2]][1:interval]
    IB_accuracy <- mean(abs(pmin(actual, IB_pred)/pmax(actual, IB_pred)) * 100)
    IB_accuracy_list <- c(IB_accuracy_list, IB_accuracy)
    
    pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=res1$lag, message =res1$message, pred_length=interval, pred_method='holt.winter')
    HW_pred <- pred1$pred[[2]]
    HW_accuracy <- mean(abs(pmin(actual, HW_pred)/pmax(actual, HW_pred)) * 100)
    HW_accuracy_list <- c(HW_accuracy_list, HW_accuracy)
    
    pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=res1$lag, message =res1$message, pred_length=interval, pred_method='stlf.arima')
    sa_pred <- pred1$pred[[2]]
    sa_accuracy <- mean(abs(pmin(actual, sa_pred)/pmax(actual, sa_pred)) * 100)
    sa_accuracy_list <- c(sa_accuracy_list, sa_accuracy)
    
    
    pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=res1$lag, message =res1$message, pred_length=interval, pred_method='stlf.ets')
    se_pred <- pred1$pred[[2]]
    se_accuracy <- mean(abs(pmin(actual, se_pred)/pmax(actual, se_pred)) * 100)
    se_accuracy_list <- c(se_accuracy_list, se_accuracy)
    
    pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=res1$lag, message =res1$message, pred_length=interval, pred_method='seasonal.arima')
    arima_pred <- pred1$pred[[2]]
    arima_accuracy <- mean(abs(pmin(actual, arima_pred)/pmax(actual, arima_pred)) * 100)
    arima_accuracy_list <- c(arima_accuracy_list, arima_accuracy)
    
    pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=res1$lag, message =res1$message, pred_length=interval, pred_method='NN.NNAR')
    nn_pred <- pred1$pred[[2]]
    nn_accuracy <- mean(abs(pmin(actual, nn_pred)/pmax(actual, nn_pred)) * 100)
    nn_accuracy_list <- c(nn_accuracy_list, nn_accuracy)
    
    combine_f <- Forecast_comb(actualed, com_forecasts, Averaging_scheme = "cls")
    fhat <- cbind(IB_pred, HW_pred, sa_pred, se_pred, arima_pred, nn_pred)
    TT <- NROW(fhat)
    predicted <- t(combine_f$weights %*% t(fhat))
    accuracy <- mean(abs(pmin(actual, predicted)/pmax(actual, predicted)) * 100)
    accuracy_list <- c(accuracy_list, accuracy)
    
    indiv_forecasts <- cbind(IB_pred[1], HW_pred[1], sa_pred[1], se_pred[1], arima_pred[1], nn_pred[1])
    com_forecasts <- rbind(com_forecasts, indiv_forecasts)
    actualed <- c(actualed, actual[1])
  }
  
  # lll <- length(PN_ts)-range-interval

  me_accuracy_mean <- mean(me_accuracy_list, na.rm = TRUE)
  IB_accuracy_mean <- mean(IB_accuracy_list, na.rm = TRUE)
  HW_accuracy_mean <- mean(HW_accuracy_list, na.rm = TRUE)
  sa_accuracy_mean <- mean(sa_accuracy_list, na.rm = TRUE)
  se_accuracy_mean <- mean(se_accuracy_list, na.rm = TRUE)
  arima_accuracy_mean <- mean(arima_accuracy_list, na.rm = TRUE)
  nn_accuracy_mean <- mean(nn_accuracy_list, na.rm = TRUE)
  accuracy_mean <- mean(accuracy_list, na.rm = TRUE)
  accuracy_results <- c(me_accuracy_mean, IB_accuracy_mean, HW_accuracy_mean, sa_accuracy_mean, se_accuracy_mean, 
                        arima_accuracy_mean, nn_accuracy_mean, accuracy_mean)
  accuracy_results <- matrix(accuracy_results, nrow = 1)
  colnames(accuracy_results) <- c("MA" ,"IB", "HW", "sa", "se", "arima", "nn", "combined_f")
  
  weights <- matrix(combine_f$weights, nrow = 1)
  colnames(weights) <- c("IB", "HW", "sa", "se", "arima", "nn")
  weights_dt <- reshape2::melt(data.frame(weights))
  colnames(weights_dt) <- c("methods", "weights")
  weights_plot <- ggplot(weights_dt, aes(x="", y=weights, fill=methods)) + theme_bw() +
    geom_bar(aes(fill=methods), position="dodge", stat="identity")
  
  date_range1 <- PN_usage[[1]][(range+1):(length(PN_ts)-interval)]
  date_range2 <- PN_usage[[1]][(range+interval+1):(length(PN_ts)-interval)]
  accurary_dt <- data.frame(date=date_range1, IB=IB_accuracy_list, MA=me_accuracy_list)
  accurary_dt <- reshape2::melt(accurary_dt, id="date")
  colnames(accurary_dt)[3] <- "accuracy"
  accurary_dt <- rbind(accurary_dt, data.frame(date=date_range2, variable="combined_f", accuracy=accuracy_list))
  accuracy_plot <- ggplot(data=accurary_dt, aes(x=date, y=accuracy, colour=variable)) + theme_bw() + geom_line()
  
  
  #prediction resuls
  res1 <- Anomaly_Detection(PN_usage[1:l, ], PN_IB, alpha = 0.1)
  res1$plot <- res1$plot + ggplot2::geom_line(data=PN_usage, ggplot2::aes_string(x="Date", y="Magnitude"), alpha=0.2)
  if(!(dim(res1$anoms)[1]==0)){
    rep_PN_usage <- replace_outlier(PN_usage[1:l, ], res1$anoms, replace_method="karman")
    replace_data <- rep_PN_usage[(rep_PN_usage[[1]] %in% res1$anoms[[1]]), ]
    res1$plot <- res1$plot + ggplot2::geom_point(data=replace_data, ggplot2::aes_string(x="Date", y="Magnitude"), colour = "black", size = 3, shape = 1)
    
  } else {
    rep_PN_usage <- PN_usage[1:l, ]
  }
  res1$plot <- res1$plot + ggplot2::geom_point(data=rep_PN_usage, ggplot2::aes_string(x="Date", y="Magnitude"), colour = "blue", alpha=0.2)
  lag <- c(lag, lag_caculation(rep_PN_usage, PN_IB))
  pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=round(mean(lag)), message =res1$message, pred_length=26, pred_method='IB')
  IB_pred_results <- pred1$pred[1:interval, ]
  res1$plot <- res1$plot + ggplot2::geom_point(data=pred1$IB_trasfer, ggplot2::aes_string(x="timestamp", y="baseline"), alpha=0.6, size = 1)
  
  me_pred <-rep(mean(rep_PN_usage[(l-range+1):l, ][[2]]), interval)
  IB_pred <- pred1$pred[[2]][1:interval]
  
  pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=res1$lag, message =res1$message, pred_length=interval, pred_method='holt.winter')
  HW_pred <- pred1$pred[[2]]
  
  pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=res1$lag, message =res1$message, pred_length=interval, pred_method='stlf.arima')
  sa_pred <- pred1$pred[[2]]
  
  pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=res1$lag, message =res1$message, pred_length=interval, pred_method='stlf.ets')
  se_pred <- pred1$pred[[2]]
  
  pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=res1$lag, message =res1$message, pred_length=interval, pred_method='seasonal.arima')
  arima_pred <- pred1$pred[[2]]
  
  pred1 <- prediction_method(rep_PN_usage, PN_IB, lag=res1$lag, message =res1$message, pred_length=interval, pred_method='NN.NNAR')
  nn_pred <- pred1$pred[[2]]
  
  combine_f <- Forecast_comb(actualed, com_forecasts, Averaging_scheme = "cls")
  fhat <- cbind(IB_pred, HW_pred, sa_pred, se_pred, arima_pred, nn_pred)
  TT <- NROW(fhat)
  predicted <- t(combine_f$weights %*% t(fhat))
  
  pred_results <- data.frame(date=IB_pred_results[[1]], IB=IB_pred, combined_f=predicted, MA=me_pred)
  
  lag_plot <- ggplot(data.frame(lag), aes(lag)) + theme_bw() + geom_density()
  
  simi_message <- res1$message
  
  freq_message <- frequency_measure(PN_usage)
  
  return(list(accuracy_plot=accuracy_plot, weights_plot=weights_plot,weights=weights, accuracy_results=accuracy_results, 
              pred_results=pred_results, pred_plot=res1$plot, lag_plot=lag_plot, simi_message=simi_message, freq_message=freq_message))
  
}