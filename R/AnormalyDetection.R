Anomaly_Detection <- function(x, y, max_anoms = 0.05, direction = 'both',
                               alpha = 0.05, threshold = 'None', plot = TRUE,
                               xlabel = 'time', ylabel = 'usage',
                               title = NULL, verbose=FALSE, na.rm = FALSE){
  
  pckg = c("dtwclust", "ggplot2", "lubridate", "dplyr", "quadprog", "zoo")
  
  temp <- unlist(lapply(pckg, require, character.only=T))
  
  if (!all(temp==1) ) {
    stop("This function relies on packages \"dtwclust\" \"ggplot2\" \"lubridate\" \"dplyr\" \"quadprog\" \"zoo\".
         Use ?install.packages if they are not yet installed. \n")
  }

  if(dim(x)[1] > dim(y)[1]*4.34524){
    stop("Data quality issue: Install base has less data than Usage")
  }

  if(x[[1L]][1] < y[[1L]][1]){
    stop("Data quality issue: the starting date of Usage is earlier than IB")
  }

  if(x[[1L]][dim(x)[1]] > y[[1L]][dim(y)[1]]){
    stop("Data quality issue: the ending date of Usage is later than IB")
  }

  if(!is.data.frame(x)){
    stop("data must be a single data frame.")
  } else {
    if(ncol(x) != 2 || !is.numeric(x[[2]])){
      stop("data must be a 2 column data.frame, with the first column being a set of timestamp, and the second coloumn being numeric values.")
    }
    if (!(class(x[[1]])[1] == "Date")) {
      stop("the first column being a set of timestamp.")
    }
  }
  if (any((names(x) == c("timestamp", "count")) == FALSE)) {
    colnames(x) <- c("timestamp", "count")
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
  if (any((names(y) == c("timestamp", "count")) == FALSE)) {
    colnames(y) <- c("timestamp", "count")
  }


  if(!is.logical(na.rm)){
    stop("na.rm must be either TRUE (T) or FALSE (F)")
  }

  if(any(is.na(x$timestamp))){
    if(na.rm){
      x <- x[-which(is.na(x$timestamp)), ]
    } else {
      stop("timestamp contains NAs, please set na.rm to TRUE or remove the NAs manually.")
    }
  }

  if(any(is.na(y$timestamp))){
    if(na.rm){
      y <- y[-which(is.na(y$timestamp)), ]
    } else {
      stop("timestamp contains NAs, please set na.rm to TRUE or remove the NAs manually.")
    }
  }

  if(max_anoms > .49){
    stop(paste("max_anoms must be less than 50% of the data points (max_anoms =", round(max_anoms*length(x[[2]]), 0), " data_points =", length(x[[2]]),")."))
  } else if(max_anoms < 0){
    stop("max_anoms must be positive.")
  } else if(max_anoms == 0){
    warning("0 max_anoms results in max_outliers being 0.")
  }

  if(!direction %in% c('pos', 'neg', 'both')){
    stop("direction options are: pos | neg | both.")
  }
  if(!(0.01 <= alpha || alpha <= 0.1)){
    if(verbose) message("Warning: alpha is the statistical signifigance, and is usually between 0.01 and 0.1")
  }
  if(!threshold %in% c('None','med_max','p95','p99')){
    stop("threshold options are: None | med_max | p95 | p99.")
  }
  if(!is.logical(plot)){
    stop("plot must be either TRUE (T) or FALSE (F)")
  }
  if(!is.character(xlabel)){
    stop("xlabel must be a string")
  }
  if(!is.character(ylabel)){
    stop("ylabel must be a string")
  }
  if(!is.character(title) && !is.null(title)){
    stop("title must be a string")
  }
  if(is.null(title)){
    title <- ""
  } else {
    title <- paste(title, " : ", sep="")
  }

  period = 4
  num_obs <- length(x[[2]])

  if(max_anoms < 1/num_obs){
    max_anoms <- 1/num_obs
  }

  all_data <- list(x)
  all_anoms <- data.frame(timestamp=numeric(0), count=numeric(0))

  anomaly_direction = switch(direction,
                             "pos" = data.frame(one_tail=TRUE, upper_tail=TRUE), # upper-tail only (positive going anomalies)
                             "neg" = data.frame(one_tail=TRUE, upper_tail=FALSE), # lower-tail only (negative going anomalies)
                             "both" = data.frame(one_tail=FALSE, upper_tail=TRUE)) # Both tails. Tail direction is not actually used.

  s_h_esd_info <- detect_anoms(all_data[[1]], y, k=max_anoms, alpha=alpha, num_obs_per_period=period,
                                     one_tail=anomaly_direction$one_tail, upper_tail=anomaly_direction$upper_tail, verbose=verbose)

  s_h_esd_timestamps <- s_h_esd_info$R_idx
  coefficients <- s_h_esd_info$coefficients
  y$baseline <- coefficients[2]+ coefficients[1]*y[[2L]]
  timestamp_lag <- y$timestamp
  if(s_h_esd_info$similarity){
    month(timestamp_lag) <- month(timestamp_lag) + s_h_esd_info$lag
  }
  y_transfer <- data.frame(timestamp=timestamp_lag, baseline=y$baseline)


  if(!is.null(s_h_esd_timestamps)){
    anoms <- subset(all_data[[1]], (all_data[[1]][[1]] %in% s_h_esd_timestamps))
  } else {
    anoms <- data.frame(timestamp=numeric(0), count=numeric(0))
  }

  if(threshold != "None"){

    periodic_maxs <- tapply(x[[2]],as.Date(x[[1]]),FUN=max)

    if(threshold == 'med_max'){
      thresh <- median(periodic_maxs)
    }else if (threshold == 'p95'){
      thresh <- quantile(periodic_maxs, .95)
    }else if (threshold == 'p99'){
      thresh <- quantile(periodic_maxs, .99)
    }

    anoms <- subset(anoms, anoms[[2]] >= thresh)
  }
  all_anoms <- rbind(all_anoms, anoms)

  anom_pct <- (length(all_anoms[[2]]) / num_obs) * 100

  if(s_h_esd_info$similarity){
    message <- "trend match"
  } else {
    message <- "trend not match"
  }

  if(anom_pct == 0){
    if(verbose) message("No anomalies detected.")

    plot_title <-  title
    color_name <- paste("\"", title, "\"", sep="")
    xgraph <- ggplot2::ggplot(x, ggplot2::aes_string(x="timestamp", y="count")) + ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "gray60"), panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), text=ggplot2::element_text(size = 14))
    xgraph <- xgraph + ggplot2::geom_line(data=x, ggplot2::aes_string(colour=color_name), alpha=0.8)
    xgraph <- xgraph + ggplot2::geom_point(data=y, ggplot2::aes_string(x="timestamp", y="baseline"), alpha=0.2, size = 1)
    xgraph <- xgraph + ggplot2::geom_point(data=y_transfer, ggplot2::aes_string(x="timestamp", y="baseline"), alpha=0.4, size = 1)
    xgraph <- xgraph + ggplot2::labs(x=xlabel, y=ylabel)
    xgraph <- xgraph + ggplot2::theme(legend.position="none")

    anoms <- data.frame()
  } else {

    plot_title <-  title
    color_name <- paste("\"", title, "\"", sep="")
    xgraph <- ggplot2::ggplot(x, ggplot2::aes_string(x="timestamp", y="count")) + ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "gray60"), panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), text=ggplot2::element_text(size = 14))
    xgraph <- xgraph + ggplot2::geom_line(data=x, ggplot2::aes_string(colour=color_name), alpha=0.8)
    xgraph <- xgraph + ggplot2::geom_point(data=y, ggplot2::aes_string(x="timestamp", y="baseline"), alpha=0.1, size = 1)
    xgraph <- xgraph + ggplot2::geom_point(data=y_transfer, ggplot2::aes_string(x="timestamp", y="baseline"), alpha=0.3, size = 1)
    xgraph <- xgraph + ggplot2::labs(x=xlabel, y=ylabel, title=plot_title)
    xgraph <- xgraph + ggplot2::geom_point(data=all_anoms, ggplot2::aes_string(color=paste("\"zzz_",title,"\"",sep="")), size = 3, shape = 1)
    xgraph <- xgraph + ggplot2::theme(legend.position="none")

    all_anoms[[1]] <- as.Date(format(all_anoms[[1]], format="%Y-%m-%d"))
    anoms <- data.frame(timestamp=all_anoms[[1]], anoms=all_anoms[[2]], stringsAsFactors=FALSE)
  }


  if(plot){
    return (list(anoms = anoms, plot = xgraph, message = message, lag = s_h_esd_info$lag))
  } else {
    return (list(anoms = anoms, plot = plot.new(), message = message, lag = s_h_esd_info$lag))
  }
}
