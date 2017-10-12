get_range <- function(dfs, index = 2, y_log = F) {
  vals <- dfs[[index]]
  if(y_log) vals <- vals[vals > 0]
  vrange = range(vals, na.rm=TRUE)
  vmin = vrange[1]
  vmax = vrange[2]
  return(c(vmin, vmax))
}

savePlot <- function(myPlot) {
  pdf("myPlot.pdf")
  print(myPlot)
  dev.off()
}

simi_measure <- function(PN_usage, PN_IB){
  PN_IB <- PN_IB[index(PN_IB) >= ss.estimate(PN_IB),]
  PN_IB <- PN_IB[PN_IB[[1L]] <= PN_usage[[1L]][dim(PN_usage)[1]],]
  query <- reinterpolate(PN_IB[[2L]], new.length = round(dim(PN_IB)[1]*4.34524))
  query <- zscore(query)
  PN_usage <- PN_usage[index(PN_usage) >= ss.estimate(PN_usage),]
  PN_usage[[2L]] <- lowess(PN_usage[[2L]], f=0.2)$y
  reference <- reinterpolate(PN_usage[[2L]], new.length = length(query))
  reference <- zscore(reference)
  score <- lb_improved(reference, query, window.size = 10, norm = 'L2')
  if(score > 17) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

lag_caculation <- function(PN_usage, PN_IB) {

  numbers <- dim(subset(PN_IB, PN_IB[[1L]]<PN_usage[[1L]][1]))[1] + 1
  length <- dim(subset(PN_IB, PN_IB[[1L]]>=PN_usage[[1L]][1] & PN_IB[[1L]]<=PN_usage[[1L]][nrow(PN_usage)]))[1]
  d.lbi <- list()
  reference <- lowess(PN_usage[[2L]], f=0.2)$y
  reference <- zscore(reference)
  index <- 0
  min.lbi <- .Machine$double.xmax
  for(i in max(1,(numbers-12)):numbers){
    query <- reinterpolate(PN_IB[i:(length+i),][[2L]], dim(PN_usage)[1])
    query <- zscore(query)
    d.lbi <- lb_improved(reference, query, window.size = 10, norm = 'L2')
    if(d.lbi<min.lbi){
      min.lbi <- d.lbi
      index <- i
    }
  }
  lag <- numbers - index
  return(lag)
}

ss.estimate <- function(train){
  quantile.number <- quantile(train[[2L]], c(0.05, 0.1, 0.2, 0.25), na.rm = TRUE)
  n <- 1
  for(i in 3:nrow(train)){
    judge <- !((train[[2L]][i] < quantile.number[3] & train[[2L]][i-1] < quantile.number[3] &
                  train[[2L]][i-2] < quantile.number[3]) |
                 (train[[2L]][i] < quantile.number[4] & train[[2L]][i-1] < quantile.number[4]) |
                 (train[[2L]][i] < (mean(train[[2L]], na.rm = TRUE)/100) &
                    train[[2L]][i-1] < (mean(train[[2L]], na.rm = TRUE)/100) &
                    train[[2L]][i-2] < (mean(train[[2L]], na.rm = TRUE)/100))
    )
    if(judge){
      break
    }
    n <- i
  }
  n
}

frequency_measure <- function(PN_usage){
  
  if(mean(PN_usage[[2]])>5){
    message <- "high freq"
  } else {
    message <- "low freq"
  }
  return(message)
}


