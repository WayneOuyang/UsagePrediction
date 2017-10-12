replace_outlier <- function(data, anoms, replace_method=c("linear", "spline", "stine", "karman")) {
  
  #no anoms
  if(dim(anoms)[1]==0) {
    return(data)
  }
  
  pckg = c("stats", "stinepack", "forecast")
  temp <- unlist(lapply(pckg, require, character.only=T))
  
  if (!all(temp==1) ) {
    stop("This function relies on packages \"dtwclust\" \"ggplot2\" \"lubridate\" \"dplyr\" \"quadprog\" \"zoo\".
         Use ?install.packages if they are not yet installed. \n")
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
  
  if(!is.data.frame(anoms)){
    stop("data must be a single data frame.")
  } else {
    if(ncol(anoms) != 2 || !is.numeric(anoms[[2]])){
      stop("data must be a 2 column data.frame, with the first column being a set of timestamp, and the second coloumn being numeric values.")
    }
    if (!(class(anoms[[1]])[1] == "Date")) {
      stop("the first column being a set of timestamp.")
    }
  }
  
  
  if(length(replace_method) != 1) {
    stop("Pick only one of the following:
         c(\"linear\", \"spline\", \"stine\", \"karman\")")
  }
  
  n <- dim(data)[1]
  
  allindx <- 1:n
  data$index <- 1:n
  
  indx <- data[!(data[[1]] %in% anoms[[1]]), ][[3]]
  missindx <- data[(data[[1]] %in% anoms[[1]]), ][[3]]
  
  data.vec <- as.vector(data[[2]])
  
  if(replace_method =="linear") {
    interp <- approx(indx,data.vec[indx], allindx, rule=2)$y
    data[data[[3]] %in% missindx, ][[2]] <- interp[missindx]
    data <- data[, c(1, 2)]
  }
  
  else if(replace_method == "spline") {
    interp <- spline(indx,data.vec[indx],n = n)$y
    data[data[[3]] %in% missindx, ][[2]] <- interp[missindx]
    data <- data[, c(1, 2)]
  }
  
  else if(replace_method == "stine") {
    interp <- stinterp(indx,data.vec[indx], allindx)$y
    data[data[[3]] %in% missindx, ][[2]] <- interp[missindx]
    data <- data[, c(1, 2)]
  }
  
  else if(replace_method == "karman") {
    nit=-1
    mod <- StructTS(data.vec)$model0
    kal <- KalmanSmooth(data.vec, mod, nit )
    erg <- kal$smooth
    karima <- erg[missindx, ,drop=FALSE] %*% as.matrix(mod$Z)
    data[data[[3]] %in% missindx, ][[2]] <- karima
    data <- data[, c(1, 2)]
  }
  
  return(data)
}