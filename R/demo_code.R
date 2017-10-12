source("/Users/wayne/workspace/UsagePrediction/R/AnormalyDetection.R")
source("/Users/wayne/workspace/UsagePrediction/R/DetectAnoms.R")
source("/Users/wayne/workspace/UsagePrediction/R/util.R")
source("/Users/wayne/workspace/UsagePrediction/R/comb_forecast.R")
source("/Users/wayne/workspace/UsagePrediction/R/outlier_replace.R")
source("/Users/wayne/workspace/UsagePrediction/R/prediction_method.R")
source("/Users/wayne/workspace/UsagePrediction/R/main_function.R")

geo <- 'NA'
raw_data <- read.csv(paste0('/Users/wayne/workspace/PCG service parts/data/TopmostPN_filter_weekly_consumption(change)/TopmostPN_', geo, '_weekly_consumption_filter(change).csv'), header = TRUE)
raw_data2 <- read.csv(paste0('/Users/wayne/workspace/PCG service parts/data/TopmostPN_ib_filter(change)/TopmostPN_', geo, '_ib_withdate_filter.csv'), header = FALSE)
PN <- '00HT193'
PN_usage <- raw_data[raw_data$TopmostPN==PN,]
PN_usage$Data <- as.Date(PN_usage$Data)
PN_usage$DataNum <- as.numeric(as.vector(PN_usage$DataNum))
PN_usage <- PN_usage[, c('Data', 'DataNum')]
colnames(PN_usage) <- c('Date', 'Magnitude')
PN_IB <- raw_data2[raw_data2$V1==PN,]
PN_IB$V3 <- as.Date(PN_IB$V3)
PN_IB$V4 <- as.numeric(as.vector(PN_IB$V4))
PN_IB <- PN_IB[, c('V3', 'V4')]
colnames(PN_IB) <- c('Date', 'Magnitude')
results <- main_function(PN_usage, PN_IB, interval=26, percentage=1)
results$pred_plot <- results$pred_plot + ggplot2::geom_point(data=results$pred_results, ggplot2::aes_string(x="date", y="IB"), colour = "#009E73", alpha=0.6, size = 2)
results$pred_plot <- results$pred_plot + ggplot2::geom_line(data=results$pred_results, ggplot2::aes_string(x="date", y="IB"), colour = "#56B4E9", alpha=0.6)

