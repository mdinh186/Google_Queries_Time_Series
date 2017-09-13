res_plot = function (model, lagM){
  res = residuals(model)
  par(mfrow= c(3,1))
  plot.ts(res, xlab = as.character(model))
  acf(res, lag.max = lagM)
  pacf(res, lag.max = lagM)
}
