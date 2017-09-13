model_list = list(c(3,1,2),c(3,1,3),c(4,1,3))
n=length(trants)
mse_list=c()
cv_mse_calculator = function(model_list){
  for (param in model_list){
    mse = 0
    for (t in floor(n/2):n-1){
      train = ts(trants[1:t], freq = 365.25/7, start = 2004)
      test = trants[(t+1):n]
      train_fit = fourier_model(train, method = c("manual"), param)
      fc = forecast(train_fit$model, xreg=fourier(train, K=train_fit$K, h=n-t))
      error = as.numeric(fc$mean)-test
      mse = mse + sum(error^2)
    }
    mse_list = c(mse_list, mse)
  }
  return(mse_list)
}
