fourier_model = function (dat, method = c("auto", "manual"), non_sea_comp, drift = FALSE){
  model = list(aicc=Inf)
  method = match.arg(method)
  optim_K = 0
  for(i in 1:25)
  {
    if (method == "auto"){
      fit = auto.arima(dat, xreg=fourier(dat, K = i), seasonal = FALSE)
    }else fit = Arima(dat,order = non_sea_comp, xreg=fourier(dat, K = i), include.drift = drift)
    
    if(fit$aicc < model$aicc){
      model = fit
      optim_K = i
    }else break;
  }
  return (list(model=model, K=optim_K))
}
