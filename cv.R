
cv_ts = function(k)
 # minimum data length for fitting a model
n <- length(a10)
mae1 <- mae2 <- mae3 <- matrix(NA,12,12)
st <- tsp(a10)[1]+(k-1)/12
for(i in 1:12)
{
  xshort <- window(a10, end=st + (i-1))
  xnext <- window(a10, start=st + (i-1) + 1/12, end=st + i)
  fit1 <- tslm(xshort ~ trend + season, lambda=0)
  fcast1 <- forecast(fit1, h=12)
  fit2 <- Arima(xshort, order=c(3,0,1), seasonal=list(order=c(0,1,1), period=12), 
                include.drift=TRUE, lambda=0, method="ML")
  fcast2 <- forecast(fit2, h=12)
  fit3 <- ets(xshort,model="MMM",damped=TRUE)
  fcast3 <- forecast(fit3, h=12)
  browser()
  mae1[i,] <- abs(fcast1[['mean']]-xnext)
  mae2[i,] <- abs(fcast2[['mean']]-xnext)
  mae3[i,] <- abs(fcast3[['mean']]-xnext)
}
cv_ts(60)
