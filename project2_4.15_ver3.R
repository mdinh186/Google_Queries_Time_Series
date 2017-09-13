setwd("/Users/MyDinh/Downloads/Stat 153/googleTrend")
library(caret)
library(tseries)
library(forecast)
library(stats)
library(car)
library(tsoutliers)



#==============================================================
#helper function

#plot ACF and PACF of the residuals
res_plot = function (model){
  res = residuals(model)
  par(mfrow= c(3,1))
  plot.ts(res, xlab = as.character(model))
  acf(res, lag.max = 104)
  pacf(res, lag.max = 104)
}


# since it's long seasonality, We use Fourier transformation instead.
# http://robjhyndman.com/hyndsight/longseasonality/
# K in the post is computed by choosing K that minimizes AICc.
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

#==============================================================


#try Yeo-Johnson transformation
df1e = read.csv("q1_train.csv")
yw_param = preProcess(df1e, method = c("YeoJohnson"))
yw_param  #search for optimal lambda, return #0.14 but 0.15 is still within CI so choose 0.15 for convenience
trans = yjPower(df1, 0.15, jacobian.adjusted = FALSE) #stabilize variance
trans_diff = diff(trans) #stabilize mean
par(mfrow= c(3,1))
tsdisplay(trans_diff) #pacf: there are apparently 3 clusters at 24, 48, 52,
#There are peak at 52, 106, 156 in acf.

# Using CV to get the error of those models.
model_list = list(c(3,1,0),c(3,1,1),c(3,1,2),c(3,1,3),c(4,1,0),c(4,1,1),c(4,1,2),c(4,1,3),c(5,1,0),c(5,1,1),c(5,1,2),c(5,1,3))
n=length(trants)
mse_list=c()

aic_calculator = function(model_list,drift=FALSE){
  out = list()
  index = 1
  for (param in model_list){
    temp = fourier_model(trants, method = c("manual"), param, drift=drift)
    out[[index]]=c(temp$model$aic, temp$model$bic)
    index = index + 1
  }
  return(out)
}

aic_calculator(model_list)

cv_mse_calculator = function(model_list, drift=TRUE){
  for (param in model_list){
    mse = 0
    for (t in 522:524){
      browser()
      train = ts(trants[1:t], freq = 365.25/7, start = 2004)
      test = trants[(t+1):n]
      train_fit = fourier_model(train, method = c("manual"), param, drift=drift)
      fc = forecast(train_fit$model, xreg=fourier(train, K=train_fit$K, h=n-t))
      error = as.numeric(fc$mean)-test
      mse = mse + sum(error^2)
      browser()
    }
    mse_list = c(mse_list, mse)
  }
  return(mse_list)
}
cv_mse_calculator(model_list)
cv_mse_calculator(list(c(3,1,0)), FALSE)


# compare its statistics: training, test errors, AIC, BIC

#note: some model may have the smallest AIC but not RSME and vice versa
# when checking the model using ljung box test, sometimes the model doens't pass the test but
# it still gives the lowest RSME

#==============================================================
# Forecast
fc = forecast(m7$model, xreg=fourier(trants, K=m7$K, h=104))
plot(fc)


#==============================================================
#==============================================================
#==============================================================

## Question2

df2 = read.csv("q2_train.csv")
df2 = df2[,2]


#plot dataset:
plot(df2, type = "l") #doesn't look stationary
acf(df2)
pacf(df2)

#log transform by the min:
df2a = df2
df2_min = min(df2a)
df2a= log(df2a - df2_min + 0.001)
plot.ts(df2a) # doesn't seem to work
acf(df2a)
pacf(df2a)

# try difference on the original data:
df2b = df2
df2b= diff(df2b)

df2b_min = min(df2b)
df2b= log(df2b - df2b_min + 0.001)

plot.ts(df2b)
acf(df2b)
pacf(df2b)

#try difference on the log transform:
df2c = diff(df2a)
plot.ts(df2c)
acf(df2c)
pacf(df2c)

#try to take log transform of differenced data
df2d = df2b
df2d_min = min(df2d)
df2d= log(df2d - df2d_min + 0.001)
plot.ts(df2d) #very good except around 460
acf(df2d) #marvelous
pacf(df2d) #marvelous

