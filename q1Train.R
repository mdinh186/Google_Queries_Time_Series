setwd("/Users/MyDinh/Downloads/Stat 153/googleTrend/Data")
library(caret)
library(tseries)
library(forecast)
library(stats)
library(car)



#==============================================================
#helper function

#plot ACF and PACF of the residuals
res_plot = function (model){
  res = residuals(model)
  lagM = 104
  par(mfrow= c(3,1))
  plot.ts(res, xlab = as.character(model))
  acf(res, lag.max = lagM)
  pacf(res, lag.max = lagM)
}

# since it's long seasonality, We use Fourier transformation instead.
# http://robjhyndman.com/hyndsight/longseasonality/
# K in the post is compute by choosing K that minimizes AICc.
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

# Step 1: Identify unusual observations and transform data

df1 = read.csv("q1_train.csv")
df1 = df1[,2]
df1 = ts(df1,freq = 365.25/7, start = 2004)

# detect and replace outlier.
df1_out = tsclean(df1)
#plot dataset:
plot(df1, type = "l") #doesn't look stationary
acf(df1)

#log transform by the min:
df1a = df1
df1_min = min(df1a)
df1a= log(df1a - df1_min + 0.001)
df1a_diff = diff(df1a)
tsdisplay(df1a_diff)
# doesn't seem to work
# try difference on the original data:
df1b = df1
df1b= diff(df1b)
tsdisplay(df1b)


#try difference on the log transform:
df1ab_diff = diff(df1a)
tsdisplay(df1ab_diff)# doesn't seem to change anything


#try Yeo-Johnson transformation
df1e = read.csv("q1_train.csv")
yw_param = preProcess(df1e, method = c("YeoJohnson"))
yw_param  #search for optimal lambda, return #0.14 but 0.15 is still within CI so choose 0.15 for convenience
trans = yjPower(df1, 0.15, jacobian.adjusted = FALSE) #stabilize variance
trans_diff = diff(trans) #stabilize mean
par(mfrow= c(3,1))
tsdisplay(trans_diff) #pacf: there are apparently 3 clusters at 24, 48, 52,
#There are peak at 52, 106, 156 in acf.

# it seems that yeo_johnson transformation didn't do well on outlier version 

# Now check if we have stationary process
# Check using Ljung- Box test:
Box.test(trans_diff) #p-value = 0.005694 indicates that the series is stationary

# Check periodic diagram to see if there is any unusual frequency
n = length(trans_diff)
I = abs(fft(trans_diff))^2/n #perioddic
P = (4/n)*I[1:(n/2)]
f = 0:round((n/2 -1))/n
plot(f, P,type = "l", xlab = "Frequency", ylab = "Scaled Periodogram")
#plot seems to have a lot of spike.


#look at monthy mean:
mdf1 = df1e
mdf1$Date = as.Date(mdf1$Date)
mdf1$month = format(mdf1$Date, "%m")
monthly_mean = aggregate(mdf1$activity, by = list(mdf1$month), FUN = mean)
plot(monthly_mean,type="b", main="Monthly Means Plot for Trend", xlab="Month", ylab="Mean")
abline(h = mean(monthly_mean$x))



# There maybe multiple seasonality patterns ( weekly, monthly): suggest to use tbat


#==============================================================
# Fit Model & Diagnostic.

#first initial guess using auto.arima
# may need to turn off approximation for long time series.


modelARIMA = auto.arima(trans, stepwise = FALSE, approximation = FALSE, max.D = 7)
res_plot(modelARIMA) 


#improve model using fourier transform
trants = ts(trans, freq = 365.25/7, start = 2004)
m1 = fourier_model(trants, method = "auto") 
res_plot(m1)
#ARIMA(5,1,3) with drift
# the acf of residual plot still show signficant spike at 1.
# Ljung-test to diagnosis it:
Box.test(res1, type = c("Ljung-Box")) # it seems to fails ljung-box test but it actually doens't matter.

#choose the lowest 3 AIC model.

#tuning model's parameter:
m2 = fourier_method(trants, c(4,1,3))#AIC=-112
res_plot(m2)
m3 = fourier_method(trants, c(3,1,3))#AIC=-101.91
res_plot(m3)
m4 = fourier_method(trants, c(2,1,3)) #AIC=-100.27
res_plot(m4)
#m5 = fourier_method(trants, c(5,1,3)) #AIC=-107.84
res_plot(m5)
m6 = fourier_method(trants, c(6,1,3)) #AIC=-113.5
res_plot(m6)
m7 = fourier_method(trants, c(5,1,2)) #AIC=-106.81
res_plot(m7)
m8 = fourier_method(trants, c(4,1,2)) #AIC=-104.94
res_plot(m8)
m9 = fourier_method(trants, c(3,1,2)) #AIC=-100.07
res_plot(m9)

# after checking plot, model of consideration: 7, 8, 2 (8 and 9 are probably the same)

# Using CV to get the error of those models.
# @Andy, can you implement it and check which of model 7,8,2 fit better/






# compare its statistics: training, test errors, AIC, BIC

#note: some model may have the smallest AIC but not RSME and vice versa
# when checking the model using ljung box test, sometimes the model doens't pass the test but
# it still gives the lowest RSME

#==============================================================
# Forecast
fc5 = forecast(m5, xreg=fourier(trants, K=8, h=104))
plot(fc5)

fc1 = forecast(m8, xreg=fourier(trants, K=8, h=104))
plot(fc8)

fc2 = forecast(m2, xreg=fourier(trants, K=8, h=104))
plot(fc2)

fsea = forecast(modelARIMA, h=104)
plot(fsea)
