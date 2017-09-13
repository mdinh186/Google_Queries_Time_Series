#setwd("C:/Users/Andy Hyunouk Ko/Desktop/School/Probability&Statistics/Stat153/project")
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
adf.test(trans_diff, alternative = "stationary") #p-value = 0.01(small) indicates stationarity
kpss.test(trans_diff) #null hypothesis is that data is stationary. p-value=0.1 --> data is stationary for alpha=0.05

# Check periodic diagram to see if there is any unusual frequency
n = length(trans_diff)
I = abs(fft(trans_diff))^2/n #perioddic
P = (4/n)*I[1:(n/2)]
f = 0:round((n/2 -1))/n
plot(f, P,type = "l", xlab = "Frequency", ylab = "Scaled Periodogram")
#plot seems to have a lot of spikes.


#look at monthy mean:
mdf1 = df1e
mdf1$Date = as.Date(mdf1$Date)
mdf1$month = format(mdf1$Date, "%m")
monthly_mean = aggregate(mdf1$activity, by = list(mdf1$month), FUN = mean)
plot(monthly_mean,type="b", main="Monthly Means Plot for Trend", xlab="Month", ylab="Mean")
abline(h = mean(monthly_mean$x))



# There maybe multiple seasonality patterns ( weekly, monthly): suggest to use that


#==============================================================
# Fit Model & Diagnostic.

#first initial guess using auto.arima
# may need to turn off approximation for long time series.
trants = ts(trans, freq = 365.25/7, start = 2004)

modelARIMA = auto.arima(trants, stepwise = FALSE, approximation = FALSE, D=1)
res_plot(modelARIMA) #AIC=108.38 
#ARIMA(0,1,3)(0,1,2)[52]     

#improve model using fourier transform
m1 = fourier_model(trants, method = "auto") #AIC = 67.5
res_plot(m1$model)
# the acf of residual plot still show signficant spike at 1.
# Ljung-test to diagnosis it:
Box.test(residuals(m1$model), type = c("Ljung-Box")) # it seems to fails ljung-box test but it actually doens't matter.

#choose the lowest 3 AIC model.

#tuning model's parameter:
m2 = fourier_model(trants, method = c("manual"), c(3,1,1), TRUE)#AIC=-11.81 
res_plot(m2)
m3 = fourier_model(trants, method = c("manual"), c(4,1,1), TRUE)#AIC=-10.13
res_plot(m3)
m10 = fourier_model(trants, method = c("manual"), c(5,1,1), T) #AIC=-11.24 
res_plot(m10)
m4 = fourier_model(trants, method = c("manual"), c(3,1,1)) #AIC=-10.4
res_plot(m4)
m5 = fourier_model(trants, method = c("manual"), c(4,1,1)) #AIC=-8.48 
res_plot(m5)
m12 = fourier_model(trants, method = c("manual"), c(5,1,1)) #AIC=-8.53 
res_plot(m12)
m6 = fourier_model(trants, method = c("manual"), c(3,1,2), T) #AIC=-10.84
res_plot(m6)
m7 = fourier_model(trants, method = c("manual"), c(4,1,2),T) #AIC=-7.93
res_plot(m7)
m11 = fourier_model(trants, method = c("manual"), c(5,1,2), T) #AIC=-9.33
res_plot(m11)
m8 = fourier_model(trants, method = c("manual"), c(3,1,2)) #AIC=-9.49 
res_plot(m8)
m9 = fourier_model(trants, method = c("manual"), c(4,1,2)) #AIC=-6.62 
res_plot(m9)
m16 = fourier_model(trants, method = c("manual"), c(5,1,2)) #AIC=-6.73  
res_plot(m16)
m15 = fourier_model(trants, method = c("manual"), c(3,1,3),T) #AIC=-12.62
res_plot(m15)
m13 = fourier_model(trants, method = c("manual"), c(4,1,3),T) #AIC=-13.86
res_plot(m13)
m14 = fourier_model(trants, method = c("manual"), c(5,1,3),T) #AIC=-8.73 
res_plot(m14)
m17 = fourier_model(trants, method = c("manual"), c(3,1,3)) #AIC=-17.26   
res_plot(m17)
m18 = fourier_model(trants, method = c("manual"), c(4,1,3)) #AIC=-12.59
res_plot(m18)
m19 = fourier_model(trants, method = c("manual"), c(5,1,3)) #AIC=-10.75
res_plot(m19)

# after checking plot, model of consideration: 7, 8, 2 (8 and 9 are probably the same)

# Using CV to get the error of those models.
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

