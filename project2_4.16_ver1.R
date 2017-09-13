setwd("/Users/MyDinh/Downloads/Stat 153/googleTrend/Data")
library(caret)
library(tseries)
library(forecast)
library(stats)
library(car)
library(VGAM)
library(MASS)


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
df1 = ts(df1, freq = 365.25/7, start = 2004)

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
# after checking plot, model of consideration: 7, 8, 2 (8 and 9 are probably the same)

# Using CV to get the error of those models.
model_list = list(c(3,1,0),c(3,1,1),c(3,1,2),c(3,1,3),c(4,1,0),c(4,1,1),c(4,1,2),c(4,1,3),c(5,1,0),c(5,1,1),c(5,1,2),c(5,1,3))
n=length(trants)
mse_list=c()

information_calculator = function(model_list,drift=TRUE, ){
  out = list()
  index = 1
  for (param in model_list){
    temp = fourier_model(trants, method = c("manual"), param, drift=drift)
    out[[index]]=c(temp$model$aic, temp$model$bic)
    index = index + 1
  }
  return(out)
}



model_generator = function(ns_order){
  model_list=list()
  for (ar in (ns_order[1]-1):(ns_order[1]+1)){
    if (ar>=0){
      for (ma in (ns_order[3]-1):(ns_order[3]+1)){
        if (ma>=0){
          model_list[[length(model_list)+1]] = c(ar, ns_order[2], ma)
        }
      }
    }
  }
  return(model_list)
}

cross_validate = function(drift=TRUE, data, trants, model_list, inv, ...){
  mse_list = NULL
  for (param in model_list){
    mse = 0
    min_st = 161
    n = length(trants)
    fold = (n-161)/7
    for (i in 1:6){
      train_cut = min_st + (fold*(i))
      train = ts(trants[1:train_cut], freq = 365.25/7, start = 2004)
      test_cut = (train_cut+1):(train_cut+fold)
      test = data[test_cut]
      train_fit = fourier_model(train, method = c("manual"), param, drift=drift)
      fc = forecast(train_fit$model, xreg=fourier(train, K=train_fit$K, h=length(test_cut)))
      error = inv(as.numeric(fc$mean), ...) - test
      mse = mse + mean(error^2)
    }
    mse_list = c(mse_list, mse)
  }
  return(mse_list)
}


cross_validate2 = function(drift=TRUE, data, trants, model_list, inv, ...){
  mse_list = NULL
  for (param in model_list){
    mse = 0
    min_st = 161
    n = length(trants)
    fold = (n-161)/7
    for (i in 1:6){
      train_cut = min_st + (fold*(i))
      train = ts(trants[1:train_cut], freq = 365.25/7, start = 2004)
      test_cut = (train_cut+1):(train_cut+fold)
      test = data[test_cut]
      test_start = 2004 + 1/(365.25/7)*(train_cut+1)
      npred = length(test_cut) # number of periods ahead to forecast
      train_fit1 = fourier_model(train, method = c("manual"), param, drift=drift)
      ns_order = train_fit1$model$arma[c(1,6,2)]
      p = ns_order[1]; q = ns_order[3];  
      k = train_fit1$K
      mo = tso(train, types=c("AO","LS","TC"), tsmethod = "auto.arima", remove.method = "bottom-up",args.tsmethod = list(max.p=p, max.q=q, seasonal=FALSE, xreg=fourier(train, K=k)))
      newxreg = outliers.effects(mo$outliers, length(train)+ npred)
      #newxreg = ts(newxreg[-seq_along(train),],start = test_start, frequency = 365.25/7)
      fc = forecast.Arima(mo$fit, xreg= newxreg, h = npred)
      error = inv(as.numeric(fc$mean), ...) - test
      mse = mse + mean(error^2)
      }
    mse_list = c(mse_list, mse)
  }
  return(mse_list)
}



ultimate = function(drift=TRUE, data, trans, inv, ...){
  # model_list is a list of sets of parameters (a,b,c) for the non-seasonal part of ARIMA model
  # trans is the transformation to be applied to data. It must have parameter named "dat" which isthe original data.
  # inv is the inverse of the transformation applied to original data. It must have parameter named "data" which is the transformed data to be inverted back.
  # ... arguments for trans and inv functions separated by delimiter '..'
  
  # Deal with arguments to be passed into trans and inv
  arguments <- list(...)
  if('..' %in% arguments){
    if(arguments[1] == '..'){
      inv_arg <- arguments[2:length(arguments)]
      trans_arg = NULL
    } else if( arguments[length(arguments)]=='..'){
      trans_arg <- arguments[1:(length(arguments)-1)]
      inv_arg = NULL
    } else{
      i = which(arguments == '..')
      trans_arg <- arguments[1:(i-1)]
      inv_arg <- arguments[(i+1):length(arguments)]
    }
  } else error("You must include '..' as an argument")
  
  # transform data
  if (is.null(trans_arg)){
    trants = trans(data)
  } else trants = do.call(trans, append(list(data), trans_arg))
  #n = length(trants)
  
  # inital guess or parameters
  init_model = fourier_model(trants, method = "auto", drift=drift)
  
  ns_order = init_model$model$arma[c(1,6,2)] #non-seaonsal order
  p = ns_order[1]; q = ns_order[3]
  
  trants = transform(df4) #temporarily define to test outlier, 
  init_model2 = tso(trants, tsmethod = "auto.arima", remove.method = "bottom-up",args.tsmethod = list(max.p=p, max.q=q, seasonal=FALSE, xreg=fourier(trants, K=init_model$K)))
  
  
  # generate model list
  model_list = model_generator(ns_order)
  model_list2 = list(init_model2$fit$arma[c(1,6,2)])
  
  # implement cross-validation and choose&fit final model
  if (is.null(inv_arg)){
    mse_list = cross_validate(drift=TRUE, trants, model_list, inv)
  } else{
    mse_list = do.call(cross_validate, append(list(drift=TRUE, data=data, trants=trants, model_list=model_list, inv=inv), inv_arg))
    mse_list2 = do.call(cross_validate2, append(list(drift=TRUE, data=data, trants=trants, model_list=model_list2, inv=inv), inv_arg))
  } 
   
  
  #implement 2nd cross-validation for outlier version: 
  print(paste0("MSE with outlier effect is", mse_list2[1]))
  
  
  # implement final prediction with model with best MSE result
  final_model_param = unlist(model_list[mse_list == min(mse_list)])
  print(paste0("Non-seasonal model parameter used is ","(", toString(final_model_param),")"))
  print(paste0("MSE is ", mse_list[mse_list == min(mse_list)]))
  final_model = fourier_model(trants, method = "manual", final_model_param, drift=TRUE)
  print(paste0("Number of coefficients in Fourier transformation is ", final_model$K))
  
  #predict and plot
  fc = forecast(final_model$model, xreg=fourier(trants, K=final_model$K, h=104))
  if (is.null(inv_arg)){
    final_fitted = inv(fc$mean)
  } else final_fitted = do.call(inv, append(list(fc$mean), inv_arg))
  whole_data = ts(c(df1,final_fitted), freq = 365.25/7, start = 2004)
  robj = list(whole_data, final_model, mse=min(mse_list))
  
  return(robj)
}



# compare its statistics: training, test errors, AIC, BIC

#note: some model may have the smallest AIC but not RSME and vice versa
# when checking the model using ljung box test, sometimes the model doens't pass the test but
# it still gives the lowest RSME

#==============================================================
# Different set of transformations

transform = function(dat){
  bool = (class(dat) == 'ts')
  if(bool){
    tinfo = tsp(dat)
  }
  dat = as.numeric(dat)
  dat_min = min(dat)
  if (dat_min < 0){
    dat_pre = dat - dat_min + 0.001
  }
  lam= BoxCox.lambda(dat_pre, method = "loglik")
  out <- boxcox(dat_pre~1)
  x = range(out$x[out$y > max(out$y)-qchisq(0.95,1)/2])
  print(paste0("Optimal lambda is ",lam ))
  print(paste0("CI of lambda is ",x ))
  dfbc = BoxCox(dat_pre, lam)
  if (bool){
    dfbc = ts(dfbc,freq = tinfo[3], start = tinfo[1])
  }
  return (dfbc)
}

test = transform(df1)

inversebc = function(dat, lam, dat2){
  bool = (class(dat) == 'ts')
  if(bool){
    tinfo = tsp(dat)
  }
  dat = as.numeric(dat)
  dat2 = as.numeric(dat2)
  dfmin = min(dat2)
  inv = InvBoxCox(dat, lambda = lam)
  if (dfmin <0){
    inv = inv + dfmin -0.001
  }
  if(bool){
    inv = ts(inv, freq = tinfo[3], start = tinfo[1])
  }
  return (inv)
}

#calculate confidence interval to improve rounding
ci_yj = function(dat){
  out = boxCox(dat~1, family="yjPower", plotit = F)
  x = range(out$x[out$y > max(out$y)-qchisq(0.95,1)/2])
  print(x)
}
lam1 = ci_yj(df1)
t1 = transform(df1)
ult1 = ultimate(drift=TRUE, df1, yeo.johnson, yeo.johnson, lambda=0.1,'..',lambda=0.1, inverse=TRUE)
#MSE is 1.52269650593397
tult=ultimate(drift=TRUE, df1, transform, inversebc, "..", lam=0.25, dat2=df1)
"MSE is 1.56040151019266"
m1 = list(c(3, 1, 1), c(1, 1, 1))
#information_calculator(m1)
#==============================================================
#==============================================================
#==============================================================

## Question2

df2 = read.csv("q2_train.csv")
df2 = df2[,2]
df2 = ts(df2,freq = 365.25/7, start = 2004)


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
pacf(df2d) #marvelous'

#yeo_johnson transformation
trants2 = yeo.johnson(df2, 0.35)
trants2_diff = diff(trants2)
trants2_diff2 = diff(trants2_diff)

# Stationarity Check
Box.test(trants2_diff) #p-value = 0.005694 indicates that the series is stationary
adf.test(trants2_diff, alternative = "stationary") #p-value = 0.01(small) indicates stationarity
kpss.test(trants2_diff) #null hypothesis is that data is stationary. p-value=0.1 --> data is stationary for alpha=0.05

Box.test(trants2_diff2) #p-value = 0.005694 indicates that the series is stationary
adf.test(trants2_diff2, alternative = "stationary") #p-value = 0.01(small) indicates stationarity
kpss.test(trants2_diff2)

# Predict using the same procedure as before
t2 = transform(df2)
lam2 = ci_yj(df2)
ult2=ultimate(drift=TRUE, df2, yeo.johnson, yeo.johnson, lambda=0.4, "..", lambda=0.4, inverse=TRUE)
#"MSE is 6.78455037268759"
tult2 =ultimate(drift=TRUE, df2, transform, inversebc, "..", lam=0.3, dat2=df2)
#"MSE is 6.95680854334869"
#check AIC:
m2 = list(c(0, 1, 1), c(2, 1, 0))
information_calculator(m2)

#==============================================================
#==============================================================
#==============================================================

## Question 3
df3 = read.csv("q3_train.csv")
df3 = df3[,2]
df3 = ts(df3, freq = 365.25/7, start = 2004)


t3 = transform(df3)
lam3 = ci_yj(df3)
ult3 = ultimate(drift=TRUE, df3, yeo.johnson, yeo.johnson, lambda=-0.2, "..", lambda=-0.2, inverse=TRUE)
"MSE is MSE is 4.01375229263308"
tult3 = tult=ultimate(drift=TRUE, df3, transform, inversebc, "..", lam=0.2, dat2=df3)
#MSE is 4.05510717685506

#check AIC: 
m3 = list(c(4, 1, 0), c(2, 1, 0))
information_calculator(m3)

#==============================================================
#==============================================================
#==============================================================

## Question 4
df4 = read.csv("q4_train.csv")
df4 = df4[,2]
df4 = ts(df4,freq = 365.25/7, start = 2004)

t4 = transform(df4)
lam4 = ci_yj(df4)
ult4=ultimate(drift=TRUE, df4, yeo.johnson, yeo.johnson, lambda=-0.08, "..", lambda=-0.08, inverse=TRUE)
#MSE is 5.48115818568678"
tult4 = tult=ultimate(drift=TRUE, df4, transform, inversebc, "..", lam=0.3, dat2=df4)
#MSE is 5.37862476237655
#check AIC:
m4 =list(c(0, 1, 1))
information_calculator(m4)

#==============================================================
#==============================================================
#==============================================================

## QUestion 5
df5 = read.csv("q5_train.csv")
df5 = df5[,2]
df5 = ts(df5,freq = 365.25/7, start = 2004)

t5 = transform(df5)
lam5 = ci_yj(df5)
ult5=ultimate(drift=TRUE, df5, yeo.johnson, yeo.johnson, lambda=-0.5, "..", lambda=-0.5, inverse=TRUE)
#MSE is 6.48513650163962
tult5 = tult=ultimate(drift=TRUE, df5, transform, inversebc, "..", lam=0.3, dat2=df5)
#MSE is 6.68539016540819
m5 = list(c(1, 1, 0), c(0, 1, 1))
information_calculator(m5)
