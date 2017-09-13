df4 = read.csv("q4_train.csv")
df4 = df4[,2]

#plot dataset:
plot.ts(df4, type = "l") #doesn't look stationary
corr_plot(df4) # guess with MA(3) or MA(2) component





# try difference on the original data:
df4b = df4
df4b= diff(df4b)
df4b_min = min(df4b)
df4b_pre = df4b - df4b_min + 0.001
df4b= log(df4b_pre)
corr_plot(df4b)

#try difference on the log transform:
df4c = diff(df4a)
plot.ts(df4c)
corr_plot(df4c)

#try to take log transform of differenced data
df4d = df4b
df4d_min = min(df4d)
df4d= log(df4d - df4d_min + 0.001)
corr_plot(df4d) # seems like there's a level shift



df4_pre = read.csv("q4_train.csv")
preProcess(df4_pre, method = c("YeoJohnson")) #YJ coef = -0.08, 
#choose 0.35 within CI instead for convenience
trans4 = yjPower(df4, -0.08, jacobian.adjusted = FALSE)
trans4_diff = diff(trans4)
corr_plot(trans4_diff) #it's quite stationary. 
#Guess AR(2) seasonal and AR(2) non seasonal 
# may have MA(1) seasonal and MA(1 or 2) non seasonal 
Box.test(trans4_diff, type = "Ljung-Box")# p = 0.003909, stationary
adf.test(trans4_diff, alternative = "stationary") #p-value = 0.01(small) indicates stationarity
kpss.test(trans4_diff) #null hypothesis is that data is stationary. p-value=0.1 --> data is stationary for alpha=0.05

transform = function(dat){
  lambda = BoxCox.lambda(dat, method = "loglik")
  dfbc = BoxCox(dat, lambda = lambda)
  trans_diff = diff(dfbc)
  corr_plot(trans_diff)
  return (dat)
}

#tbats model for q4: 

