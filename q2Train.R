library(tsoutliers)


#helper function: 

corr_plot = function (data){
  par(mfrow = c(3,1))
  lag = 104
  plot.ts(data)
  acf(data, lag.max = lag)
  pacf(data, lag.max = lag) 
}

#=====================================================
## Question2

df2 = read.csv("q2_train.csv")
df2 = df2[,2]

#plot dataset:
plot.ts(df2, type = "l") #doesn't look stationary
corr_plot(df2) # guess with MA(3) or MA(2) component

#log transform by the min:
df2a = df2
df2_min = min(df2a)
df2a_pre= df2a - df2_min + 0.001
df2a= log(df2a_pre)
corr_plot(df2a)

# try difference on the original data:
df2b = df2
df2b= diff(df2b)
corr_plot(df2b)

#try difference on the log transform:
df2c =diff(df2a)
corr_plot(df2c)

#try to take log transform of differenced data
df2d = df2b
df2d_min = min(df2d)
df2d= log(df2d - df2d_min + 0.001)
corr_plot(df2d) # seems like there's a level shift


df2_pre = read.csv("q2_train.csv")
preProcess(df2_pre, method = c("YeoJohnson")) #YJ coef = 0.33, 
#choose 0.35 within CI instead for convenience
trans2 = yjPower(df2, 0.35, jacobian.adjusted = FALSE)
trans2_diff = diff(trans2)
corr_plot(trans2_diff) #it's quite stationary. 
#Guess AR(2) seasonal and AR(2) non seasonal 
# may have MA(1) seasonal and MA(1 or 2) non seasonal 
Box.test(trans2_diff, type = "Ljung-Box")# p = 9.72e-10, stationary
