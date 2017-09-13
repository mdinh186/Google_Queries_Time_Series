df5 = read.csv("q5_train.csv")
df5 = df5[,2]

#plot dataset:
plot.ts(df5, type = "l") #doesn't look stationary
corr_plot(df5) # guess with MA(3) or MA(2) component





# try difference on the original data:
df5b = df5
df5b= diff(df5b)
df5b_min = min(df5b)
df5b= log(df5b - df5b_min + 0.001)
corr_plot(df5b)

#try difference on the log transform:
df5c = diff(df5a)
plot.ts(df5c)
corr_plot(df5c)

#try to take log transform of differenced data
df5d = df5b
df5d_min = min(df5d)
df5d= log(df5d - df5d_min + 0.001)
corr_plot(df5d) # seems like there's a level shift



df5_pre = read.csv("q5_train.csv")
preProcess(df5_pre, method = c("YeoJohnson")) #YJ coef = -0.08, 
#choose 0.35 within CI instead for convenience
trans5 = yjPower(df5, -0.08, jacobian.adjusted = FALSE)
trans5_diff = diff(trans5)
corr_plot(trans5_diff) #it's quite stationary. 
# Guess AR(2) seasonal and AR(2) non seasonal 
# may have MA(1) seasonal and MA(1 or 2) non seasonal 
Box.test(trans5_diff, type = "Ljung-Box")# p = 0.003909, stationary
adf.test(trans5_diff, alternative = "stationary") #p-value = 0.01(small) indicates stationarity
kpss.test(trans5_diff) #null hypothesis is that data is stationary. p-value=0.1 --> data is stationary for alpha=0.05


#tbats model for q4: 

