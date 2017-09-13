df3 = read.csv("q3_train.csv")
df3 = df3[,2]

#plot dataset:
corr_plot(df3) # guess with MA(3) or MA(2) component


#log transform by the min:
df3a = df3
df3_min = min(df3a)
df3a_pre = df3a - df3_min + 0.001
df3a= log(df3a_pre)
df3a_diff = diff(df3a)
tsdisplay(df3a_diff)
# doesn't seem to work
# try difference on the original data:
df3b = df3
df3b= diff(df3b)
tsdisplay(df3b)
#difference on log transform: 
df3c = diff(df3a)
corr_plot(df3c)


df3e= read.csv("q3_train.csv")
preProcess(df3e, method = c("YeoJohnson")) #YJ coef = -0.22, 
#choose 0.35 within CI instead for convenience
trans3 = yjPower(df3, -0.25, jacobian.adjusted = FALSE)
trans3_diff = diff(trans3)
corr_plot(trans3_diff) #it's quite stationary. 
#Guess AR(2) seasonal and AR(2) non seasonal 
# may have MA(1) seasonal and MA(1 or 2) non seasonal 
Box.test(trans3_diff, type = "Ljung-Box")# p = 0.003909, stationary
adf.test(trans3_diff, alternative = "stationary") #p-value = 0.01(small) indicates stationarity
kpss.test(trans3_diff) #null hypothesis is that data is stationary. p-value=0.1 --> data is stationary for alpha=0.05






#tbats model for q3: 

