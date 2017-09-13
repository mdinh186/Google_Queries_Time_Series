
q1tbats = tbats(trants)
fpc2 = forecast(q1tbats, h=104)
plot(fpc2, ylab="thousands of barrels per day")

#decomposition of seasonality
plot.tbats(q1tbats)

#take care of multiple seasonality
x.msts = msts(trants,seasonal.periods=c(7,30.5, 365.25))
fit = tbats(x.msts)

plot(residuals(fit))
fc = forecast(fit, h = 104)
plot(fc)

accuracy(fc)
plot(density(residuals(fit)))



