using LsqFit
using Optim

model(x, p) = p[1]*exp.(-x*p[2])
xdata = range(0, stop=5, length=300)
ydata = model(xdata, [1.0 2.0]) + 0.05*randn(length(xdata))
p0 = [0.5, 0.5]
fit = curve_fit(model, xdata, ydata, p0,lower=[0.0,0.0])

plot(xdata, [ydata,model(xdata,fit.param)],title="Fit with p="*string(round.(fit.param,sigdigits=4)))
