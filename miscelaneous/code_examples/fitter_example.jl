using LsqFit
using Optim

model(x, p) = p[1]*exp.(-x*p[2])
xdata = range(0, stop=5, length=300)
ydata = model(xdata, [1.0 2.0]) + 0.05*randn(length(xdata))
p0 = [0.5, 0.5]
fit = curve_fit(model, xdata, ydata, p0,lower=[0.0,0.0])

plot(xdata, [ydata,model(xdata,fit.param)],title="Fit with p="*string(round.(fit.param,sigdigits=4)))

####
# another one 
####

using Optim
f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

lower = [1.25, -2.1]
upper = [Inf, Inf]
initial_x = [2.0, 2.0]
inner_optimizer = NelderMead()
results = optimize(f, lower, upper, initial_x, Fminbox(inner_optimizer))

####
# another one 
####

# Import the package and define the problem to optimize
using Optimization
rosenbrock(u, p) = (p[1] - u[1])^2 + p[2] * (u[2] - u[1]^2)^2
u0 = zeros(2)
p = [1.0, 100.0]

prob = OptimizationProblem(rosenbrock, u0, p)

# Import a solver package and solve the optimization problem
using OptimizationOptimJL
sol = solve(prob, NelderMead())

# Import a different solver package and solve the optimization problem a different way
using OptimizationBBO
prob = OptimizationProblem(rosenbrock, u0, p, lb = [-1.0, -1.0], ub = [1.0, 1.0])
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited())