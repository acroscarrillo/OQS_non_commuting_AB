include("../src/src.jl")
using Optimization

df_in = DataFrame(CSV.File("data/MI_at_L_A_half_wb_zoom.csv")) #format MI|MI_err|L|p|L_A

temp_F(params) = fss_cost(params, df_in; g_noise=false)

lower = [0.1, 0.1]
upper = [10.0, 10.0]
initial_x = [1.0,6.0]
inner_optimizer = GradientDescent()
results = optimize(temp_F, lower, upper, initial_x, Fminbox(inner_optimizer))

p_c, nu = Optim.minimizer(results)
display( "Cost="*string( temp_F( (p_c, nu) ) )*", p_c="*string(p_c)*", nu="*string(nu) )