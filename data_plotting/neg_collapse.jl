include("../src/src.jl")
using Optimization, OptimizationOptimJL

df_in = DataFrame(CSV.File("data/neg_at_L_A_half_zoom_wb.csv")) #format neg|neg_err|L|p|L_A

lower = [0.8, 0.8]
upper = [1.2, 6]
initial_x = [1.0,2.0]
temp_F(params) = log( fss_cost_neg(params, df_in) )
inner_optimizer = NelderMead()
results = optimize(temp_F, lower, upper, initial_x,  Fminbox(inner_optimizer),Optim.Options(g_tol=1e-8))
Optim.minimizer(results)

runs = 100
p_c_array, nu_array = zeros(runs), zeros(runs)
for n in tqdm(Vector(1:runs))
    df_perturbed = perturb_df_neg(df_in)
    temp_F(params) = log( fss_cost_neg(params, df_perturbed) )
    inner_optimizer = NelderMead()
    results = optimize(temp_F, lower, upper, initial_x, Fminbox(inner_optimizer), Optim.Options(g_tol=1e-12))
    p_c_array[n], nu_array[n] =  Optim.minimizer(results)
end

p_c, nu, p_c_err, nu_err = mean(p_c_array), mean(nu_array),  std(p_c_array)/sqrt(runs), std(nu_array)/sqrt(runs)

display( "Cost="*string( temp_F( (p_c, nu) ) )*", p_c="*string(p_c)*"±"*string(p_c_err)*", nu="*string(nu)*"±"*string(nu_err) )





temp_F(params) = log( fss_cost_neg(params, df_in) )

lower = [0.94, 1]
upper = [0.95, 7]
initial_x = [0.95,2.0]
inner_optimizer = NelderMead()
results = optimize(temp_F, lower, upper, initial_x, Fminbox(inner_optimizer), Optim.Options(g_tol=1e-12))

p_c, nu = Optim.minimizer(results)
display( "Cost="*string( temp_F( (p_c, nu) ) )*", p_c="*string(p_c)*", nu="*string(nu) )
