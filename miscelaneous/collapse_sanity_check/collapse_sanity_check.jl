include("../../src/src.jl")
using Optim
using NPZ

old_data = npzread("miscelaneous/collapse_sanity_check/old_collapse_data.npy")

df_in = DataFrame(old_data, ["MI", "MI_err", "L", "p", "L_A"]) #format MI|MI_err|L|p|L_A

lower = [0.1, 0.1]
upper = [10.0, 10.0]
initial_x = [1.0,1.0]
temp_F(params) = fss_cost(params, df_in)
results = optimize(temp_F, lower, upper, initial_x, NelderMead(),Optim.Options(g_tol=1e-8))
Optim.minimizer(results)

runs = 100
p_c_array, nu_array = zeros(runs), zeros(runs)
for n in tqdm(Vector(1:runs))
    df_perturbed = perturb_df(df_in)
    temp_F(params) = fss_cost(params, df_perturbed)
    results = optimize(temp_F, lower, upper, initial_x, Fminbox(inner_optimizer))
    p_c_array[n], nu_array[n] =  Optim.minimizer(results)
end

p_c, nu, p_c_err, nu_err = mean(p_c_array), mean(nu_array),  std(p_c_array)/sqrt(runs), std(nu_array)/sqrt(runs)

display( "Cost="*string( temp_F( (p_c, nu) ) )*", p_c="*string(p_c)*"±"*string(p_c_err)*", nu="*string(nu)*"±"*string(nu_err) )