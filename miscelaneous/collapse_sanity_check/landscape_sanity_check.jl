include("../../src/src.jl")
using Optimization
using NPZ

old_data = npzread("miscelaneous/collapse_sanity_check/old_collapse_data.npy")

df_in = DataFrame(old_data, ["MI", "MI_err", "L", "p", "L_A"]) #format MI|MI_err|L|p|L_A

temp_F(params) = fss_cost(params, df_in)

p_c_array = Vector(range(0.75, stop=1.25, length=100))
nu_array  = Vector(range(4, stop=10, length=100))

loss = zeros((length(p_c_array),length(nu_array)))
for (j,p_c) in ProgressBar(enumerate(p_c_array))
    for (k,nu) in enumerate(nu_array)
        loss[j,k] = log(temp_F((p_c,nu)))
    end
end

heatmap(nu_array,p_c_array,loss,ylabel=L"p_c",xlabel=L"\nu",colorbar_title=L"\log(Cost)")