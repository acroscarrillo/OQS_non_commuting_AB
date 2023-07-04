include("../src/src.jl")
using LaTeXStrings
df_in = DataFrame(CSV.File("data/MI_at_L_A_half_wb_zoom.csv")) #format MI|MI_err|L|p|L_A



p_c_array = Vector(range(0.8, stop=1.75, length=200))
nu_array  = Vector(range(4, stop=8, length=200))


runs = 100
p_c_min = zeros(runs)
nu_min = zeros(runs)
for n in tqdm(Vector(1:runs))
    df_perturbed = perturb_df(df_in)
    temp_F(params) = fss_cost(params, df_perturbed)
    loss = zeros((length(p_c_array),length(nu_array)))
    for (j,p_c) in enumerate(p_c_array)
        for (k,nu) in enumerate(nu_array)
            loss[j,k] = log(temp_F((p_c,nu)))
        end
    end
    p_c_min_index, nu_min_index = argmin(loss)[1], argmin(loss)[2]
    p_c_min[n], nu_min[n] = p_c_array[p_c_min_index],  nu_array[nu_min_index]
end

p_c, nu = mean(p_c_min), mean(nu_min)
p_c_err, nu_err = std(p_c_min)/sqrt(runs), std(nu_min)/sqrt(runs)

heatmap(nu_array,p_c_array,loss,ylabel=L"p_c",xlabel=L"\nu",colorbar_title=L"\log(Cost)")
heatmap(nu_array,p_c_array,loss,ylabel=L"p_c",xlabel=L"\nu",colorbar_title=L"\log(Cost)")

savefig("figs/MI_collapse_loss_landscape/loss_landscape.pdf")
savefig("figs/MI_collapse_loss_landscape/loss_landscape.png")