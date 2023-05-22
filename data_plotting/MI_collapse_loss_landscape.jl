include("../src/src.jl")

df_in = DataFrame(CSV.File("data/MI_at_L_A_half_wb_zoom.csv")) #format MI|MI_err|L|p|L_A

temp_F(params) = fss_cost(params, df_in; g_noise=false)

p_c_array = Vector(range(0, stop=2, length=1000))
nu_array  = Vector(range(0.5, stop=10, length=1000))

loss = zeros((length(p_c_array),length(nu_array)))
for (j,p_c) in ProgressBar(enumerate(p_c_array))
    for (k,nu) in enumerate(nu_array)
        loss[j,k] = log(temp_F((p_c,nu)))
    end
end

heatmap(nu_array,p_c_array,loss,ylabel=L"p_c",xlabel=L"\nu",colorbar_title=L"\log(Cost)")

savefig("figs/MI_collapse_loss_landscape/loss_landscape.pdf")
savefig("figs/MI_collapse_loss_landscape/loss_landscape.png")