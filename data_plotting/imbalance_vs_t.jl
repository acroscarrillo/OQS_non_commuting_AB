include("../src/src.jl")


df = DataFrame(CSV.File("data/imb_vs_t.csv"))

t_array = unique(df.t)
p_array = unique(df.p)
imb, imb_err, t = reshape([], 0, length(t_array)),reshape([], 0, length(t_array)),reshape([], 0, length(t_array))
for val in p_array
    df_subset = df[df.p .== val,:]
    temp_imb, temp_imb_err, temp_t = transpose(df_subset.imb), transpose(df_subset.imb_err), transpose(df_subset.t)
    imb, imb_err, t = vcat(imb,temp_imb), vcat(imb_err,temp_imb_err), vcat(t,temp_t)
end

plot(transpose(float.(t)), transpose(float.(imb)), yerr = transpose(float.(imb_err)),
        labels="p = ".*string.(transpose(p_array)), xlabel="t", ylabel="imb",dpi=300,xaxis=:lin, yaxis=:log)

savefig("figs/imbalance/imb_vs_t.pdf")
savefig("figs/imbalance/imb_vs_t.png")