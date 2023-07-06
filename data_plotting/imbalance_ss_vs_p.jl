include("../src/src.jl")


df = DataFrame(CSV.File("data/imbalance_ss.csv"))

p_array = unique(df.p)
L_array = unique(df.L)
imb, imb_err, p = reshape([], 0, length(p_array)),reshape([], 0, length(p_array)),reshape([], 0, length(p_array))
for val in L_array
    df_subset = df[df.L .== val,:]
    temp_imb, temp_imb_err, temp_p = transpose(df_subset.imb), transpose(df_subset.imb_err), transpose(df_subset.p)
    imb, imb_err, p = vcat(imb,temp_imb), vcat(imb_err,temp_imb_err), vcat(p,temp_p)
end

plot(transpose(float.(p)), transpose(float.(imb)), yerr = transpose(float.(imb_err)),
        labels="L = ".*string.(transpose(L_array)), xlabel="p", ylabel="imb",dpi=300)

savefig("figs/imbalance/imb_vs_p.pdf")
savefig("figs/imbalance/imb_vs_p.png")