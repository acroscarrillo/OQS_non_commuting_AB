include("../src/src.jl")

using Plots

df = DataFrame(CSV.File("data/neg_at_L_A_half.csv"))

p_array = unique(df.p)
L_array = unique(df.L)
neg, neg_err, L = reshape([], 0, length(L_array)), reshape([], 0, length(L_array)), 
                reshape([], 0, length(L_array))
for val in p_array
    df_subset = df[df.p .== val,:]
    temp_neg, temp_neg_err, temp_L = transpose(df_subset.neg), transpose(df_subset.neg_err), transpose(df_subset.L)
    neg, neg_err, L = vcat(neg,temp_neg), vcat(neg_err,temp_neg_err), vcat(L,temp_L)
end

plot(transpose(float.(L)), transpose(float.(neg)), yerr = transpose(float.(neg_err)),
        labels="p = ".*string.(transpose(p_array)), xlabel="L", ylabel="neg",dpi=300,
        xaxis=:lin, yaxis=:lin)

savefig("figs/neg_vs_L/neg_vs_L_nb.pdf")
savefig("figs/neg_vs_L/neg_vs_L_nb.png")