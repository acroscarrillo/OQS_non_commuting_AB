include("../src/src.jl")

using Plots

df = DataFrame(CSV.File("data/neg_vs_L_A_wb.csv"))

p_array = unique(df.p)
L_A_array = unique(df.L_A)
neg, neg_err, L_A = reshape([], 0, length(L_A_array)), reshape([], 0, length(L_A_array)), 
                reshape([], 0, length(L_A_array))
for val in p_array
    df_subset = df[df.p .== val,:]
    temp_neg, temp_neg_err, temp_L_A = transpose(df_subset.neg), transpose(df_subset.neg_err), transpose(df_subset.L_A)
    neg, neg_err, L_A = vcat(neg,temp_neg), vcat(neg_err,temp_neg_err), vcat(L_A,temp_L_A)
end

plot(transpose(float.(L_A)), transpose(float.(neg)), yerr = transpose(float.(neg_err)),
        labels="p = ".*string.(transpose(p_array)), xlabel="L_A", ylabel="neg",dpi=300,
        xaxis=:log, yaxis=:log)

savefig("figs/neg_vs_L_A/neg_vs_L_A_wb.pdf")
savefig("figs/neg_vs_L_A/neg_vs_L_A_wb.png")