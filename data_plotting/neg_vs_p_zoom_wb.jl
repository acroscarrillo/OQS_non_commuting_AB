include("../src/src.jl")

df = DataFrame(CSV.File("data/neg_at_L_A_half_zoom_wb.csv"))

p_array = unique(df.p)
L_array = unique(df.L)
neg, neg_err, p = reshape([], 0, length(p_array)),reshape([], 0, length(p_array)),reshape([], 0, length(p_array))
for val in L_array
    df_subset = df[df.L .== val,:]
    temp_neg, temp_neg_err, temp_p = transpose(df_subset.neg), transpose(df_subset.neg_err), transpose(df_subset.p)
    neg, neg_err, p = vcat(neg,temp_neg), vcat(neg_err,temp_neg_err), vcat(p,temp_p)
end

plot(transpose(float.(p)), transpose(float.(neg)), yerr = transpose(float.(neg_err)),
        labels="L = ".*string.(transpose(L_array)), xlabel="p", ylabel="neg",dpi=300)

        
savefig("figs/neg_vs_p/neg_vs_p_zoom_wb.pdf")
savefig("figs/neg_vs_p/neg_vs_p_zoom_wb.png")