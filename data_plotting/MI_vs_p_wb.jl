include("../src/src.jl")

df = DataFrame(CSV.File("data/MI_at_L_A_half_wb.csv"))

p_array = unique(df.p)
L_array = unique(df.L)
MI, MI_err, p = reshape([], 0, length(p_array)),reshape([], 0, length(p_array)),reshape([], 0, length(p_array))
for val in L_array
    df_subset = df[df.L .== val,:]
    temp_MI, temp_MI_err, temp_p = transpose(df_subset.MI), transpose(df_subset.MI_err), transpose(df_subset.p)
    MI, MI_err, p = vcat(MI,temp_MI), vcat(MI_err,temp_MI_err), vcat(p,temp_p)
end

plot(transpose(float.(p)), transpose(float.(MI)), yerr = transpose(float.(MI_err)),
        labels="L = ".*string.(transpose/Users/alex/Documents/GitHub/OQS_non_commuting_AB/figs/MI_vs_p.png(L_array)), xlabel="p", ylabel="MI",dpi=300)

savefig("figs/MI_vs_p.pdf")
savefig("figs/MI_vs_p.png")