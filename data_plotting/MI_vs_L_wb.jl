include("../src/src.jl")

df = DataFrame(CSV.File("data/MI_at_L_A_half_wb.csv"))

p_array = unique(df.p)
L_array = unique(df.L)
MI, MI_err, L = reshape([], 0, length(L_array)),reshape([], 0, length(L_array)),reshape([], 0, length(L_array))
for val in p_array
    df_subset = df[df.p .== val,:]
    temp_MI, temp_MI_err, temp_L = transpose(df_subset.MI), transpose(df_subset.MI_err), transpose(df_subset.L)
    MI, MI_err, L = vcat(MI,temp_MI), vcat(MI_err,temp_MI_err), vcat(L,temp_L)
end

plot(transpose(float.(L)), transpose(float.(MI)), yerr = transpose(float.(MI_err)),
        labels="p = ".*string.(transpose(p_array)), xlabel="L", ylabel="MI",dpi=300,
        xaxis=:log, yaxis=:log)

savefig("figs/MI_vs_L.pdf")
savefig("figs/MI_vs_L.png")