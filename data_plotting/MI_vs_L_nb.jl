@time begin
using CSV
using Plots 
using DataFrames

df = DataFrame(CSV.File("data/MI_at_L_A_half_nb.csv"))

p_array = unique(df.p)
L_array = unique(df.L)
MI, MI_err, L = zeros(0, length(L_array)), zeros(0, length(L_array)), zeros(0, length(L_array))
for val in p_array
    df_subset = df[df.p .== val,:]
    temp_MI, temp_MI_err, temp_L = transpose(df_subset.MI), transpose(df_subset.MI_err), transpose(df_subset.L)
    global MI, MI_err, L = vcat(MI,temp_MI), vcat(MI_err,temp_MI_err), vcat(L,temp_L)
end

plot(transpose(float.(L)), transpose(float.(MI)), yerr = transpose(float.(MI_err)),
        labels="p = ".*string.(transpose(p_array)), xlabel="L", ylabel="MI",dpi=300,
        xaxis=:log, yaxis=:log)

# savefig("figs/MI_vs_L/MI_vs_L_nb.pdf")
savefig("figs/MI_vs_L/35.png")
end

@time begin
    using CSV
    using Plots 
    using DataFrames
    
    df = DataFrame(CSV.File("data/MI_at_L_A_half_nb.csv"))
    
    p_array = unique(df.p)
    L_array = unique(df.L)
    MI, MI_err, L = zeros(0, length(L_array)), zeros(0, length(L_array)), zeros(0, length(L_array))
    for val in p_array
        df_subset = df[df.p .== val,:]
        temp_MI, temp_MI_err, temp_L = transpose(df_subset.MI), transpose(df_subset.MI_err), transpose(df_subset.L)
        global MI, MI_err, L = vcat(MI,temp_MI), vcat(MI_err,temp_MI_err), vcat(L,temp_L)
    end
    
    plot(transpose(float.(L)), transpose(float.(MI)), yerr = transpose(float.(MI_err)),
            labels="p = ".*string.(transpose(p_array)), xlabel="L", ylabel="MI",dpi=300,
            xaxis=:log, yaxis=:log)
    
    # savefig("figs/MI_vs_L/MI_vs_L_nb.pdf")
    savefig("figs/MI_vs_L/35.png")
end