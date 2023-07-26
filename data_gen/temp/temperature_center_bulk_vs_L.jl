include("../../src/src.jl")
using ProgressBars
using LaTeXStrings

L_array = Vector(10:10:1000)
# L_array = Vector(10:1:100)
n_samples = Vector(1090:-10:100)
p_array = [0.4,4]

data = zeros(length(p_array)*sum(n_samples), 5) #std, mean, L, p, run
counter_data = 1
for (i,L) in tqdm(enumerate(L_array))
    # temp_array = zeros(n_samples[i], 4) #std, mean, L, p
    for (j,p) in enumerate(p_array)
        for i in 1:n_samples[i]
            # A, B = pwr_law_mat(L, p), pwr_law_mat(L, p)
            A = pwr_law_mat(L, p)
            B = I(L)-A
            C = correlation_ness(A, B)
            # temp_array[i,:] .= bulk_bias_variance(C,L÷5)...,L,p
            data[counter_data, :] .= bulk_bias_variance(C,L÷5)..., L, p, i
            counter_data += 1
        end
        # data[counter_data, :] .= mean(temp_array[:,1]), std(temp_array[:,1])/sqrt(length(temp_array[:,1])), L, p
        # counter_data += 1
    end
end

df_temp = DataFrame(data, ["std","mean","L","p","run"]) 
CSV.write("data/bias_bulk_vs_L.csv", df_temp)

df_std = DataFrame(std_mean=Float32[], L=Int32[],p=Float32[]) 
for p_temp in p_array
    for L_temp in L_array
        push!(df_std, [mean(  df_temp[(df_temp.p .== p_temp) .& (df_temp.L .== L_temp), :].std  ), L_temp, p_temp])
    end
end
CSV.write("data/bias_bulk_vs_L_STD.csv", df_std)

temps_y = hcat(df_std[df_std.p .== 4.0,:].std_mean, df_std[df_std.p .== Float32(0.4),:].std_mean)
temps_x = hcat(L_array,L_array)

plot(temps_x, temps_y, ylabel=L"\overline{\sigma}", xlabel=L"L",dpi=300,label = [L"p=4" L"p=0.4"],xaxis=:log,yaxis=:log)

savefig("figs/temp/bias_bulk_vs_L.pdf")
savefig("figs/temp/bias_bulk_vs_L.png")