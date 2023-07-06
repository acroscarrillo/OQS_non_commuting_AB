include("../../src/src.jl")

using ProgressBars
using Statistics
using CSV 

L = 500
t_array = Vector(0:0.1:10)
p_array = [0.0,4.0]
runs = 100
data_array = zeros( length(t_array)*length(p_array), 5 ) # imb, imb_err, L, p, t

C_0 = Neel_corr(L)

counter = 1
for (_,p) in ProgressBar(enumerate(p_array))
    temp = zeros(runs, length(t_array))
    for (_,j) in ProgressBar(enumerate(1:runs))
        for (k, t) in enumerate(t_array)
            A = pwr_law_mat(L,p)
            # B = pwr_law_mat(L,p)
            B = I(L)-A
            C = correlation(A, B, C_0, t)
            temp[j,k] = imbalance(C)
        end
    end
    mean_t, err_t = mean(temp,dims=1), std(temp,dims=1)/sqrt(runs)
    for (n, t) in ProgressBar(enumerate(t_array))
        data_array[counter,:] .=  mean_t[n], err_t[n], L, p, t
        counter += 1
    end
end

df_wb = DataFrame(data_array, ["imb","imb_err","L","p","t"])
CSV.write("data/imb_vs_t.csv", df_wb)