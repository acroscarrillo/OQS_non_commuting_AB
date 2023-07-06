include("../../src/src.jl")

using ProgressBars
using Statistics
using CSV 


L_array = Vector(400:400:2000)
p_array = Vector(0:0.2:4)
runs_array = Vector(500:-100:100)
data_array = zeros( length(L_array)*length(p_array), 4 )  # imb, imb_err, L, p

counter = 1
for (_,p) in ProgressBar(enumerate(p_array))
    for (i, L) in ProgressBar(enumerate(L_array))
        temp = zeros(runs_array[i])
        for j=1:runs_array[i]
            A = pwr_law_mat(L,p)
            # B = pwr_law_mat(L,p)
            B = I(L)-A
            C = correlation_ness(A, B)
            temp[j] = imbalance(C)
        end
        data_array[counter,:] .= mean(temp), std(temp)/sqrt(runs_array[i]), L, p
        counter += 1
    end
end

df_wb = DataFrame(data_array, ["imb","imb_err","L","p"])
CSV.write("data/imbalance_ss.csv", df_wb)