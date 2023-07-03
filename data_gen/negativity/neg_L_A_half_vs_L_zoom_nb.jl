include("../../src/src.jl")

using ProgressBars
using Statistics
using CSV 

L_array = Vector(400:400:2000)
p_array = Vector(0.90:0.01:1.1)
runs_array = Vector(1000:-200:200)
data_array = zeros( length(L_array)*length(p_array), 5 ) # MI, MI_err, L, p, L_A

counter = 1
for (_,p) in ProgressBar(enumerate(p_array))
    for (i, L) in enumerate(L_array)
        subsys_A = Vector(1:L÷2)
        subsys_B = Vector(L÷2+1:L)
        temp = zeros(runs_array[i])
        for j=1:runs_array[i]
            A = pwr_law_mat(L,p)
            # B = pwr_law_mat(L,p)
            B = I(L)-A
            C = correlation_ness(A, B)
            temp[j] = negativity(C,subsys_A,subsys_B)
        end
        data_array[counter,:] .= mean(temp), std(temp)/sqrt(runs_array[i]), L, p, L÷2
        counter += 1
    end
end

df_wb = DataFrame(data_array, ["neg","neg_err","L","p","L_A"])
CSV.write("data/neg_at_L_A_half_zoom_nb.csv", df_wb)