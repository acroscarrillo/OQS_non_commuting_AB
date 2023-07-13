include("../../src/src.jl")

using ProgressBars
using Statistics
using CSV 

L_array = Vector(100:100:1000)
p_array = [0.0, 0.2, 0.6, 1.0, 1.6]
runs_array = Vector(2000:-200:200)
data_array = zeros( length(L_array)*length(p_array), 5 ) # MI, MI_err, L, p, L_A

counter = 1
for (_,p) in ProgressBar(enumerate(p_array))
    for (i, L) in ProgressBar(enumerate(L_array))
        subsys_A, subsys_B = Vector(1:L÷3), Vector(2*L÷3+1:L)
        temp = zeros(runs_array[i])
        for j=1:runs_array[i]
            A = pwr_law_mat(L,p)
            # B = pwr_law_mat(L,p)
            B = I(L)-A
            C = correlation_ness(A, B)
            C_sub = sub_correlation(C, vcat(subsys_A,subsys_B))
            L_AB = length(vcat(subsys_A,subsys_B))
            subsys_A, subsys_B = Vector(1:L÷3), Vector(L÷3+1:L_AB) # due to the trace
            temp[j] = negativity(C_sub,subsys_A,subsys_B)
        end
        data_array[counter,:] .= mean(temp), std(temp)/sqrt(runs_array[i]), L, p, L÷2
        counter += 1
    end
end

df_wb = DataFrame(data_array, ["neg","neg_err","L","p","L_A"])
CSV.write("data/neg_at_L_A_half_wb.csv", df_wb)