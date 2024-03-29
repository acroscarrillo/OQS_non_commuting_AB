include("../../src/src.jl")

using ProgressBars
using Statistics
using CSV 

L = 100
L_A_array = Vector(13:3:L÷2)
p_array = [0.0, 0.2, 0.6, 1.0, 1.6]
runs = 1000
data_array = zeros( length(L_A_array)*length(p_array), 5 ) # MI, MI_err, L, p, L_A

counter = 1
for (_,p) in ProgressBar(enumerate(p_array))
    for (_, L_A) in ProgressBar(enumerate(L_A_array))
        subsys_A,subsys_B = Vector(1:L_A), Vector(L_A+1:L)
        temp = zeros(runs)
        for j=1:runs 
            A = pwr_law_mat(L,p)
            # B = pwr_law_mat(L,p)
            B = I(L)-A
            C = correlation_ness(A, B)
            C_sub = sub_correlation(C, vcat(subsys_A,subsys_B)) # trace boundary out
            L_AB = length(vcat(subsys_A,subsys_B)) # new effective L after trace
            subsys_A, subsys_B = Vector(1:L_A), Vector(L_A+1:L_AB) # due to the trace
            temp[j] = negativity(C_sub,subsys_A,subsys_B)
        end
        data_array[counter,:] .= mean(temp), std(temp)/sqrt(runs), L, p, L_A
        counter += 1
    end
end

df_wb = DataFrame(data_array, ["neg","neg_err","L","p","L_A"])
CSV.write("data/neg_vs_L_A_nb.csv", df_wb)