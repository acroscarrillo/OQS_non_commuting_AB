include("../../src/src.jl")

using ProgressBars
using Statistics

L_array = Vector(10:10:80)
p_array = Vector(0.5:0.025:1)
runs_array = Vector(2000:-250:250)
data_array = zeros( length(L_array)*length(p_array), 5 ) # MI, MI_err, L, p, L_A

counter = 1
for (j,p) in ProgressBar(enumerate(p_array))
    for (i, L) in enumerate(L_array)
        L_A = L÷3
        subsys_A = Vector(1:L_A)
        subsys_B = Vector(2*L_A:L)
        temp = zeros(runs_array[i])
        for j=1:runs_array[i]
            A = pwr_law_mat(L,p)
            B = pwr_law_mat(L,p)
            C = correlation_steady_state(A, B)
            temp[j] = MI_NESS(C,subsys_A,subsys_B)
        end
        data_array[counter,:] .= mean(temp), std(temp)/sqrt(runs_array[i]), L, p, L_A
        counter += 1
    end
end

df_wb = DataFrame(data_array, ["MI","MI_err","L","p","L_A"])
CSV.write("data/MI_at_L_A_half_wb_zoom.csv", df_wb)