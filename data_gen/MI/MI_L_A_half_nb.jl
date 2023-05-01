include("../src/src.jl")

L_array = Vector(10:10:50)
p_array = Vector(0:0.1:2)
runs_array = [1000,1000,1000,200,200]
data_array = zeros( length(L_array)*length(p_array), 5 ) # MI, MI_err, L, p, L_A

counter = 1
for (j,p) in enumerate(p_array)
    display( (j-1) / length(p_array))
    for (i, L) in enumerate(L_array)
        L_A = L÷2
        subsys_A = Vector(1:L_A)
        temp = zeros(runs_array[i])
        for j=1:runs_array[i]
            temp[j] = MI_NESS(zeros(L,L),A(L,p),A(L,p),subsys_A)
        end
        data_array[counter,:] .= mean(temp), std(temp)/sqrt(runs_array[i]), L, p, L_A
        counter += 1
    end
end

df_nb = DataFrame(data_array, ["MI","MI_err","L","p","L_A"])
CSV.write("data/MI_at_L_A_half_nb.csv", df_nb)