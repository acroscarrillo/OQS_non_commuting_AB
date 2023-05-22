include("../../src/src.jl")

using ProgressBars

L=1000
n_samples = 100
p_array = Vector(0:0.1:2) 

df_gap = DataFrame(zeros(0,8),["tau","lamb_min","lamb_delta","tau_err","lamb_min_err","lamb_delta_err","L","p"] )
for p in tqdm(p_array)
    tau_array, lamb_min_array, lamb_delta_array = zeros(n_samples), zeros(n_samples), zeros(n_samples)
    for j=1:n_samples
        A, B = pwr_law_mat(L, p), pwr_law_mat(L, p)
        M_lambs = M_spectrum(A,B)

        tau_array[j] = 1/M_lambs[1]
        lamb_min_array[j] = M_lambs[1]
        lamb_delta_array[j] = M_lambs[2]-M_lambs[1]
    end
    data = hcat( mean(tau_array), mean(lamb_min_array), mean(lamb_delta_array), std(tau_array)/sqrt(n_samples), std(lamb_min_array)/sqrt(n_samples), std(lamb_delta_array)/sqrt(n_samples), L, p )
    df_2_append = DataFrame(data, ["tau","lamb_min","lamb_delta","tau_err","lamb_min_err","lamb_delta_err","L","p"] )
    df_gap = vcat(df_gap, df_2_append)
end

CSV.write("data/spectral_gap.csv", df_gap)