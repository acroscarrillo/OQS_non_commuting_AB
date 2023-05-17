include("../../src/src.jl")

using ProgressBars

L=10
n_samples = 1000
p_array = [0.4,4]
h = zeros(L,L)

df_spectrum = DataFrame(zeros(0,3), ["lamb","L","p"] )
for _=tqdm(1:n_samples)
    for p in p_array
        A = pwr_law_mat(L, p,false)
        B = pwr_law_mat(L, p,false)
        M = master_op(h,A, B, false)
        data = hcat( eigvals(M), L*ones(L^2), p*ones(L^2) )

        df_2_append = DataFrame(data, ["lamb","L","p"] )
        df_spectrum = vcat(df_spectrum, df_2_append)
    end
end

CSV.write("data/M_spectrum.csv", df_spectrum)