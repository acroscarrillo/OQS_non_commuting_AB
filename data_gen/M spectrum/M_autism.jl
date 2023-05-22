include("../../src/src.jl")

using ProgressBars

L=1000
n_samples = 100
p_array = [0.4,4]
h = zeros(L,L)

df_spectrum = DataFrame(zeros(0,3), ["lamb","L","p"] )
for _=tqdm(1:n_samples)
    for p in p_array
        A, B = pwr_law_mat(L, p,false), pwr_law_mat(L, p,false)
        M_lambs = M_spectrum(A,B)
        data = hcat( M_lambs, L*ones(L), p*ones(L) )

        df_2_append = DataFrame(data, ["lamb","L","p"] )
        df_spectrum = vcat(df_spectrum, df_2_append)
    end
end

CSV.write("data/M_spectrum.csv", df_spectrum)