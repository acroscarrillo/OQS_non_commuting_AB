include("../../src/src.jl")
using ProgressBars

L=100
n_samples = 10000
p_array = [0.4,4]

temp_array = zeros(n_samples*length(p_array), 3)
counter =1
for (j,p) in enumerate(p_array)
    for i in tqdm(1:n_samples)
        A, B = pwr_law_mat(L, p), pwr_law_mat(L, p)
        C = correlation_ness(A, B)
        temp_array[counter,:] .= central_occ_bias(C),L,p
        counter += 1
    end
end

df_temp = DataFrame(temp_array, ["temp","L","p"]) 
CSV.write("data/temps.csv", df_wb)

temps_4_hist = vcat(filter(row -> row.p == 4, df_temp).temp',filter(row -> row.p == .4, df_temp).temp')

histogram(temps_4_hist', fillalpha=0.5, normalize=true, xlabel=L"\tilde{n}", ylabel=L"f(\tilde{n})",dpi=300,label = [L"p=4" L"p=0.4"])

savefig("figs/temp/temp.pdf")
savefig("figs/temp/temp.png")