include("../../src/src.jl")
using ProgressBars

L=1000
n_samples = 1000
p_array = [0.4,4]

temp_array = zeros(n_samples*length(p_array), 4) #std, mean, L, p
counter =1
for (j,p) in enumerate(p_array)
    for i in tqdm(1:n_samples)
        # A, B = pwr_law_mat(L, p), pwr_law_mat(L, p)
        A = pwr_law_mat(L, p)
        B = I(L)-A
        C = correlation_ness(A, B)
        temp_array[counter,:] .= bulk_bias_variance(C,L÷5)...,L,p
        counter += 1
    end
end

df_temp = DataFrame(temp_array, ["mean","std","L","p"]) 
CSV.write("data/bias_bulk.csv", df_temp)

temps_4_hist = vcat(filter(row -> row.p == 4, df_temp).std',filter(row -> row.p == .4, df_temp).std')

histogram(temps_4_hist', fillalpha=0.5, normalize=true, xlabel=L"\tilde{n}", ylabel=L"f(\tilde{n})",dpi=300,label = [L"p=4" L"p=0.4"])

savefig("figs/temp/bias_bulk.pdf")
savefig("figs/temp/bias_bulk.png")