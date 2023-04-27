include("../src/src.jl")

L=30
n_samples = 1000
p_array = [0.4,4]

temp_array = zeros(n_samples, length(p_array))
for i=1:n_samples
    display((i-1)/n_samples)
    for j=1:length(p_array)
        temp_array[i,j] = temp_bulk(   C_NESS(zeros(L,L), A(L,p_array[j]), A(L,p_array[j]))    )
    end
end

df_temp = DataFrame(temp_array, ["temp","p"]) #This has a wrong format
CSV.write("data/temps.csv", df_wb)

histogram(temp_array, fillalpha=0.5,normalize=true, xlabel="T", ylabel="P(T)",dpi=300)

savefig("figs/temp/temp.pdf")
savefig("figs/temp/temp.png")