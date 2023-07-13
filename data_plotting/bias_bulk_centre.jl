include("../src/src.jl")

df = DataFrame(CSV.File("data/bias_bulk.csv"))

temps_4_hist = vcat(filter(row -> row.p == 4, df).std',filter(row -> row.p == .4, df).std')

histogram(temps_4_hist', fillalpha=0.5, normalize=true, xlabel="σ", ylabel="f(σ)",dpi=300,label = ["p=4" "p=0.4"])

savefig("figs/temp/std_bias_bulk.pdf")
savefig("figs/temp/std_bias_bulk.png")


temps_4_hist = vcat(filter(row -> row.p == 4, df).mean',filter(row -> row.p == .4, df).mean')

histogram(temps_4_hist', fillalpha=0.5, normalize=true, xlabel="mean", ylabel="f(mean)",dpi=300,label = ["p=4" "p=0.4"])

savefig("figs/temp/mean_bias_bulk.pdf")
savefig("figs/temp/mean_bias_bulk.png")