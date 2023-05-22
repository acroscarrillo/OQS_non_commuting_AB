using LaTeXStrings

df_in = DataFrame(CSV.File("data/M_spectrum.csv")) 


data = vcat(filter(row -> row.p == 4, df_spectrum).lamb',filter(row -> row.p == .4, df_spectrum).lamb')

data = filter(row -> row.p == 4, df_spectrum).lamb'

histogram(data',ylabel=L"f(\lambda_M)",xlabel=L"\lambda_M",label = [L"p=4" L"p=0.4"],fillalpha=.4,normalize=true)

histogram2d(real.(data'),imag.(data'),ylab= "Im(M)",xlab= "Re(M)",fillalpha=0.4,normalize=true,bins=(range(-2, 2,100),range(-2, 2,100)),title= "p=4" )

@df filter(row -> row.p == 4, df_in) histogram(:lamb)

@df filter(row -> row.p == .4, df_in) histogram(:lamb)

savefig("figs/M_spectrum/spectrum.pdf")
savefig("figs/M_spectrum/spectrum.png")