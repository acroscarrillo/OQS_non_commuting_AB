df_in = DataFrame(CSV.File("data/M_spectrum.csv")) 


data = vcat(filter(row -> row.p == 4, df_spectrum).lamb',filter(row -> row.p == .4, df_spectrum).lamb')

data = filter(row -> row.p == 4, df_spectrum).lamb'

histogram(abs.(data'),fillalpha=4,normalize=true)

histogram2d(real.(data'),imag.(data'),ylab=L"Im(M)",xlab=L"Re(M)",fillalpha=0.4,normalize=true,bins=(range(-2, 2,100),range(-2, 2,100)),title=L"p=4" )

@df filter(row -> row.p == 4, df_in) histogram(:lamb)

@df filter(row -> row.p == .4, df_in) histogram(:lamb)

