using LaTeXStrings

df_in = DataFrame(CSV.File("data/M_spectrum.csv")) 

L = Int( unique(df_in,:L).L[1] ) #This is obviously not general
runs = Int( maximum(  unique(df_in,:sample_n).sample_n  ) )

df_level_spacing = DataFrame(zeros(0,3), ["spacing","L","p"] )
for n=tqdm(1:runs)
    df_temp = filter(row -> row.sample_n == n, df_in)
    p_array = unique(df_in,:p).p
    for p in p_array
        lambs = filter(row -> row.sample_n == n, df_in).lamb
        level_spacings = zeros(L^2-1)
        for i=1:L^2-1
            level_spacings[i] = (lambs[i+1]-lambs[i])
        end

        data = hcat( level_spacings, L*ones(L^2-1), p*ones(L^2-1))

        df_2_append = DataFrame(data,["spacing","L","p"] )
        df_level_spacing = vcat(df_level_spacing, df_2_append)
    end
end

data = vcat(filter(row -> row.p == 4, df_level_spacing).spacing',filter(row -> row.p == .4, df_level_spacing).spacing')

histogram(data',ylabel=L"f(\lambda_M)",xlabel=L"s",label = [L"p=4" L"p=0.4"],fillalpha=.4,normalize=true)

#savefig("figs/M_spectrum/spectrum.pdf")
#savefig("figs/M_spectrum/spectrum.png")