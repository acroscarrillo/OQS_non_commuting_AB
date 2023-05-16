using ProgressBars
using Statistics
using LaTeXStrings

L_array = Vector(10:1000:10000)
p_array = [0]
runs = 10

norm_array = zeros(length(L_array)*length(p_array),4) #avg_norm, avg_norm_err, L, p

counter = 1
for (j,p) in ProgressBar(enumerate(p_array))
    for (k,L) in ProgressBar(enumerate(L_array))
        temp_array = zeros(runs)
        for n=1:runs
            A, B = randn(L,L), randn(L,L)
            A, B = A + A', B + B'
            temp_array[n] = norm( A*B - B*A )
        end
        norm_array[counter, :] .= mean(temp_array), std(temp_array)/sqrt(runs), L, p
        counter += 1
    end
end

df_norm = DataFrame(norm_array, ["avg_norm", "avg_norm_err", "L", "p"])

temp = copy(norm_array)
temp[:,1] = temp[:,1] ./ (temp[:,3].*temp[:,3])
temp[:,1] = temp[:,1] .* temp[:,1]
df_norm_by_L = DataFrame(temp, ["avg_norm", "avg_norm_err", "L", "p"])

@df df_norm plot(:p, :avg_norm, group=:L, seriestype=:scatter, markersize = 2, xlab = L"p", ylab= L"||[A,B]||")

deleteat!(df_norm_by_L,1)
@df df_norm_by_L plot(:L, :avg_norm, seriestype=:scatter, markersize = 2, xlab = L"L", ylab= L"||[A,B]||/L^2",xaxis=log,yaxis=log)

savefig("miscelaneous/other/norm_randn.png")

