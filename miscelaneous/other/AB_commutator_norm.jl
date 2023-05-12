using ProgressBars
using Statistics
using LaTeXStrings

L_array = Vector(10:100:1000)
p_array = Vector(0:0.2:10)
runs = 100

norm_array = zeros(length(L_array)*length(p_array),4) #avg_norm, avg_norm_err, L, p

counter = 1
for (j,p) in ProgressBar(enumerate(p_array))
    for (k,L) in enumerate(L_array)
        temp_array = zeros(runs)
        for n=1:runs
            A, B = pwr_law_mat(L,p), pwr_law_mat(L,p)
            temp_array[n] = norm( A*B - B*A )
        end
        norm_array[counter, :] .= mean(temp_array), std(temp_array)/sqrt(L), L, p
        counter += 1
    end
end

df_norm = DataFrame(norm_array, ["avg_norm", "avg_norm_err", "L", "p"])

temp = norm_array
temp[:,1] = temp[:,1] ./ temp[:,3]
df_norm_by_L = DataFrame(temp, ["avg_norm", "avg_norm_err", "L", "p"])

@df df_norm plot(:p, :avg_norm, group=:L, seriestype=:scatter, markersize = 2, xlab = L"p", ylab= L"||[A,B]||")

@df df_norm_by_L plot(:p, :avg_norm, group=:L, seriestype=:scatter, markersize = 2, xlab = L"p", ylab= L"||[A,B]||/L")

savefig("miscelaneous/other/norm.png")

