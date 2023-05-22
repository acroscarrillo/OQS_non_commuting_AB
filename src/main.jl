using LinearAlgebra
using Metal
using DataFrames
using Arpack

⨷(a, b) = kron(a, b)

function pwr_law_mat(L, p,real=true)
    if real
        temp = zeros(Float64, (L, L))
        for col in 1:L
            for row in 1:L
                temp[row, col] = randn() / (abs(row - col) + 1)^p
            end
        end
        temp_sym = temp * (temp)'
        max_lamb = eigs(temp_sym, nev = 1, ritzvec=false, which=:LR)[1][1]
        return temp_sym /max_lamb
    else
        temp = zeros(ComplexF64, (L, L))
        for col in 1:L
            for row in 1:L
                temp[row, col] = (randn()+im*randn()) / (abs(row - col) + 1)^p
            end
        end
        temp_sym = temp * (temp)'
        max_lamb = eigs(temp_sym, nev = 1, ritzvec=false, which=:LR)[1][1]
        return temp_sym /max_lamb
    end
end



function correlation_ness(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    A_dim, B_dim = size(A)[1], size(B)[1]
    if !(A_dim == B_dim)
        throw(ArgumentError("Dimensions of A and B mismatch"))
    end
    
    L = size(A)[1]
    theta = (A + B) / 2
    lamb, vec_array = eigen(Hermitian(theta))
    A_tilda = (vec_array)' * (A*vec_array)
    C_NESS = A_tilda ./ (lamb .+ lamb')
    return vec_array * C_NESS * vec_array'
end



function correlation(A,B,C_0,t)
    A_dim, B_dim, C_0_dim = size(A)[1], size(B)[1], size(C_0)[1]
    if !(A_dim == B_dim == C_0_dim)
        throw(ArgumentError("Dimensions of A, B and C_0 mismatch"))
    end
    
    L = size(A)[1]
    theta = (A + B) / 2
    lamb, vec_array = eigen(theta)
    A_tilda = (vec_array)' * (A*vec_array)
    C_0_tilda = (vec_array)' * (C_0*vec_array)
    C = zeros(L,L)
    for n=1:L
        for m=1:L
            temp_ness = A_tilda[n,m] ./ (lamb[n] + lamb[m])
            temp_C_0 = C_0_tilda[n,m]
            exponent = exp( -(lamb[n] + lamb[m])*t )
            C += (vec_array[:, n]*vec_array[:, m]') * ( (temp_C_0-temp_ness)*exponent  +  temp_ness)
        end
    end
    return C
end




function sub_correlation(C_tot, subsys)
    L_sub = size(subsys)[1]
    C_sub = zeros(ComplexF64, (L_sub, L_sub))
    for col in 1:L_sub
        for row in 1:L_sub
            C_sub[row, col] = C_tot[subsys[row], subsys[col]]
        end
    end
    return C_sub
end



function vn_entropy(C)
    lambs = real(eigvals(C))
    entropy = 0
    for i in 1:length(lambs)
        lamb = lambs[i]
        if 0 < lamb < 1
            entropy += -(1 - lamb) * log2(1 - lamb) - lamb * log2(lamb)
        end
    end
    return entropy
    # return real( sum( -(1 .- lambs).*log2.(1 .- lambs) .- lambs.*log2.(lambs) ) ) THIS DOESNT WORK WHEN AN EIGVAL IS 0 OR 1 HENCE THE LOOP
end



function mutual_info(C, subsys_A, subsys_B)
    tot_sub_sys = vcat(subsys_A, subsys_B)
    return vn_entropy(sub_correlation(C, subsys_A)) +
           vn_entropy(sub_correlation(C, subsys_B)) -
           vn_entropy(sub_correlation(C, tot_sub_sys))
end
function mutual_info(C, subsys_A)
    L = size(C)[1]
    L_A = length(subsys_A)
    subsys_B = zeros(L - L_A)
    counter = 1
    for i in 1:L
        aux = 0
        for j in 1:L_A
            if i != subsys_A[j]
                aux += 1
            end
        end
        if aux == L_A
            subsys_B[counter] = i
            counter += 1
        end
    end
    return MI_NESS(C, subsys_A, Int.(subsys_B))
end
function mutual_info(C)
    L = size(C)[1]
    subsys_A = Vector(1:(L÷2))
    return mutual_info(C, subsys_A)
end




function central_occ_bias(C)
    L = size(C)[1]
    xi = real(C[L ÷ 2, L ÷ 2]) #it is a real number already as C is hermitian. This makes temp a Float instead of ComplexFloat
    return log(1 / xi - 1)
end




function fss_cost(params, df_in::DataFrame; g_noise=false)
    p_c, nu = params[1], params[2]

    df = DataFrame()
    if g_noise # add gaussian noise accordingly
        df.y =
            (df_in.MI .+ sqrt.(df_in.MI_err) .* randn(size(df_in.MI))) .*
            (df_in.L .^ (1 / nu))
        df.x = (df_in.p .- p_c) .* (df_in.L .^ (1 / nu))
        df.d = df_in.MI_err .* sqrt.(df_in.L)
    else # no noise
        df.y = df_in.MI .* (df_in.L .^ (1 / nu))
        df.x = (df_in.p .- p_c) .* (df_in.L .^ (1 / nu))
        df.d = df_in.MI_err .* sqrt.(df_in.L)
    end
    sort!(df, :x)  # sort in ascending x_i

    O_val = 0
    for i in 2:(nrow(df) - 1) # as each loop requires n.n.
        y_m, y, y_p = df.y[i - 1], df.y[i], df.y[i + 1]
        x_m, x, x_p = df.x[i - 1], df.x[i], df.x[i + 1]
        d_m, d, d_p = df.d[i - 1], df.d[i], df.d[i + 1]

        ###################
        #CAREFUL WITH THIS#   this just handles the p=p_c situation which although 
        # mathematically defined, it is numerically unstable so we put it by hand
        ###################
        temp = x_p - x_m
        if temp == 0
            frac_p, frac_m = 1 / 2, 1 / 2
        else
            frac_p = (x_p - x) / temp
            frac_m = (x - x_m) / temp
        end
        ####################

        y_bar = y_m * frac_p + y_p * frac_m
        Delta_sqrd = d^2 + (d_m * frac_p)^2 + (d_p * frac_m)^2

        O_val += (y - y_bar)^2#/Delta_sqrd
    end

    return O_val
end



function M_spectrum(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    A_dim, B_dim = size(A)[1], size(B)[1]
    if !(A_dim == B_dim)
        throw(ArgumentError("Dimensions of A and B mismatch"))
    end
    
    L = size(A)[1]
    theta = (A + B) / 2
    lambs = eigvals(Hermitian(theta))
    M_spec = Vector(lamb .+ lamb')
    return lamb
end

