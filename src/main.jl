using LinearAlgebra
using Metal
using DataFrames
using Statistics
using Arpack
using SparseArrays

⨷(a, b) = kron(a, b)

function pwr_law_mat(L, p,isreal=true)
    if isreal
        temp = zeros(Float64, (L, L))
        for col in 1:L
            for row in 1:L
                temp[row, col] = randn() / (abs(row - col) + 1)^p
            end
        end
        temp_sym = temp * (temp)'
        max_lamb = eigs(temp_sym, nev = 1, ritzvec=false, which=:LR)[1][1]::Float64
        return temp_sym /max_lamb
    else
        temp = zeros(ComplexF64, (L, L))
        for col in 1:L
            for row in 1:L
                temp[row, col] = (randn()+im*randn()) / (abs(row - col) + 1)^p
            end
        end
        temp_sym = temp * (temp)'
        max_lamb = real(eigs(temp_sym, nev = 1, ritzvec=false, which=:LR)[1][1])
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
    lamb, vec_array = eigen(Hermitian(theta))
    A_tilda = (vec_array)' * (A*vec_array)
    lamb_plus_lamb = lamb .+ lamb'
    C_NESS = A_tilda ./ lamb_plus_lamb
    C_0_tilda = (vec_array)' * (C_0*vec_array)
    exp_t = exp.(-lamb_plus_lamb*t)
    C = (C_0_tilda .- C_NESS) .* exp_t .+ C_NESS
    return vec_array * C * vec_array'
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
    return mutual_info(C, subsys_A, Int.(subsys_B))
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

function bulk_bias_variance(C,bulk_size)
    L = size(C)[1]
    xi = real(diag(C))[L÷2-bulk_size÷2:L÷2+bulk_size÷2] #it is a real number already as C is hermitian. This makes temp a Float instead of ComplexFloat
    bias = log.(1 ./ xi .- 1)
    return std(bias), mean(bias)
end


function perturb_df_MI(df_in)
    df = DataFrame()
    df.MI = df_in.MI .+ sqrt.(df_in.MI_err) .* randn(size(df_in.MI))
    df.MI_err, df.L, df.p, df.L_A = df_in.MI_err, df_in.L, df_in.p, df_in.L_A
    return df
end


function perturb_df_neg(df_in)
    df = DataFrame()
    df.neg = df_in.neg .+ sqrt.(df_in.neg_err) .* randn(size(df_in.neg))
    df.neg_err, df.L, df.p, df.L_A = df_in.neg_err, df_in.L, df_in.p, df_in.L_A
    return df
end


function fss_cost_MI(params, df_in::DataFrame)
    p_c, nu = params[1], params[2]

    df = DataFrame()
    df.y = df_in.MI .* (df_in.L .^ (1 / nu))
    df.x = (df_in.p .- p_c) .* (df_in.L .^ (1 / nu))
    df.d = df_in.MI_err .* sqrt.(df_in.L)

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

        O_val += (y-y_bar)^2#/Delta_sqrd
    end

    return O_val
end




function fss_cost_neg(params, df_in::DataFrame)
    p_c, nu = params[1], params[2]

    df = DataFrame()
    df.y = df_in.neg .* (df_in.L .^ (1 / nu))
    df.x = (df_in.p .- p_c) .* (df_in.L .^ (1 / nu))
    df.d = df_in.neg_err .* sqrt.(df_in.L)

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

        O_val += (y-y_bar)^2#/Delta_sqrd
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
    M_spec = vec(lambs .+ lambs')
    return  sort(M_spec)
end





function M_spectrum(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    A_dim, B_dim = size(A)[1], size(B)[1]
    if !(A_dim == B_dim)
        throw(ArgumentError("Dimensions of A and B mismatch"))
    end
    
    L = size(A)[1]
    theta = (A + B) / 2
    lambs = eigvals(Hermitian(theta))
    M_spec = vec(lambs .+ lambs')
    return  sort(M_spec)
end


function negativity(C, subsys_A, subsys_B)
    L, L_A, L_B = size(C)[1], length(subsys_A), length(subsys_B)
    if !(L_A + L_B == L)
        throw(ArgumentError("L_A + L_B must be equal to L. Please provide an LxL Correlation matrix such that this condition is fulfiled. If you wish to work with subsystems, please provide an already traced out correlation matrix."))
    end
    G = 2*C - I(L)
    L_trans =  blockdiag(-im*sparse(I(L_A)), sparse(I(L_B)))
    R_trans = blockdiag( im*sparse(I(L_A)), -sparse(I(L_B)))
    G_plus = Matrix(L_trans * G * R_trans)
    G_min = Matrix(L_trans' * G * R_trans')
    G_T = 0.5*( I(L)-inv(I(L) + G_plus * G_min)*(G_plus+G_min) )
    mu = eigvals(G_T)
    lamb = eigvals(C)
    negativity = 0
    for i=1:L
        mu_part = log( sqrt(mu[i]) + sqrt(1-mu[i]) )
        lamb_part = 0.5*log(lamb[i]^2 + (1-lamb[i])^2)
        negativity +=  mu_part + lamb_part
    end
    return real(negativity)
end



function imbalance(C)
    L = size(C)[1]
    N_odd, N_even = 0, 0
    for i=1:L÷2
        N_even += C[2*i,2*i]
        N_odd += C[2*i-1,2*i-1]
    end
    return abs(N_odd - N_even)/(N_odd + N_even)
end