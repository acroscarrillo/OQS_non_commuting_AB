using LinearAlgebra
using LinearSolve
using SparseArrays
using Statistics
using Metal
using ProgressBars
using DataFrames
using CSV
using Plots

function A(L::Int, p::Real)
    temp = zeros(Float64,(L,L))
    for col = 1:L
        for row = 1:L
            temp[row, col] = randn()/(abs(row-col)+1)^p
        end
    end
    temp_sym = temp * transpose(temp)
    return  temp_sym/max(eigvals(temp_sym)...)
end 

function M(h,A,B)
    h_dim, A_dim, B_dim = size(h)[1], size(A)[1], size(B)[1]
    if !(h_dim == A_dim == B_dim)
        throw(ArgumentError("Dimensions of h, A and B mismatch"))
    end
    Id = I(h_dim)
    left = Id - im*h - 0.5*transpose( A - B )
    right = Id + im*h - 0.5*( A - B )
    return kron(left, sparse(Id)) + kron(sparse(Id), right) #sparse
    # return kron(left, Id) + kron(Id, right) #dense
end

function C_NESS(h,A,B,gpu=false)
    if !gpu
        prob = LinearProblem( M(h,A,B) , ComplexF64.(vec(A)))
        sol = solve(prob)
        return reshape(sol.u, size(h))
    else
        A_f = lu(MtlArray(Float32.(M(h,A,B))))
        sol = Matrix(A_f.U) \ ( Matrix(A_f.L) \ vec(A) )
        return reshape(sol, size(h))
    end
end

function C_sub(C_tot,subsys)
    L_sub = size(subsys)[1]
    C_sub = zeros(ComplexF64,(L_sub,L_sub))
    for col = 1:L_sub
        for row = 1:L_sub
            C_sub[row,col] = C_tot[subsys[row],subsys[col]]
        end
    end
    return C_sub
end

function S(C)
    lambs = real(eigvals(C))
    entropy = 0
    for i=1:length(lambs)
        lamb = lambs[i]
        if 0<lamb<1
            entropy += -(1 - lamb)*log2(1 - lamb) - lamb*log2(lamb) 
        end
    end
    return entropy
    # return real( sum( -(1 .- lambs).*log2.(1 .- lambs) .- lambs.*log2.(lambs) ) ) THIS DOESNT WORK WHEN AN EIGVAL IS 0 OR 1 HENCE THE LOOP
end

function MI_NESS(h,A,B,subsys_A,subsys_B=nothing)
    C_tot = C_NESS(h,A,B)
    L = size(h)[1]
    L_A = length(subsys_A)
    if subsys_B==nothing
        subsys_B = zeros(L-L_A)
        counter = 1
        for i=1:L
            aux = 0
            for j=1:L_A
                if i!=subsys_A[j]
                    aux += 1
                end
            end
            if aux == L_A
                subsys_B[counter] = i
                counter += 1
            end
        end
        return S(C_sub(C_tot,subsys_A)) + S(C_sub(C_tot,Int.(subsys_B))) - S(C_tot)
    else
        tot_sub_sys = vcat(subsys_A,subsys_B)
        return return S(C_sub(C_tot,subsys_A)) + S(C_sub(C_tot,subsys_B)) - S(C_sub(C_tot,tot_sub_sys))
    end
end

function temp_bulk(C)
    L = size(C)[1]
    xi = real(C[L÷2,L÷2]) #it is a real number already as C is hermitian. This makes temp a Float instead of ComplexFloat
    return log(1/xi-1)
end

function fss_cost(params, df_in::DataFrame; g_noise=false)
    p_c, nu = params[1], params[2]

    df = DataFrame()
    if g_noise # add gaussian noise accordingly
        df.y = (df_in.MI .+ sqrt.(df_in.MI_err).*randn(size(df_in.MI))) .* (df_in.L.^(1/nu))    
        df.x = (df_in.p .- p_c) .* (df_in.L.^(1/nu))
        df.d = df_in.MI_err.*sqrt.(df_in.L)
    else # no noise
        df.y = df_in.MI.*(df_in.L.^(1/nu))    
        df.x = (df_in.p .- p_c) .* (df_in.L.^(1/nu))
        df.d = df_in.MI_err.*sqrt.(df_in.L)
    end
    sort!(df,:x)  # sort in ascending x_i

    O_val = 0
    for i=2:nrow(df)-1 # as each loop requires n.n.
        y_m, y, y_p = df.y[i-1], df.y[i], df.y[i+1]
        x_m, x, x_p = df.x[i-1], df.x[i], df.x[i+1]
        d_m, d, d_p = df.d[i-1], df.d[i], df.d[i+1]
        
        ###################
        #CAREFUL WITH THIS#   this just handles the p=p_c situation which although 
        # mathematically defined, it is numerically unstable so we put it by hand
        ###################
        temp = x_p-x_m
        if temp==0
            frac_p, frac_m = 1/2, 1/2
        else
            frac_p = (x_p-x)/temp
            frac_m = (x-x_m)/(x_p-x_m)
        end 
        ####################
        
        y_bar = y_m*frac_p + y_p*frac_m
        Delta_sqrd = d^2 + (d_m*frac_p)^2 + (d_p*frac_m)^2
                
        O_val += (y-y_bar)^2#/Delta_sqrd
    end

    return O_val 
end