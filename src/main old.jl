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




function derivative_action(C, t, h, A, B)
    return im * (transpose(h) * C - C * transpose(h)) - 0.5 * ((A + B) * C + C * (A + B)) + A
end
function derivative_action(C, t, A, B)
    return - 0.5 * ((A + B) * C + C * (A + B)) + A
end



function master_op(h, A, B,use_sparse::Bool=true)
    h_dim, A_dim, B_dim = size(h)[1], size(A)[1], size(B)[1]
    if !(h_dim == A_dim == B_dim)
        throw(ArgumentError("Dimensions of h, A and B mismatch"))
    end
    Id = I(h_dim)
    left = - im * h - 0.5 * (A + B)'
    right = im * h - 0.5 * (A + B)
    if use_sparse
        return left ⨷ sparse(Id) + sparse(Id) ⨷ right #sparse
    else
        return left ⨷ Id + Id ⨷ right #dense
    end
end
function master_op(A, B,use_sparse::Bool=true)
    A_dim, B_dim = size(A)[1], size(B)[1]
    if !(A_dim == B_dim)
        throw(ArgumentError("Dimensions of A and B mismatch"))
    end
    Id = I(A_dim)
    left = - 0.5 * (A + B)'
    right = - 0.5 * (A + B)
    if use_sparse
        return left ⨷ sparse(Id) + sparse(Id) ⨷ right #sparse
    else
        return left ⨷ Id + Id ⨷ right #dense
    end
end




function correlation_steady_state(h, A, B, ness_guess::Bool)
    if ness_guess
        prob = LinearProblem(master_op(h, A, B), -ComplexF64.(vec(A));u0=-ComplexF64.(vec(A)))
        sol = solve(prob)
        return reshape(sol.u, size(A))
    else
        prob = LinearProblem(master_op(h, A, B), -ComplexF64.(vec(A)))
        sol = solve(prob)
        return reshape(sol.u, size(A))
    end
end
function correlation_steady_state(h, A, B)
    prob = LinearProblem(master_op(h, A, B), -ComplexF64.(vec(A)))
    sol = solve(prob)
    return reshape(sol.u, size(A))
end
function correlation_steady_state(A, B, gpu::Bool)
    if !gpu
        prob = LinearProblem(master_op(A, B), -ComplexF64.(vec(A)))
        sol = solve(prob)
        return reshape(sol.u, size(A))
    else
        A_f = lu(MtlArray(Float32.(master_op(A, B))))
        sol = Matrix(A_f.U) \ (Matrix(A_f.L) \ vec(A))
        return reshape(sol, size(A))
    end
end
function correlation_steady_state(A, B) #if gpu not specified, use based on L
    L = size(A)[1]
    if L < 40
        return correlation_steady_state(A, B, false)
    else
        return correlation_steady_state(A, B, true)
    end
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