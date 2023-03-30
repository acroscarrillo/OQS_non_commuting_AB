using LinearAlgebra
using LinearSolve


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
    return kron(left, Id) + kron(Id, right)
end

function C_NESS(h,A,B)
    prob = LinearProblem( M(h,A,B) , ComplexF64.(vec(A)))
    sol = solve(prob)
    return reshape(sol.u, size(h))
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
    lambs = eigvals(C) 
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