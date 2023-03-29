using LinearAlgebra
using LinearSolve


function A(L::Int, p::Real)
    temp = zeros(Float64,(L,L))
    for col = 1:L
        for row = 1:L
            temp[row, col] = randn()/(abs(row-col)+1)^p
        end
    end
    return temp * transpose(temp)
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

function NESS(h,A,B)
    prob = LinearProblem( M(h,A,B) , vec(A))
    sol = solve(prob)
    return reshape(sol.u, size(h))
end