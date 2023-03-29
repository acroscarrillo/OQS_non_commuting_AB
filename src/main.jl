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
    left = I(h_dim)-im*h-0.5*transpose(A-B)
    right = I(h_dim)+im*h-0.5*(A-B)
    return kron(left,right)
end

function NESS(h,A,B)
    prob = LinearProblem( M(h,A,B) , vec(A))
    sol = solve(prob)
    return reshape(sol.u, size(h))
end

# aux_mat = np.zeros((L,L),dtype=np.complex128)
#         for i in range(L): 
#             for j in range(L):
#                 if j+i<=L-1:
#                     aux_mat[j,j+i]= ( np.random.normal(0,1) + 1j*np.random.normal(0,1) )/(i+1)**b#np.exp(i*b/L)#(i+1)**b#(i/L+1)**b
#                     aux_mat[j+i,j]= ( np.random.normal(0,1) + 1j*np.random.normal(0,1) )/(i+1)**b#np.exp(i*b/L)#(i+1)**b#(i/L+1)**b
#     #                 aux_mat[j,j+i]= np.random.normal(0,1) #ATENTION!!! MATRIX "A" IS RESTRICTED TO REAL ENTRIES, CHANGE THIS IF NEEDED!!! 
#         A = np.dot(aux_mat, aux_mat.conj().T)

#         lamb = LA.eigvalsh(A)

#         return A/(2*max(lamb))