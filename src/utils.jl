using DataFrames
using CSV
using Plots

function rand_sym(N::Int)
    temp = randn(N,N)
    return temp*transpose(temp)
end

function rand_sym_pos(N::Int)
    temp = randn(N,N)
    temp_sym =  temp*transpose(temp)
    return temp_sym/max(eigvals(temp_sym)...)
end