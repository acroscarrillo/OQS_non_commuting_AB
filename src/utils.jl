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

function sort_2_plot(df, x, y, y_err, fixed)
    data_2_plot = Vector(4) 
    for label in fixed
        df_tmp = filter(row -> row.label == "Arizona", df3)
    end
    return x_2_plot, y_2_plot, labels_2_plot
end