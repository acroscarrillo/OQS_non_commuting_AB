
df = DataFrame( CSV.File("data/MI_at_L_A_half_nb_dynamics.csv") ) 

@df df_norm plot(:t, :MI ; :MI_err, group=:p, markersize = 2)
savefig("miscelaneous/other/norm.png")