
df_2_plot = DataFrame( CSV.File("data/MI_at_L_A_half_nb_dynamics.csv") ) 

@df df_2_plot plot(:t, :MI ; :MI_err, group=:p, markersize = 2,xlab=L"t",ylab=L"\mathcal{I}",leg_title=L"p",leg_pos="bottom",legend=:bottomright)

df_formatted = filter(row -> row.p_0 <2, df_2_plot)

@df df_formatted plot(:t, :MI ; :MI_err, group=:p, markersize = 2,xlab=L"t",ylab=L"\mathcal{I}",leg_title=L"p",leg_pos="bottom",legend=:bottomright)

savefig("figs/MI_vs_t/MI_dynamics_p_0_0.png")