using StatsPlots

#["tau","lamb_min","lamb_delta","tau_err","lamb_min_err","lamb_delta_err","L","p"]
df_2_plot = DataFrame( CSV.File("data/spectral_gap.csv") ) 

@df df_2_plot plot(:p, :tau ; :tau_err, markersize = 2,xlab=L"p",ylab=L"\tau",yaxis=:ln)

@df df_2_plot plot(:p, :lamb_min ; :lamb_min_err, markersize = 2,xlab=L"p",ylab=L"\lambda_{min}",yaxis=:ln)

@df df_2_plot plot(:p, :lamb_delta ; :lamb_delta_err, markersize = 2,xlab=L"p",ylab=L"\Delta \lambda")


savefig("figs/MI_vs_t/MI_dynamics_p_0_0.png")