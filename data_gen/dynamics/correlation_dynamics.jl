include("../../src/src.jl")

using SciMLBase
using OrdinaryDiffEq
using Plots
using ProgressBar

# ["MI","MI_err","L","p","p_0","L_A","t"]
df = DataFrame( CSV.File("data/MI_at_L_A_half_nb_dynamics.csv") ) 

L = 10
p_array = Vector(0:0.2:2)
p_0 = 5
runs = 100
tspan = (0.0,50.0)

for (j,p) in tqdm(enumerate(p_array))
    for _=1:runs
        A, B = pwr_law_mat(L, p), pwr_law_mat(L, p)
        f(C,p,t) = derivative_action(C, t, A, B)
        C_0 = pwr_law_mat(L, p_0)

        prob = ODEProblem(f,C_0,tspan)

        sol = solve(prob)
        sol.t,mutual_info.(sol.u)
    end
end

