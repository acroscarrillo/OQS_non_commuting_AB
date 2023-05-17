include("../../src/src.jl")

using SciMLBase
using OrdinaryDiffEq
using ProgressBars

# ["MI","MI_err","L","p","p_0","L_A","t","runs"]
df = DataFrame( CSV.File("data/MI_at_L_A_half_nb_dynamics.csv") ) 

L = 10
p_array = Vector(0:0.2:2)
p_0 = 0
runs = 100
tspan = (0.0,50.0)

dt = 1
t_array = Vector(tspan[1]:dt:tspan[end])

for (_,p) in tqdm(enumerate(p_array))
    temp = zeros(length(t_array),runs)
    for n=1:runs
        A, B = pwr_law_mat(L, p), pwr_law_mat(L, p)
        f(C,p,t) = derivative_action(C, t, A, B)
        C_0 = pwr_law_mat(L, p_0)

        prob = ODEProblem(f,C_0,tspan)

        sol = solve(prob)
        temp[:,n] .= mutual_info.(sol(t_array).u)
    end
    data_array = hcat(    mean(temp,dims=2), std(temp,dims=2)/sqrt(runs), L*ones(length(t_array)), p*ones(length(t_array)), p_0*ones(length(t_array)),  (L÷2)*ones(length(t_array)), t_array, runs*ones(length(t_array))   )
    df_2_append = DataFrame(data_array, ["MI","MI_err","L","p","p_0","L_A","t","runs"])
    df = vcat(df, df_2_append)
end

CSV.write("data/MI_at_L_A_half_nb_dynamics.csv", df)