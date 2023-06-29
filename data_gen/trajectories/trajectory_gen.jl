include("../../src/main_trajectories.jl")

using Plots
using StatsBase
using ProgressBars

function H_eff_evolve_to_t(psi_0, theta_array, c, t)
    N = length(psi_0)
    psi_t = (I(N) + c[1]' * c[1] *(exp(-t * theta_array[1] / 2) - 1)) * psi_0
    for i in 2:L
        psi_t = (I(N) + c[i]' * c[i] *(exp(-t * theta_array[i] / 2) - 1)) * psi_t
    end
    return psi_t
end


function jump(psi, rates_array, jump_op_array)
    weights = [rates_array[i] * norm(jump_op_array[i]*psi) for i=1:length(rates_array)]
    jump_op = sample(jump_op_array, Weights(weights))
    return jump_op*psi/norm(jump_op*psi)
end

L = 8
p = 1
t_0, dt, t_f = 0, .01, 20

t_array = Vector(t_0:dt:t_f)
rand_numbs = rand(length(t_array))

A, B = pwr_law_mat(L, p), pwr_law_mat(L, p)
Omega = A + B   
theta_array = eigvals(Omega)

c = [ComplexF64.(create(L, j)) for j in 1:L]
jump_op_array = vcat(c,Transpose.(c))
rates_array = vcat(theta_array,theta_array)

# psi_0 = randn(2^L) + im * randn(2^L)
# psi_0 = pure_rand_gaussian(L,5,c)
psi_0 = zeros(2^L)
psi_0[rand(1:2^L)] = 1
psi_0 = psi_0 / norm(psi_0)

psi_t = [ComplexF64.(zeros(2^L)) for t in t_array]
psi_t[1] = psi_0
norm_t = zeros(length(t_array))


jumps_counter = 1
for j in tqdm(2:length(t_array))
    # display("norm="*string(norm(psi_t[j-1]))*" ,η="*string(rand_numbs[jumps_counter])*", jumps="*string(jumps_counter)*", t="*string(t_array[j]))
    ocupation_t = [psi_t[j]'*c[i]'*c[i]*psi_t[j] for i=1:L]
    norm_t[j-1] = norm(psi_t[j-1])
    if norm_t[j-1]^2 >= rand_numbs[jumps_counter]
        psi_t[j] = H_eff_evolve_to_t(psi_t[j-1], theta_array, c, dt)
    else
        psi_t[j] =  jump(psi_t[j-1], rates_array, jump_op_array)
        jumps_counter += 1
    end
end

psi_t = psi_t ./ norm.(psi_t)

correlation_t = [ComplexF64.(zeros(L,L)) for t in t_array]
for j in tqdm(1:length(t_array))
    temp_correlation = ComplexF64.(zeros(L,L))
    for n=1:L
        for m=1:L
            temp_correlation[n,m] = psi_t[j]' * c[n]'*c[m] * psi_t[j]
        end
    end
    correlation_t[j] = temp_correlation
end

gaussness_t = zeros(length(t_array))
for j in tqdm(1:length(t_array))
    rho_gauss = correlation_to_rho(correlation_t[j],c)
    rho_t = psi_t[j] * psi_t[j]'
    # display(rho_gauss-rho_t)
    gaussness_t[j] = 0.5*norm(rho_gauss-rho_t,2)
end


N = sum([c[i]'*c[i] for i=1:L])
N_exp = zeros(length(t_array))
for j=1:length(t_array)
    N_exp[j] = real(psi_t[j]' * N * psi_t[j])
end

plot(t_array,gaussness_t)