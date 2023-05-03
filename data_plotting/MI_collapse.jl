include("../src/src.jl")
using LsqFit

model(x, p) = p[1]*exp.(-x*p[2])
xdata = range(0, stop=10, length=20)
ydata = model(xdata, [1.0 2.0]) + 0.01*randn(length(xdata))
p0 = [0.5, 0.5]
fit = curve_fit(model, xdata, ydata, p0,lower=[0.0,0.0])
fit.param

df = DataFrame(CSV.File("data/MI_at_L_A_half_wb_zoom.csv")) #format MI|MI_err|L|p|L_A

##############################
# call fitter, generate data #
##############################
# res = minimize(O_w_b_alt, np.array([1,1,1]), method='Powell',options={'disp': True, "maxiter": 1000}, bounds=((0,2),(0,1),(0,2)),tol=1e-20)

# # p_c, nu, zeta = 0.97,5.2,0# res.x
# p_c, nu, zeta = res.x

# print(p_c,nu,zeta)

data_in=np.load("data/MI_vs_p_with_boundaries_1to13range.npy") #format MI|MI_err|L|p
data  = np.zeros((len(data_in),4), dtype=np.float64)
for i in range(len(data_in)):
    data[i,0] = (data_in[i,2]**(zeta/nu))*data_in[i,0]                 # y_i
    data[i,1] = ( data_in[i,3] - p_c )*(data_in[i,2]**(1/nu))          # x_i
    data[i,2] = data_in[i,1]                                           # d_i
    data[i,3] = data_in[i,2]                                           # L

data_sorted = data[data[:, 1].argsort()] # sort in ascending x_i

L_list = [400,800,1200,1600,2000]
data_2_plot = []
for L in L_list:
    aux_array = np.zeros((int(len(data_sorted)/len(L_list)),3))
    counter = 0
    for d in data_sorted:
        if d[3] == L:
            aux_array[counter,0] = d[0]
            aux_array[counter,1] = d[1]
            aux_array[counter,2] = d[2]
            counter += 1
    data_2_plot.append(aux_array)

with open("data/collapse_data_w_b_1to13range.ob", 'wb') as fp:
    pickle.dump(data_2_plot, fp)

with open ("data/collapse_data_w_b_1to13range.ob", 'rb') as fp:
    data_2_plot = pickle.load(fp)


res = minimize(O_w_b_alt, np.array([1,1,1]), method='Powell',options={'disp': True, "maxiter": 1000}, bounds=((0,2),(0,1),(0,2)),tol=1e-20)

# p_c, nu, zeta = 0.97,5.2,0# res.x
p_c, nu, zeta = res.x

print(p_c,nu,zeta)

#############
# plot data #
#############
pt = 0.0138889 
fig, ax  = plt.subplots(figsize = (246*pt,200*pt))

plt.rcParams['axes.facecolor'] = 'white'
fig.set_facecolor("white")
ax.set_facecolor('white')
plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.dpi'] = 200

for i in range(len(L_list)):
    ax.scatter(data_2_plot[i][:,1],data_2_plot[i][:,0],label=r"$L=$"+str(L_list[i]),s=3)

# ax.set_yscale('log')
ax.set_ylabel(r"$\mathcal{I}$")
ax.set_xlabel(r"$L^{1/\nu}(p-p_c)$")
plt.title(r"$p_c=$"+str(p_c)+", "+r"$\nu=$"+str(nu),size=10)
ax.legend(fontsize=8);
plt.tight_layout()

plt.show()