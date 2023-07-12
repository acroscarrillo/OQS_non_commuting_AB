# import model.py methods (A, C, MI, ...) 
# which contains some main packages (numpy, numba,..)
import sys
import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit

# os.environ["PATH"] += os.pathsep + "C:/Users/Meli/AppData/Local/Programs/MiKTeX/miktex/bin/x64"

dir = os.path.dirname(__file__)
sys.path.insert(0, os.path.abspath(os.path.join(dir,'..',"..")))
from src import *
import pickle #to save python objects

# import plotting stuff
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, NullLocator, ScalarFormatter, LogFormatter

from numpy import genfromtxt

plt.rcParams["font.family"] = "Times New Roman"
# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=10)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
plt.rc('legend', fontsize=10)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title



def fit_log_log(x, y, subset=None):
    """A simple function to perform a linear fit to a set of data in log-log scale.
    
    Sorts the data according to the x axis, then performs linear regression on the logged
    x and y data for the specified subset of indices. If the subset is not specified, all
    data is used.

    Args:
        x: data for the x axis
        y: data for the y axis
        subset: an object which can be used for numpy array indexing, indicating the
            subset of data sorted data to use

    Returns:
        A tuple of the gradient and y intercept
    """
    sort_permutation = np.argsort(x)
    if type(subset) == None:
        x_data = np.log(x[sort_permutation])
        y_data = np.log(y[sort_permutation])
    else:
        x_data = np.log(x[sort_permutation][subset])
        y_data = np.log(y[sort_permutation][subset])
    res = stats.linregress(x_data, y_data)
    return res[0], res[1]


###############
# import data #
###############
df_L = pd.read_csv("data/neg_at_L_A_half_wb.csv") 
df_L_A = pd.read_csv("data/neg_vs_L_A_wb.csv") 


#############
# plot data #
#############
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 600
plt.rcParams["text.usetex"] = True
pt = 0.0138889 
fig = plt.figure(figsize = (246*pt,250*pt))
gs = fig.add_gridspec(2,1)
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

color_ls = ["C0","C1","C2","C3","C4","C5","C6","C7","C8","C9"]

# MI vs L
p_array = df_L["p"].unique()
for i,p in enumerate(p_array):
    df_temp = df_L[df_L["p"] == p]

    L_array, neg_array, neg_err_array= df_temp["L"], df_temp["neg"], df_temp["neg_err"]

    g, y_intercept = fit_log_log( np.array(L_array), np.array(neg_array) )
    ax1.plot(L_array, np.exp(y_intercept)*(L_array**g),lw=0.6,c=color_ls[i],ls="dashed")
    ax1.errorbar(L_array,neg_array,yerr=neg_err_array,label=r'$\Delta$'+"="+str(np.round(g,2)),ms=0.8,marker="o", fmt=".",lw=0.6,color=color_ls[i])

ax1.legend(fontsize=6,loc="upper left",framealpha=0.5)
ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.set_ylabel(r"$\mathcal{E}$",rotation=0)
# ax1.set_yticks([10,0.001])
# ax1.yaxis.set_label_coords(-0.2,0.6)
# ax1.xaxis.set_label_coords(-0.2,-0.15)
plt.text(-2, 3.3, '(b)')
ax1.set_xlabel(r"$L$")

# MI vs L_A
p_array = df_L_A["p"].unique()
for i,p in enumerate(p_array):
    df_temp = df_L_A[df_L_A["p"] == p]
    print(df_temp)
    LA_array, neg_array, neg_err_array =  df_temp["L_A"], df_temp["neg"], df_temp["neg_err"]

    g, y_intercept = fit_log_log( np.array(LA_array[1:len(LA_array)//4+1]),np.array(neg_array[1:len(LA_array)//4+1]) )
    # ax2.plot(LA_array[1:len(LA_array)//2+1], np.exp(y_intercept)*(LA_array[1:len(LA_array)//2+1]**g),lw=0.6,c=color_ls[i],ls="dashed",label=r"$p=$"+str(round(p,2))+", "+ r'$\Delta$'+"="+str(np.round(g,2)))
    # ax2.errorbar(LA_array[1:len(LA_array)//2+1],neg_array[1:len(LA_array)//2+1],fmt=".",yerr=neg_err_array[1:len(LA_array)//2+1],c=color_ls[i],alpha=0.4,ms=0.8,marker="o",lw=1)
    ax2.errorbar(LA_array,neg_array,fmt=".",yerr=neg_err_array,c=color_ls[i],alpha=0.4,ms=0.8,marker="o",lw=1)


ax2.legend(fontsize=4,loc="lower left",framealpha=0.5)
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_ylabel(r"$\mathcal{E}$",rotation=0)
ax2.set_xlabel(r"$L_A$")
# ax2.xaxis.set_minor_formatter(LogFormatter(),minor_thresholds=None)
# minor_thresholds=(2, 0.5)
# ax2.yaxis.set_minor_formatter(LogFormatter(),minor_thresholds=None)

# ax2.xaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
# ax2.yaxis.set_major_formatter(LogFormatter(la3belOnlyBase=True))



# ax2.yaxis.set_minor_formatter(NullFormatter())
#ax2.xaxis.set_major_locator(NullLocator())
#ax2.xaxis.set_minor_locator(NullLocator())


# ax2.yaxis.set_ticks([10**(-6),20],[r"$10^{-6}$", r"$2\times 10^1$"])
# ax2.xaxis.set_ticks([10.0,100, 300], [r"$10$", "", r"$3\times 10^2$"])
# ax2.yaxis.set_label_coords(-0.15,0.6)
# ax2.xaxis.set_label_coords(0.5,-0.15)
plt.text(-0.8, 3.3, '(c)')

# send it!
plt.tight_layout()
plt.show()
plt.savefig("paper_figs/images/MI_no_boundaries_fig.pdf")