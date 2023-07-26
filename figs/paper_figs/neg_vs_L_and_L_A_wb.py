# import model.py methods (A, C, MI, ...) 
# which contains some main packages (numpy, numba,..)
import sys
import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit
import string


# os.environ["PATH"] += os.pathsep + "C:/Users/Meli/AppData/Local/Programs/MiKTeX/miktex/bin/x64"

dir = os.path.dirname(__file__)
sys.path.insert(0, os.path.abspath(os.path.join(dir,'..',"..")))
from src import *
import pickle #to save python objects

# import plotting stuff
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, NullLocator, ScalarFormatter, LogFormatter
import matplotlib.ticker

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
plt.rcParams['figure.dpi'] = 450
plt.rcParams['savefig.dpi'] = 600
plt.rcParams["text.usetex"] = True
pt = 0.0138889 

# fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(246*pt,110*pt),constrained_layout=True)

fig = plt.figure(figsize = (246*pt,110*pt),layout="compressed")
gs = fig.add_gridspec(1, 2, wspace=0.6,left=0.15,right=.95,bottom=0.3,top=0.9)
ax = gs.subplots()

# fig.subplots_adjust(right=1.01,left = 0)

labels = ["a)","b)"]
for n, a in enumerate(ax):
    a.text(-0.45, 0.85, labels[n], transform=a.transAxes, 
            size=10, weight='bold')

# fig = plt.figure(figsize = (246*pt,110*pt))
# gs = fig.add_gridspec(1,2)
# ax[0] = fig.add_subplot(gs[0])
# ax[1] = fig.add_subplot(gs[1])

color_ls = ["C0","C1","C2","C3","C4","C5","C6","C7","C8","C9"]

# MI vs L
p_array = df_L["p"].unique()
for i,p in enumerate(p_array):
    df_temp = df_L[df_L["p"] == p]

    L_array, neg_array, neg_err_array= df_temp["L"], df_temp["neg"], df_temp["neg_err"]

    g, y_intercept = fit_log_log( np.array(L_array), np.array(neg_array) )
    ax[0].plot(L_array, np.exp(y_intercept)*(L_array**g),lw=0.6,c=color_ls[i],ls="dashed")
    ax[0].errorbar(L_array,neg_array,yerr=neg_err_array,label=r'$\Delta$'+"="+str(np.round(g,2)),ms=0.8,marker="o", fmt=".",lw=0.6,color=color_ls[i])

ax[0].legend(fontsize=5,loc="upper left",framealpha=0.5)
ax[0].set_yscale("log")
ax[0].set_xscale("log")
ax[0].set_ylabel(r"$\mathcal{E}$",rotation=0)
# ax[0].set_yticks([10,0.001])
# ax[0].yaxis.set_label_coords(-0.2,0.6)
# ax[0].xaxis.set_label_coords(-0.2,-0.15)
plt.text(-0.1, 0, '(a)')
ax[0].set_xlabel(r"$L$")

# MI vs L_A
p_array = df_L_A["p"].unique()
for i,p in enumerate(p_array):
    df_temp = df_L_A[df_L_A["p"] == p]
    print(df_temp)
    LA_array, neg_array, neg_err_array =  df_temp["L_A"], df_temp["neg"], df_temp["neg_err"]

    g, y_intercept = fit_log_log( np.array(LA_array[1:len(LA_array)]),np.array(neg_array[1:len(LA_array)]) )
    # ax[1].plot(LA_array[1:len(LA_array)//2+1], np.exp(y_intercept)*(LA_array[1:len(LA_array)//2+1]**g),lw=0.6,c=color_ls[i],ls="dashed",label=r"$p=$"+str(round(p,2))+", "+ r'$\Delta$'+"="+str(np.round(g,2)))
    ax[1].plot(LA_array, np.exp(y_intercept)*(LA_array**g),lw=0.6,c=color_ls[i],ls="dashed",label=r"$p=$"+str(round(p,2))+", "+ r'$\Delta$'+"="+str(np.round(g,2)))
    # ax[1].errorbar(LA_array[1:len(LA_array)//2+1],neg_array[1:len(LA_array)//2+1],fmt=".",yerr=neg_err_array[1:len(LA_array)//2+1],c=color_ls[i],alpha=0.4,ms=0.8,marker="o",lw=1)
    ax[1].errorbar(LA_array,neg_array,fmt=".",yerr=neg_err_array,c=color_ls[i],alpha=1,ms=0.8,marker="o",lw=1)


ax[1].legend(fontsize=5,loc="upper left",framealpha=0.5)
ax[1].set_xscale("log")
ax[1].set_yscale("log")
ax[1].set_ylabel(r"$\mathcal{E}$",rotation=0)
ax[1].set_xlabel(r"$L_A$")
# ax[1].xaxis.set_minor_formatter(LogFormatter(),minor_thresholds=None)
# minor_thresholds=(2, 0.5)
# ax[1].yaxis.set_minor_formatter(LogFormatter(),minor_thresholds=None)

# ax[1].xaxis.set_major_formatter(LogFormatter(labelOnlyBase=False))
# ax[1].yaxis.set_major_formatter(LogFormatter(labelOnlyBase=True))


locmaj = matplotlib.ticker.LogLocator(base=10,numticks=12) 
ax[1].xaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=12)
ax[1].xaxis.set_minor_locator(locmin)
ax[1].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

locmaj = matplotlib.ticker.LogLocator(base=10,numticks=12) 
ax[1].yaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 ),numticks=12)
ax[1].yaxis.set_minor_locator(locmin)
ax[1].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

# ax[1].yaxis.set_minor_formatter(NullFormatter())
#ax[1].xaxis.set_major_locator(NullLocator())
#ax[1].xaxis.set_minor_locator(NullLocator())


# ax[1].yaxis.set_ticks([10**(-2),30],[r"$10^{-2}$", r"$3\times 10^1$"])
# ax[1].xaxis.set_ticks([10.0,100, 300], [r"$10$", "", r"$3\times 10^2$"])

# ax[1].yaxis.set_label_coords(-0.15,0.6)
# ax[1].xaxis.set_label_coords(0.5,-0.15)
# plt.text(-0.8, 3.3, '(b)')
plt.text(0, 10, '(b)')
# fig.tight_layout()


fig.subplots_adjust(wspace=100, hspace=1110)

# send it!
plt.show()
plt.savefig("figs/paper_figs/negativity_comparison.pdf")
plt.savefig("figs/paper_figs/negativity_comparison.png")