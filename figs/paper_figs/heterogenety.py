# import model.py methods (A, C, MI, ...) 
# which contains some main packages (numpy, numba,..)
import sys
import os

dir = os.path.dirname(__file__)
sys.path.insert(0, os.path.abspath(os.path.join(dir,'..',"..")))
import pickle #to save python objects
import numpy as np
from numpy import genfromtxt


# import plotting stuff
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, NullLocator, ScalarFormatter, LogFormatter

plt.rcParams["font.family"] = "Times New Roman"
# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=10)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rcParams['text.usetex'] = True
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 600
plt.rcParams["text.usetex"] = True


runs = 1000
L = 1000

#SELECT WHICH DATA TO PLOT#
hists = genfromtxt('data/bias_bulk.csv', delimiter=',')[1:,:]

beta_area_law_array = np.zeros(len(hists)//2)
beta_vol_law_array = np.zeros(len(hists)//2)
area_counter, volume_counter = 0, 0
for row in hists:
    if row[3] == 4.0:
        beta_area_law_array[area_counter] = row[0]
        area_counter += 1
    elif row[3] == 0.4:
        beta_vol_law_array[volume_counter] = row[0]
        volume_counter += 1


#############
# plot data #
#############
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 600
pt = 0.0138889 
fig,ax = plt.subplots(1,1,figsize = (246*pt,125*pt))

ax.hist(beta_area_law_array,label=r"$p=4$",color="C0",bins="auto",density=True,alpha=0.8)
ax.hist(beta_vol_law_array,label=r"$p=0.4$",color="C1",bins="auto",density=True,alpha=0.8)

ax.legend(loc="upper left",fontsize=7,framealpha=0.4)

ax.set_ylabel(r"$f(\sigma)$",rotation=0,fontsize=11)
ax.set_xlabel(r"$\sigma$",fontsize=12)
ax.set_ylim([-1,55])
ax.set_yticks([0.0, 50.0])
ax.set_xticks([0.0, 1.0, 2.0])
ax.set_xticklabels(['0.0', '1.0', '2.0'])
ax.set_yticklabels(['0.0', '50'])
ax.yaxis.set_label_coords(-0.06,0.45)
ax.xaxis.set_label_coords(0.525,-0.08)

######################
# standard dev inset #
# ######################

# L_array = np.arange(50,1000,50)

# area_dev = np.zeros(len(L_array))
# vol_dev = np.zeros(len(L_array))
# counter = 0
# for L in L_array:
#     name_2_load = "data/temps/temp_area_4_L_"+str(L)
#     area_array = np.load(name_2_load + ".npy")
#     area_dev[counter] = np.std(area_array)#/(np.sum(area_array)/len(area_array))
    
#     name_2_load = "data/temps/temp_vol_04_L_"+str(L)
#     vol_array = np.load(name_2_load + ".npy")
#     vol_dev[counter] = np.std(vol_array)#/(np.sum(vol_array)/len(vol_array))
    
#     counter += 1


# These are in unitless percentages of the figure size. (0,0 is bottom left)
left, bottom, width, height = [0.45, 0.4525, 0.49, 0.46]
ax_inset = fig.add_axes([left, bottom, width, height])

std_vs_L = genfromtxt('data/bias_bulk_vs_L_STD.csv', delimiter=',')[1:,:]

std_area = np.zeros(len(std_vs_L)//2)
std_vol = np.zeros(len(std_vs_L)//2)
L_array = np.zeros(len(std_vs_L)//2)
area_counter, volume_counter = 0, 0
for row in std_vs_L:
    if row[2] == 4.0:
        std_area[area_counter] = row[0]
        L_array[area_counter] = row[1]
        area_counter += 1
    elif row[2] == 0.4:
        std_vol[volume_counter] = row[0]
        volume_counter += 1

# g_vol, y_intercept_vol = fit_log_log( L_array, vol_dev )

ax_inset.plot(L_array, std_area, label="Area",lw=1,color="C0")
ax_inset.plot(L_array, std_vol, label="Volume",lw=1,color="C1")

ax_inset.set_ylabel(r"$\overline{\sigma}$",size=10,rotation=0)
ax_inset.set_xlabel(r"$L$",size=10)

ax_inset.set_xticks([100,1000])
ax_inset.set_yticks([10,0.1])

ax_inset.set_yscale("log")
ax_inset.set_xscale("log")

# g_area, y_intercept_area = fit_log_log( L_array, area_dev )
# ax_inset.plot(L_array, area_dev, label="Area",lw=1,color="C0")

# print(g_vol, y_intercept_vol)
# print(g_area, y_intercept_area)

ax_inset.set_facecolor('white')

ax_inset.xaxis.set_tick_params(labelsize=8)
ax_inset.yaxis.set_tick_params(labelsize=8)


# ax_inset.set_xticklabels(ax_inset.get_xticks())


# # ax_inset.yaxis.set_minor_formatter(NullFormatter())
# # ax_inset.set_xticks([100,1000])
ax_inset.yaxis.set_label_coords(-0.13,0.4)
ax_inset.xaxis.set_label_coords(0.6,-0.13)
# ax_inset.set_xlim([47,1000])

plt.tight_layout()
plt.savefig("figs/paper_figs/heterogenety.pdf")
plt.savefig("figs/paper_figs/heterogenety.png")
plt.show()