#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 10:09:40 2020

@author: jaime
"""
import dill
import numpy as np
import cv2 as cv
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns

import trompy as tp
# from sympy import *
# from sympy.geometry import *

def get_good_xys(df, threshold=0.9, verbose=False):
    """
    Returns xy coordinates from deeplabcut only when high likelihood.

    Parameters
    ----------
    df : Dataframe
        Includes xy coordinates and likelihood column.
    threshold : Float, optional
        Likelihood threshold for inclusion. The default is 0.9.
    verbose : Bool, optional
        Prints number of frames excluded. The default is False.

    Returns
    -------
    None.

    """
    
    L = df['likelihood'] > threshold
    F = df.index[df['likelihood'] > threshold].tolist()
    df_hi = df[L]
    
    x = [int(xval) for xval in df_hi['x']]
    y = [int(yval) for yval in df_hi['y']]
    if verbose:
        print(f"There are {sum(L)} good frames out of a total of {len(df)}.")
    
    return (x, y, F)

def plot_cues(ax, xydict, scale=8):
    casx, casy = [coord/scale for coord in xydict['cas_pos']]
    maltx, malty = [coord/scale for coord in xydict['malt_pos']]
    
    ax.plot(casx, casy, marker='+', color='white')
    ax.plot(maltx, malty, marker='+', color='grey')



def calc_dist(pt1, pt2):
    return np.sqrt((pt1[0] - pt2[0])**2 + (pt1[1] - pt2[1])**2)


MASTERFOLDER = "/home/jaime/Github/PPP_video/"

pickle_in = MASTERFOLDER + "PPP_video_data.pickle"
with open(pickle_in, "rb") as dill_file:
    PPP_video_data = dill.load(dill_file)

plot_heat_maps = True

if plot_heat_maps:
    f1, ax = plt.subplots(ncols=14, figsize=(16, 2))
    for idx, data in enumerate(PPP_video_data):
        (x, y, F) = get_good_xys(data['nose'], threshold=0.95)
        
        axis = ax[idx]
        sns.kdeplot(x, y, cmap="Reds", shade=True, bw = 25, shade_lowest=False, n_levels = 50, ax=axis)
        # ax = sns.kdeplot(x, y, kernel="gau", bw = 25, cmap="Reds", n_levels = 50, shade=True, shade_lowest=False, gridsize=100)
        axis.set_frame_on(False)
        plot_cues(axis, data, scale=1)
        axis.set_xlim(0, 640)
        axis.set_ylim(480, 0)
        axis.set_xticks([])
        axis.set_yticks([])
        for pos in ['left', 'right', 'top', 'bottom']:
            axis.spines[pos].set_visible(True)
        
        axis.set_title(data['rat'])
        # axis.axis('off')
        # plt.show()
    
    f1.savefig(MASTERFOLDER + "heatmaps.pdf")

# Initializes dictionaries to store data relating to casein and maltodextrin
cas_dist_med, malt_dist_med = {}, {}

for d in [cas_dist_med, malt_dist_med]:
    d['NR'] = []
    d['PR'] = []
    
for data in PPP_video_data:
    
    diet = data['diet']
    
    (x, y, F) = get_good_xys(data['nose'], threshold=0.95)
    cas_dist, malt_dist = [], []
    for X, Y in zip(x, y):
        cas_dist.append(calc_dist((X, Y), data['cas_pos']))
        malt_dist.append(calc_dist((X, Y), data['malt_pos']))
    
    cas_dist_med[diet].append(np.median(cas_dist))
    malt_dist_med[diet].append(np.median(malt_dist))

NR = [cas_dist_med['NR'], malt_dist_med['NR']]
PR = [cas_dist_med['PR'], malt_dist_med['PR']]

fig, ax = plt.subplots()
tp.barscatter(tp.data2obj2D([NR, PR]), ax=ax, paired=True,
              barfacecoloroption="individual",
              barfacecolor=["grey", "w", "blue", "xkcd:light blue"],
              barlabels=['Cas', 'Malt', 'Cas', 'Malt'],
              grouplabel=['NR', 'PR'],
              grouplabeloffset=-0.03
              )

ax.set_ylabel("Distance from sipper (px)")

# fig, ax = plt.subplots(ncols=2)
# fig.subplots_adjust(wspace=0.05)

# bins = range(0, 500, 20)

# cas_hist = np.histogram(cas_dist, bins=bins)
# malt_hist = np.histogram(malt_dist, bins=bins)
# # malt_bins = [-b for b in bins]
# ax[0].set_xlim([500, 0])
# ax[1].set_xlim([0, 500])
# ax[0].plot(bins[:-1], cas_hist[0])
# ax[1].plot(bins[:-1], malt_hist[0])

# ax[1].set_yticks([])
# ax[1].spines['left'].set_visible(False)

# for axis in [ax[0], ax[1]]:
#     axis.set_ylim([-20, 4200])
#     axis.spines['right'].set_visible(False)
#     axis.spines['top'].set_visible(False)
    
# ax[0].set_xlabel('Distance from Casein cue (px)')
# ax[1].set_xlabel('Distance from Maltodextrin cue (px)')
# ax[0].set_ylabel('Frequency')




# # https://zbigatron.com/generating-heatmaps-from-coordinates/

