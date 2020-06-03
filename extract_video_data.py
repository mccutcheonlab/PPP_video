#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 10:42:33 2020

@author: jaime
"""


import tdt
import cv2 as cv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import dill

import trompy as tp


def make_ppp_video_dict(rat, diet, box, tank, cam,
                        vidfile, 
                        posfile, scorer,
                        cas_ttl, malt_ttl,
                        casxy, maltxy):
    """
    Function to make dictionary with position and cue information for rats in PPP
    experiment.

    Parameters
    ----------
    tank : String
        Path to TDT tank.
    posfile : String
        Path to hd5 file generated by DeepLabCut.
    scorer : String
        Scorer string from DeepLabCut file. Can be generated automatically in future.
    cas_ttl : String
        TTL for casein trial times from TDT file.
    malt_ttl : String
        TTL for maltodextring trial times.
    casxy : Tuple
        Position of casein cue light.
    maltxy : Tuple
        Position of maltodextrin cue light.

    Returns
    -------
    data : Dictionary
        Includes position data, cue data and other information needed for further analysis.

    """
    
    epocs = (tdt.read_block(tank, evtype=['epocs'])).epocs

    cap = cv.VideoCapture(vidfile)
    fps = int(cap.get(cv.CAP_PROP_FPS))
    success,image = cap.read()
    print(success, vidfile)
    cv.imwrite(Path(vidfile).stem + ".jpg", image)
    cap.release()

    df = pd.read_hdf(posfile)
    
    return {'rat': rat, 'diet': diet, 'box': box,
            'nose': df[scorer, 'nose'],
            'L_ear': df[scorer, 'L_ear'],
            'R_ear': df[scorer, 'R_ear'],
            'fps': fps,
            'cas_cue': getattr(epocs, cas_ttl).onset,
            'malt_cue': getattr(epocs, malt_ttl).onset,
            'cas_pos': casxy,
            'malt_pos': maltxy}

def process_row(row,
                fileprefix,
                scorer,
                pos):
    """
    

    Parameters
    ----------
    row : List
        DESCRIPTION.
    fileprefix : TYPE
        DESCRIPTION.
    scorer : TYPE
        DESCRIPTION.
    pos : List of list of tuples
        Position (x, y) of left and right sippers in each box.

    Returns
    -------
    output : TYPE
        DESCRIPTION.

    """
    
    # gets data from metafile row
    rat = row[2]
    diet = row[6]
    box = row[5][0]
    tank = TANKFOLDER + row[0]
    
    print(f"Processing rat {row[2]}...")

    # uses box to work out file with tracking data and set cue position

    if box == '1':
        cam = "Cam2"
        Lxy = pos[0][0] # (428, 115)
        Rxy = pos[0][1] # (530, 114)

    else:
        cam = "Cam1"
        Lxy = pos[1][0] # (359, 92)
        Rxy = pos[1][1] # (461, 97)
    
    # works out whether casein on left or right so that cue pos and ttls can be set
    if 'cas' in row[8]:
        cas_ttl = row[18]
        malt_ttl = row[19]
        casxy = Lxy
        maltxy = Rxy
        print(rat, 'Hoo yeah')
    else:
        cas_ttl = row[19]
        malt_ttl = row[18]
        casxy = Rxy
        maltxy = Lxy
    
    filestem = f"{VIDFOLDER}{fileprefix}{row[0]}"
    vidfile = f"{filestem}_{cam}.avi"
    posfile = f"{filestem}_{cam}{scorer}.h5"
    
    print(box, vidfile)
    
    output = make_ppp_video_dict(rat, diet, box, tank, cam,
                        vidfile, posfile, scorer,
                        cas_ttl, malt_ttl,
                        casxy, maltxy)
    
    return output
    

MASTERFOLDER = "/home/jaime/Github/PPP_video/"
VIDFOLDER = MASTERFOLDER + "videos/"
TANKFOLDER = MASTERFOLDER + "tanks/"
METAFILEFOLDER = MASTERFOLDER + "metafiles/"
SCORER = "DLC_resnet50_operantbox_LEMar29shuffle1_550000"
SCORER = "DLC_resnet50_PPP-DLCMay29shuffle1_600000"

PPP_video_data = []

# Analysis of PPP1
tp.metafilemaker(METAFILEFOLDER+"PPP1.xlsx", METAFILEFOLDER+"PPP1_metafile", fileformat='txt')
rows, header = tp.metafilereader(METAFILEFOLDER+"PPP1_metafile.txt")

rats = ['PPP1.1', 'PPP1.2', 'PPP1.3', 'PPP1.4', 'PPP1.5', 'PPP1.6', 'PPP1.7']
POS = [[(428, 115), (530, 114)], [(359, 92), (461, 97)]]

for row in rows:
    if row[2] in rats:
        if row[3] == 's10':
            out = process_row(row,
                              "PPP1-171017-081744_",
                              SCORER,
                              POS)
            PPP_video_data.append(out)

# Analysis of PPP3
tp.metafilemaker(METAFILEFOLDER+"PPP3.xlsx", METAFILEFOLDER+"PPP3_metafile",
                  sheetname="PPP3_metafile", fileformat='txt')
rows, header = tp.metafilereader(METAFILEFOLDER+"PPP3_metafile.txt")

rats = ['PPP3.2', 'PPP3.3', 'PPP3.4', 'PPP3.5','PPP3.8']
POS = [[(361, 51), (454, 52)], [(308, 294), (206, 306)]]

for row in rows:
    if row[2] in rats:
        if row[3] == 's10':
            out = process_row(row,
                              "PPP3-180705-135030_",
                              SCORER,
                              POS)
            PPP_video_data.append(out)

# Analysis of PPP4
tp.metafilemaker(METAFILEFOLDER+"PPP4.xlsx", METAFILEFOLDER+"PPP4_metafile",
                  sheetname="PPP4_metafile", fileformat='txt')
rows, header = tp.metafilereader(METAFILEFOLDER+"PPP4_metafile.txt")

rats = ['PPP4.1', 'PPP4.4']
POS = [[(275, 388), (175, 389)], [(309, 295), (207, 307)]]

for row in rows:
    if row[2] in rats:
        if row[3] == 's10':
            out = process_row(row,
                              "PPP4_1-191004-084000_",
                              SCORER,
                              POS)
            PPP_video_data.append(out)


pickle_out = MASTERFOLDER + "PPP_video_data.pickle"
with open(pickle_out, "wb") as dill_file:
    dill.dump(PPP_video_data, dill_file)