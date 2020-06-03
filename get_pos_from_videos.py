# -*- coding: utf-8 -*-
"""
Created on Tue May 26 11:58:18 2020

@author: jaime
"""

import deeplabcut


folder = "/home/jaime/Github/PPP_video/"
path_config_file = folder + "config.yaml"

videofile_path = [folder+"videos/"]

deeplabcut.analyze_videos(path_config_file, videofile_path, videotype='.avi')
