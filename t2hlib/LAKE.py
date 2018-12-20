#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 18:55:22 2018
@author: kyungdoehan
"""

import pandas as pd
import os
import numpy as np
#%% Importing lake coordinates and levels
class lake_data_read:
    def __init__(self, Ma, nz, ny, nx, grid, rch, dx, dy, lakeK, top, bot, lakeb):
        self.btop = top - bot[0, :, :]
        self.dir = os.getcwd()
        self.df = pd.read_csv(self.dir + "/tisc_output/lakes_" + str(Ma) + "0Ma.csv",\
                 sep=",", header=None)
        self.df.columns = ["x", "y", "level", "type", "id", "nlakes", "evap", "inflow"]
        self.x = self.df.loc[:, "x"].values.astype('int')
        self.y = self.df.loc[:, "y"].values.astype('int')
        self.level = self.df.loc[:, "level"].values.astype('float32')
        self.type = self.df.loc[:, "type"].values.astype('str')
        self.id = self.df.loc[:, "id"].values.astype('int')
        self.nlakes = self.df.loc[:, "nlakes"].values.astype('int')[0]
        self.evap = self.df.loc[:, "evap"].values.astype('float32')
        self.inflow = self.df.loc[:, "inflow"].values.astype('float32')
        self.stages = [0 for i in range(self.nlakes)]
        self.levap = [0 for i in range(self.nlakes)]
        self.linflow = [0 for i in range(self.nlakes)]
        self.bdlknc = np.zeros((nz, ny, nx), dtype=np.float32)     
        for i in range(1, self.nlakes+1):
            for j in range(len(self.id)):
                if self.id[j] == i:
                    self.stages[i-1] = self.level[j]
                    self.levap[i-1] = self.evap[j]
                    self.linflow[i-1] = self.inflow[j]
        self.lakarr = np.zeros((nz, ny, nx), dtype=np.int16)
        self.ghb = [0 for i in range(len(self.id))]
        for lake in range(len(self.x)):
            self.ycoord = np.abs(self.y[lake] - int((grid-1)/2))
            self.xcoord = np.abs(self.x[lake] + int((grid-1)/2))
#            if self.ycoord != 0 and self.xcoord != 0:
#                if self.ycoord != ny-1 and self.xcoord != 0:
            self.lakarr[0, self.ycoord, self.xcoord] = self.level[lake]
            self.ghb[lake] = [0, self.ycoord, self.xcoord, self.level[lake],\
                             dx*dy*lakeK/(lakeb)]
        self.ghb = {0: self.ghb}
def lake(Ma, nz, ny, nx, grid, rch, dx, dy, lakeK, top, bot, lakeb):
    return lake_data_read(Ma, nz, ny, nx, grid, rch, dx, dy, lakeK, top, bot, lakeb)     
