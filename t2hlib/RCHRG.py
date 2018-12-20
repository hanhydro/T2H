#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 17:44:53 2018

@author: kyungdoehan
"""

import numpy as np
#%% Precipitation function from TISC 
# return array "precip"
class precipitation:
    def __init__(self, rech, Rflag, topo, nx, ny, lake_x, lake_y,\
                 CXrain, CYrain, lake_level, grid, xmax, xmin, ymax, ymin,\
                 dx, dy, Krain, ilak):
        if Rflag == 1:
            self.prcp_t = topo.copy()
            self.prcp = np.zeros((ny, nx), dtype=np.float32)
            self.rchrg = np.zeros((ny, nx), dtype=np.float32)
            self.ror = np.zeros((ny, nx), dtype=np.float32)
            # lake level correction since topo has no lake levels
            if ilak == 1:
                for i in range(len(lake_x)):
                    self.ycoord = np.abs(lake_y[i] - int((grid-1)/2))
                    self.xcoord = np.abs(lake_x[i] + int((grid-1)/2))
                    self.prcp_t[self.ycoord, self.xcoord] = lake_level[i] # top + laketops
                    topo[self.ycoord, self.xcoord] = lake_level[i]
            # Nichols (2000) Coefficents
            # 0.2032 m/yr < precip < 0.3048 m/yr = 0.008
            # 0.3048 m/yr < precip < 0.4604 m/yr = 0.130
            # 0.4604 m/yr < precip < 0.5080 m/yr = 0.144
            # 0.5080 m/yr < precip < 0.8636 m/yr = 0.158
            # precip > 0.8636 m/yr = 0.626
            
            for i in range(nx):
                for j in range(ny):
                    self.prcp_t[j, i] = max((rech + self.prcp_t[j, i] * Krain / 1E3), 0)
                    self.prcp[j, i] = self.prcp_t[j, i] - rech
                    if self.prcp[j, i] >= 0.8636:
                        self.prcp[j, i] = 0.626 * self.prcp[j, i]
                    elif self.prcp[j, i] >= 0.5080 and self.prcp[j, i] < 0.8636:
                        self.prcp[j, i] = 0.158 * self.prcp[j,i]
                    elif self.prcp[j, i] >= 0.4604 and self.prcp[j, i] < 0.5080:
                        self.prcp[j, i] = 0.144 * self.prcp[j, i]
                    elif self.prcp[j, i] >= 0.3048 and self.prcp[j, i] < 0.4604:
                        self.prcp[j, i] = 0.130 * self.prcp[j, i]
                    elif self.prcp[j, i] >= 0.2032 and self.prcp[j, i] < 0.3048:
                        self.prcp[j, i] = 0.008 * self.prcp[j, i]
                    else:
                        self.prcp[j, i] = 0
                        
            #lakes also will have recharge
#            for i in range(len(lake_x)):
#                self.ycoord = np.abs(lake_y[i] - int((grid-1)/2))
#                self.xcoord = np.abs(lake_x[i] + int((grid-1)/2))
#                self.ror[self.ycoord, self.xcoord] = 0 # laketop rech = 0
            for i in range(nx):
                for j in range(ny):
                    if j == 0:
                        self.rchrg[j, i] = 0
                    elif i == 0:
                        self.rchrg[j, i] = 0
                    elif j == ny-1:
                        self.rchrg[j, i] = 0
                    elif i == nx-1:
                        self.rchrg[j, i] = 0
                    else:
                        self.rchrg[j, i] = self.prcp[j,i]
            
            if CXrain > 0:
                for i in range(nx):
                    for j in range(ny):
                        self.prcp[j, i] = max(self.prcp[j,i], 1+(xmin+i*dx-(xmax+xmin)/2)/CXrain)
            if CYrain > 0:
                for i in range(nx):
                    for j in range(ny):
                        self.prcp[j, i] = max(self.prcp[j,i], 1+(ymax-j*dy-(ymax+ymin)/2)/CYrain)
            print("RECHARGE console> precipitation input applied \n")
        elif Rflag == 2:
            print("precipitation option 2 is not available yet")
        elif Rflag == 3:
            print("precipitation option 3 is not available yet")

def preci1(rech, Rflag, topo, nx, ny, lake_x, lake_y,\
           CXrain, CYrain, lake_level, grid, xmax, xmin, ymax, ymin,\
           dx, dy, Krain, ilak):
    return precipitation(rech, Rflag, topo, nx, ny, lake_x, lake_y,\
                         CXrain, CYrain, lake_level, grid, xmax,\
                         xmin, ymax, ymin, dx, dy, Krain, ilak)
# TISC function "max_water_in_air_column" is needed