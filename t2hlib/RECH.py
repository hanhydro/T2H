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
    def __init__(self, rech, Rflag, top, nx, ny, lake_x, lake_y,\
                 CXrain, CYrain, lake_level, grid, xmax, xmin, ymax, ymin,\
                 dx, dy, Krain):
        if Rflag == 1:
            self.prcp = top
            # lake level correction since topo has no lake levels
            for i in range(len(lake_x)):
                self.ycoord = np.abs(lake_y[i] - int((grid-1)/2))
                self.xcoord = np.abs(lake_x[i] + int((grid-1)/2))
                self.prcp[self.ycoord, self.xcoord] = lake_level[i] # top + laketops
            for i in range(nx):
                for j in range(ny):
                    self.prcp[j, i] = max((rech+self.prcp[j,i]*Krain/1000), 0)
            if CXrain > 0:
                for i in range(nx):
                    for j in range(ny):
                        self.prcp[j, i] = max(self.prcp[j,i],\
                                   1+(xmin+i*dx-(xmax+xmin)/2)/CXrain)
            if CYrain > 0:
                for i in range(nx):
                    for j in range(ny):
                        self.prcp[j, i] = max(self.prcp[j,i],\
                                   1+(ymax-j*dy-(ymax+ymin)/2)/CYrain)
            print("RECHARGE console> precipitation input applied \n")
        elif Rflag == 2:
            print("precipitation option 2 is not available yet")
        elif Rflag == 3:
            print("precipitation option 3 is not available yet")

def precip(rech, Rflag, top, nx, ny, lake_x, lake_y,\
           CXrain, CYrain, lake_level, grid, xmax, xmin, ymax, ymin,\
           dx, dy, Krain):
    return precipitation(rech, Rflag, top, nx, ny, lake_x, lake_y,\
                         CXrain, CYrain, lake_level, grid, xmax,\
                         xmin, ymax, ymin, dx, dy, Krain)
# TISC function "max_water_in_air_column" is needed