#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 19:04:13 2018

@author: kyungdoehan
"""

import numpy as np
#%% Making square arrays of x, y, z of the overall topography
class XYZ_data:
    def __init__(self, a, x, y, z):
        self.X = np.zeros((a, a))
        self.Y = np.zeros((a, a))
        self.Z = np.zeros((a, a))
        for i in range(a):
            for j in range(a):
                self.X[j, i] = x[i + j * a]
                self.Y[j, i] = y[i + j * a]
                self.Z[j, i] = z[i + j * a]
def XYZ(grid, x, y, z):
    return XYZ_data(grid, x, y, z)
#%%
class delz_ratio:
    def __init__(self, i):
        self.dzratio = np.exp(np.arange(1, i + 1) / 10)
        self.dzratio = self.dzratio / np.sum(self.dzratio)
 
def dzratio(i):
    return delz_ratio(i)       
#%%     
class bottom:
    def __init__(self, inz, ifixed, j, top, dat_var, dat_new, dzratio):
        self.tot_b = top - dat_var
        self.bot = np.zeros((inz, j, j))
        for irow in range(j):
            for icol in range(j):
                self.bot[:, irow, icol] = top[irow, icol] - \
                np.cumsum(self.tot_b[irow, icol] * dzratio)
        self.bot_fixed = np.zeros((ifixed, j, j))
        self.bot_fixed[0, :, :] = self.bot[inz - 1, :, :] + dat_new / ifixed
        for i in range(ifixed - 1):
            self.bot_fixed[i+1, :, :] = self.bot_fixed[i, :, :]+dat_new/ifixed
        self.bot = np.vstack((self.bot, self.bot_fixed))
def bot(inz, ifixed, j, top, dat_var, dat_new, dzratio):
    return bottom(inz, ifixed, j, top, dat_var, dat_new, dzratio)
#%%
class delz:
    def __init__(self, top, bot, nz, ny, nx):
        self.dzs = np.zeros((nz, ny, nx), dtype=np.float32)
        self.dzs[0, :, :] = top - bot[0, :, :]
        for ilay in range(nz-1):
            self.dzs[ilay+1, :, :] = bot[ilay, :, :] - bot[ilay+1, :, :]
def dzs(top, bot, nz, ny, nx):
    return delz(top, bot, nz, ny, nx)
#%%
class nodes:
    def __init__(self, bot, dzs, nz, ny, nx):
        self.node = np.zeros((nz, ny, nx), dtype=np.float32)
        for irow in range(ny):
            for icol in range(nx):
                self.node[:, irow, icol] = bot[:, irow, icol] + 0.5 * dzs[:, irow, icol]
def node(bot, dzs, nz, ny, nx):
    return nodes(bot, dzs, nz, ny, nx)
    