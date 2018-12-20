#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PROPERTIES module for T2H beta
"""
import numpy as np
#%% Default hydraulic conductivity array for all cells
class conductivity:
    def __init__(self, nz, ny, nx, Kconst, hratio, node, top, lamb, rho_rock1):
            self.hk = np.ones((nz, ny, nx), dtype=np.float32)
            self.hk = self.hk * Kconst
            self.h_ini = hratio * node
            self.dens = np.ones((nz, ny, nx), dtype=np.float32)
            self.dens = self.dens * rho_rock1
#            for k in range(nz-1):
#                self.hk[k+1, :, :] = self.hk[k+1] * 10**(-lamb*(top - node[k+1, :, :]))
            
def cond(nz, ny, nx, Kconst, hratio, node, top, lamb, rho_rock1):
    return conductivity(nz, ny, nx, Kconst, hratio, node, top, lamb, rho_rock1)
#%% Sediment hydrologic properties for 1-sed models
class sediment:
    def __init__(self, dzs_cumsum, sed_b, sed_y, sed_x,\
                 grid, hk, nz, perm_sed, imodel, node, top, ny, nx, Kconst, lamb, rho_rock1, dens, rho_sed, bot, dy, dx, dzs):
        if imodel == 1 or imodel == 3:
            self.sed_len = len(sed_x)
            for cell in range(self.sed_len):
                for ilay in range(nz):
                    self.ycoord = np.abs(sed_y[cell] - int((grid - 1) / 2))
                    self.xcoord = np.abs(sed_x[cell] + int((grid - 1) / 2))
                    if sed_b[cell] >= dzs_cumsum[ilay, self.ycoord, self.xcoord]:
                        hk[ilay, self.ycoord, self.xcoord] = perm_sed
                        dens[ilay, self.ycoord, self.xcoord] = rho_sed
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        self.dens1 = 0
                        if dens[k, j, i] == rho_rock1:
                            if k != 1:
                                for l in range(k-1):
                                    self.dens1 = self.dens1 + dens[l, j, i] * dzs[l, j, i]
                            elif k == 1:
                                self.dens1 = dens[0, j, i]
                            self.dens1 = (self.dens1 + rho_rock1 * 0.5 * dzs[k, j, i])/(top[j, i] - node[k,j,i])
                            #print(self.dens1)
                            hk[k, j, i] = Kconst * np.exp(-lamb * self.dens1 * (top[j, i] - node[k,j,i]))
                        elif dens[k, j, i] == rho_sed:
                            if k != 1:
                                for l in range(k-1):
                                    self.dens1 = self.dens1 + dens[l, j, i] * dzs[l, j, i]
                            elif k == 1:
                                self.dens1 = dens[0, j, i]
                            self.dens1 = (self.dens1 + rho_sed * 0.5 * dzs[k, j, i])/(top[j, i] - node[k,j,i])    
                            hk[k, j, i] = perm_sed * np.exp(-lamb * self.dens1 * (top[j, i] - node[k,j,i]))
        else:
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        hk[k, j, i] = Kconst * np.exp(-lamb*rho_rock1*(top[j, i] - node[k, j, i]))
                        
def sed(dzs_cumsum, sed_b, sed_y, sed_x, grid, hk, nz, perm_sed, imodel, node, top, ny, nx, Kconst, lamb, rho_rock1, dens, rho_sed, bot, dy, dx, dzs):
    return sediment(dzs_cumsum, sed_b, sed_y, sed_x, grid, hk, nz, perm_sed, imodel, node, top, ny, nx, Kconst, lamb, rho_rock1, dens, rho_sed, bot, dy, dx, dzs)
#%% Vertical and horizontal flow barriers using fault geometry
class faults:
    def __init__(self, hydchr, grid, nz, ny, nx, fault_z, \
                 fault_y, fault_x, fault_x1, fault_y1, fault_z1, node, dat):
        self.hfb_pair = []
        self.hfb_info = np.zeros((nz, ny, nx), dtype=np.float32)
        self.ihfbs = -1
        for fault in range(len(fault_x)-1):
            self.ihfbs = self.ihfbs + 1
            if fault_y[self.ihfbs] == fault_y[self.ihfbs + 1]:
                # row by row operation
                if fault_z[self.ihfbs] <= dat or fault_z[self.ihfbs+1] <= dat:
                    print("cannot implement fault depth\
                          exceeding model depth \n")
                else:
                    # finding nearest node in depth-wise
                    self.ycoord = np.abs(fault_y[self.ihfbs] - int((grid-1)/2))
                    self.xcoord = np.abs(fault_x[self.ihfbs] + int((grid-1)/2))
                    self.ycoord1 = np.abs(fault_y[self.ihfbs+1]\
                                          - int((grid-1)/2))
                    self.xcoord1 = np.abs(fault_x[self.ihfbs+1]\
                                          + int((grid-1)/2))
                    self.val_1 = np.abs(np.abs(fault_z[self.ihfbs])\
                                        - np.abs(node[:, self.ycoord,\
                                                      self.xcoord]))
                    self.val_2 = np.abs(np.abs(fault_z[self.ihfbs+1])\
                                        - np.abs(node[:, self.ycoord1,\
                                                      self.xcoord1]))
                    self.ilay1 = np.where(self.val_1 == self.val_1.min())
                    self.ilay2 = np.where(self.val_2 == self.val_2.min())
                    self.ilay1 = np.max(self.ilay1)
                    self.ilay2 = np.max(self.ilay2)
                    # determine left or right
                    if self.ilay1 == self.ilay2:
                        if node[self.ilay1, self.ycoord, self.xcoord]\
                        < fault_z[self.ihfbs]:
                        # on the left (NW to SE fault)
                            if node[self.ilay2, self.ycoord1, self.xcoord1]\
                            > fault_z[self.ihfbs + 1]:
                                self.hfb_pair.append([self.ilay1, self.ycoord,\
                                                      self.xcoord,\
                                                      self.ycoord,\
                                                      self.xcoord1,\
                                                      hydchr])
                                self.hfb_info[self.ilay1, self.ycoord,\
                                              self.xcoord] = 1
                                self.hfb_info[self.ilay1,\
                                              self.ycoord,\
                                              self.xcoord1] = 1
                                                            
                        elif node[self.ilay1, self.ycoord, self.xcoord]\
                        > fault_z[self.ihfbs]:
                # on the right (NE to SW fault)
                            if node[self.ilay2, self.ycoord1, self.xcoord1]\
                            < fault_z[self.ihfbs+1]:
                                self.hfb_pair.append([self.ilay1, self.ycoord,\
                                                      self.xcoord, self.ycoord,\
                                                      self.xcoord1, hydchr])
                                self.hfb_info[self.ilay1, self.ycoord,\
                                              self.xcoord] = 1
                                self.hfb_info[self.ilay1, self.ycoord,\
                                              self.xcoord1] = 1
                    elif self.ilay1 != self.ilay2:
                        # Case I: fault above left node and plunging right (x+)
                        if  node[self.ilay1, self.ycoord, self.xcoord]\
                        < fault_z[self.ihfbs]:
                            if node[self.ilay2, self.ycoord1, self.xcoord1]\
                            < fault_z[self.ihfbs]:
                                for i in range(self.ilay2-self.ilay1):
                                    self.hfb_pair.append([self.ilay1+i-1, self.ycoord,\
                                                      self.xcoord, self.ycoord,\
                                                      self.xcoord1, hydchr])
                                    self.hfb_info[self.ilay1+i-1, self.ycoord,\
                                              self.xcoord] = 1
                                    self.hfb_info[self.ilay1+i-1, self.ycoord,\
                                              self.xcoord1] = 1  
                        # Case II: fault below left node                                                  
                        elif node[self.ilay1, self.ycoord, self.xcoord]\
                        >= fault_z[self.ihfbs]:                                
                            if node[self.ilay2, self.ycoord1, self.xcoord1]\
                            < fault_z[self.ihfbs]:
                                for i in range (self.ilay2-self.ilay1):
                                    self.hfb_pair.append([self.ilay1+i+1, self.ycoord,\
                                                      self.xcoord, self.ycoord,\
                                                      self.xcoord1, hydchr])
                                    self.hfb_info[self.ilay1+i+1, self.ycoord,\
                                              self.xcoord] = 1
                                    self.hfb_info[self.ilay1+i+1, self.ycoord,\
                                              self.xcoord1] = 1
                        # Case III: fault abobe left node and plunging left (x-)                                                 
                        elif node[self.ilay1, self.ycoord, self.xcoord]\
                        < fault_z[self.ihfbs] and self.ilay1 > self.ilay2:                                
                            if node[self.ilay2, self.ycoord1, self.xcoord1]\
                            < fault_z[self.ihfbs]:
                                for i in range (self.ilay1-self.ilay2):
                                    self.hfb_pair.append([self.ilay1-i, self.ycoord,\
                                                      self.xcoord, self.ycoord,\
                                                      self.xcoord1, hydchr])
                                    self.hfb_info[self.ilay1-i, self.ycoord,\
                                              self.xcoord] = 1
                                    self.hfb_info[self.ilay1-i, self.ycoord,\
                                              self.xcoord1] = 1  
                        # Case IV: fault below left node and plunging left (x-)                                                 
                        elif node[self.ilay1, self.ycoord, self.xcoord]\
                        >= fault_z[self.ihfbs] and self.ilay1 > self.ilay2:                                
                            if node[self.ilay2, self.ycoord1, self.xcoord1]\
                            < fault_z[self.ihfbs]:
                                for i in range (self.ilay1-self.ilay2):
                                    self.hfb_pair.append([self.ilay1-i, self.ycoord,\
                                                      self.xcoord, self.ycoord,\
                                                      self.xcoord1, hydchr])
                                    self.hfb_info[self.ilay1-i, self.ycoord,\
                                              self.xcoord] = 1
                                    self.hfb_info[self.ilay1-i, self.ycoord,\
                                              self.xcoord1] = 1 
                                                
            # y-propagating, horizontal
            if fault_x1[self.ihfbs] == fault_x1[self.ihfbs + 1]:
                # col by col operation
                if fault_z1[self.ihfbs] <= dat or\
                fault_z1[self.ihfbs+1] <= dat:
                    print("cannot implement fault depth exceeding model depth \n")
                else:
                    # finding nearest node in depth-wise
                    self.ycoord = np.abs(fault_y1[self.ihfbs]\
                                         - int((grid-1)/2))
                    self.xcoord = np.abs(fault_x1[self.ihfbs]\
                                         + int((grid-1)/2))
                    self.ycoord1 = np.abs(fault_y1[self.ihfbs+1]\
                                          - int((grid-1)/2))
                    self.xcoord1 = np.abs(fault_x1[self.ihfbs+1]\
                                          + int((grid-1)/2))
                    self.val_1 = np.abs(np.abs(fault_z1[self.ihfbs])\
                                   - np.abs(node[:, self.ycoord, self.xcoord]))
                    self.val_2 = np.abs(np.abs(fault_z1[self.ihfbs+1])\
                                   - np.abs(node[:, self.ycoord1,\
                                                 self.xcoord1]))
                    self.ilay1 = np.where(self.val_1 == self.val_1.min())
                    self.ilay2 = np.where(self.val_2 == self.val_2.min())
                    self.ilay1 = np.max(self.ilay1)
                    self.ilay2 = np.max(self.ilay2)
                        # determine left or right
                    if self.ilay1 == self.ilay2:
                        if node[self.ilay1, self.ycoord, self.xcoord]\
                        < fault_z1[self.ihfbs]:
                            # on the left (NW to SE fault)
                            if node[self.ilay2, self.ycoord1, self.xcoord1]\
                            > fault_z1[self.ihfbs+1]:
                                self.hfb_pair.append([self.ilay2, self.ycoord,\
                                                      self.xcoord, self.ycoord1,\
                                                      self.xcoord, hydchr])
                                self.hfb_info[self.ilay2, self.ycoord,\
                                              self.xcoord] = 1
                                self.hfb_info[self.ilay2, self.ycoord1,\
                                              self.xcoord] = 1
                        elif node[self.ilay1, self.ycoord, self.xcoord]\
                        > fault_z1[self.ihfbs]:
                        # on the right (NE to SW fault)
                            if node[self.ilay2, self.ycoord1, self.xcoord1]\
                            < fault_z1[self.ihfbs+1]:
                                self.hfb_pair.append([self.ilay2, self.ycoord,\
                                                      self.xcoord, self.ycoord1,\
                                                      self.xcoord, hydchr])
                                self.hfb_info[self.ilay2, self.ycoord,\
                                              self.xcoord] = 1
                                self.hfb_info[self.ilay2, self.ycoord1,\
                                              self.xcoord] = 1
                    elif self.ilay1 != self.ilay2:  
                        # Case I: fault above left node and plunging right (x+)
                        if  node[self.ilay1, self.ycoord, self.xcoord]\
                        < fault_z[self.ihfbs]:
                            if node[self.ilay2, self.ycoord1, self.xcoord1]\
                            < fault_z[self.ihfbs]:
                                for i in range(self.ilay2-self.ilay1):
                                    self.hfb_pair.append([self.ilay1+i-1, self.ycoord,\
                                                      self.xcoord, self.ycoord1,\
                                                      self.xcoord, hydchr])
                                    self.hfb_info[self.ilay1+i-1, self.ycoord,\
                                              self.xcoord] = 1
                                    self.hfb_info[self.ilay1+i-1, self.ycoord1,\
                                              self.xcoord] = 1  
                        # Case II: fault below left node                                                  
                        elif node[self.ilay1, self.ycoord, self.xcoord]\
                        >= fault_z[self.ihfbs]:                                
                            if node[self.ilay2, self.ycoord1, self.xcoord1]\
                            < fault_z[self.ihfbs]:
                                for i in range (self.ilay2-self.ilay1):
                                    self.hfb_pair.append([self.ilay1+i+1, self.ycoord,\
                                                      self.xcoord, self.ycoord1,\
                                                      self.xcoord, hydchr])
                                    self.hfb_info[self.ilay1+i+1, self.ycoord,\
                                              self.xcoord] = 1
                                    self.hfb_info[self.ilay1+i+1, self.ycoord1,\
                                              self.xcoord] = 1
                        # Case III: fault abobe left node and plunging left (x-)                                                 
                        elif node[self.ilay1, self.ycoord, self.xcoord]\
                        < fault_z[self.ihfbs] and self.ilay1 > self.ilay2:                                
                            if node[self.ilay2, self.ycoord1, self.xcoord1]\
                            < fault_z[self.ihfbs]:
                                for i in range (self.ilay1-self.ilay2):
                                    self.hfb_pair.append([self.ilay1-i, self.ycoord,\
                                                      self.xcoord, self.ycoord1,\
                                                      self.xcoord, hydchr])
                                    self.hfb_info[self.ilay1-i, self.ycoord,\
                                              self.xcoord] = 1
                                    self.hfb_info[self.ilay1-i, self.ycoord1,\
                                              self.xcoord] = 1  
                        # Case IV: fault below left node and plunging left (x-)                                                 
                        elif node[self.ilay1, self.ycoord, self.xcoord]\
                        >= fault_z[self.ihfbs] and self.ilay1 > self.ilay2:                                
                            if node[self.ilay2, self.ycoord1, self.xcoord1]\
                            < fault_z[self.ihfbs]:
                                for i in range (self.ilay1-self.ilay2):
                                    self.hfb_pair.append([self.ilay1-i, self.ycoord,\
                                                      self.xcoord, self.ycoord1,\
                                                      self.xcoord, hydchr])
                                    self.hfb_info[self.ilay1-i, self.ycoord,\
                                              self.xcoord] = 1
                                    self.hfb_info[self.ilay1-i, self.ycoord1,\
                                              self.xcoord] = 1                                             
                                          
def hfb(hydchr, grid, nz, ny, nx, fault_z,\
        fault_y, fault_x, fault_x1, fault_y1, fault_z1, node, dat):
    return faults(hydchr, grid, nz, ny, nx, fault_z, \
                 fault_y, fault_x, fault_x1, fault_y1, fault_z1, node, dat)
# Vertical flow barriers
class vertical:
    def __init__(self, nz, ny, nx, fault_x, fault_y, fault_z,\
                 node, grid, hydchr, dzratio):
        self.vkcb_data = np.ones((nz, ny, nx), dtype=np.float32)
        self.vkcb = np.ones((len(fault_x)), dtype=np.float32)
        self.layer = np.zeros((ny, nx), dtype=np.float32)
        for faults in range(len(fault_x)):
            self.ycoord = np.abs(fault_y[faults] - int((grid-1)/2))
            self.xcoord = np.abs(fault_x[faults] + int((grid-1)/2))
            self.val_1 = np.abs(fault_z[faults]\
                                - node[:, self.ycoord, self.xcoord])
            self.ilay1 = np.where(self.val_1 == self.val_1.min())
            self.ilay1 = np.max(self.ilay1)
            # node is above the fault
            if node[self.ilay1, self.ycoord, self.xcoord] >= fault_z[faults]:
                self.vkcb_data[self.ilay1, self.ycoord, self.xcoord] = hydchr
                self.vkcb[faults] = 1
                self.layer[self.ycoord, self.xcoord]\
            = node[self.ilay1, self.ycoord, self.xcoord]
            elif node[self.ilay1, self.ycoord, self.xcoord]\
            < fault_z[faults]: # node is below the fault
                if self.ilay1 > 0:
                    self.vkcb_data[self.ilay1-1, self.ycoord,\
                                   self.xcoord] = hydchr
                    self.vkcb[faults] = 1
                    self.layer[self.ycoord, self.xcoord]\
            = node[self.ilay1, self.ycoord, self.xcoord]
        
def vfb(nz, ny, nx, fault_x, fault_y, fault_z, node, grid, hydchr, dzratio):
    return vertical(nz, ny, nx, fault_x, fault_y, fault_z,\
                    node, grid, hydchr, dzratio)
    
class layer_type:
    def __init__(self, switch, nz, inz):
        self.laytyp = np.ones([nz])
        if switch == 1:
            for i in range(nz):
                self.laytyp[i] = 1
                if i > inz-1:
                    self.laytyp[i] = 0
                
def laytyp(switch, nz, inz):
    return layer_type(switch, nz, inz)
#%% Lake treatment added for the DRN package    
class drn:
    def __init__(self, top, bot, dy, dx, ny, nx, hk, x, y, ibound, grid):
        self.btop = top - bot[0, :, :]
        self.conductance = []
        self.con = dx * dy * hk[0, :, :] / (self.btop / 2) # only surface grid-blocks
        self.df = []
        self.condarr = np.zeros((ny, nx), dtype=np.float32)
        for j in range(grid):
            for i in range(grid):
                if j != 0 and i != 0 and j != ny-1 and i != nx-1:
                    self.condarr[j, i] = self.con[j,i]
                    self.conductance.append([0, j, i, top[j, i], self.con[j, i]])
            
#        for lake in range(len(lake_x)):
#            self.lx[lake] = int(np.abs(lake_x[lake] + int((grid-1)/2)))
#            self.ly[lake] = int(np.abs(lake_y[lake] - int((grid-1)/2)))
#        for j in range(grid):
#            for i in range(grid):
#                if j == self.ly[j] and i == self.lx[i]:
#                    self.condarr[j, i] = 0
#                elif j == 0 or i == 0 or j == ny-1 or i == nx-1:
#                    self.condarr[j, i] = 0
#                else:
#                    self.condarr[j, i] = self.con[j, i]
#                    self.conductance.append([0, j, i, top[j, i], self.con[j, i]])
        
def cd(top, bot, dy, dx, ny, nx, hk, x, y, ibound, grid):
    return drn(top, bot, dy, dx, ny, nx, hk, x, y, ibound, grid)
#%% ibound < 0: const. head; = =0: inactive; >0: active
class BC:
    def __init__(self, bot, nz, ny, nx, grid):
        self.ibound = np.ones((nz, ny, nx), dtype=np.int)
        self.ibound[:, :, 0] = 0
        self.ibound[:, 0, :] = 0
        self.ibound[:, ny - 1, :] = 0
        self.ibound[:, :, nx - 1] = 0
        self.ibound[nz-1, :, :] = 0     
        
def boundary(bot, nz, ny, nx, grid):
    return BC(bot, nz, ny, nx, grid)

class precipitation:
    def __init__(self, top):
        self.mintop = np.min(top)
        self.maxtop = np.max(top)
        
def frech(top):
    return precipitation(top)