#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:43:11 2018
@author: kyungdoehan

"""
#%% Required packages
import os
import numpy as np
#import multiprocessing as mp
from t2hlib import READDATA, LAKE, PACKAGES, PROPERTIES, GRID, PLOT, RCHRG
from os.path import expanduser
#%% Setting up directories
home = expanduser('~')
flopypth = os.path.join(home+'/anaconda3/lib/python3.6/site-packages')
#%% Initial variables
Ma = 12.5
nz = 40
nz_fixed = 10
inz = 30
dx = 1000
dy = 1000
dat = -10000
dat_var = -3000
idat = dat - dat_var
rech = 0.4 # background runoff
target = 75
Kconst = 100
perm_sed = Kconst * 2
hratio = 1.0
hydchr = Kconst / 10 # 100 m, thickness
iskip = 4
ivtk = 2
h_tol = 1E-4
Rflag = 1
CXrain = 0
CYrain = 0
Krain = 1
lakeK = 100
lakeb = 10.0
lamb = 1E-8
# sediment and rock density (arbitrary)
rho_sed = 1700
rho_rock1 = 3000
#%% Reading data file by specified Mas
class __main__:
    def __init__(self, Ma, nz, nz_fixed, inz, dx, dy, dat, dat_var, idat\
                 , rech, perm_sed, target_row, Kconst, hratio, hydchr,\
                 iskip, ivtk, h_tol, Rflag, CXrain, CYrain, Krain,\
                 lakeK, lakeb, lamb, rho_rock1, rho_sed):
        """ data checker """
        datdef = READDATA.data(Ma)
        self.topoT = datdef.topo
        self.bdtopoT = datdef.bdtopo
        self.sedT = datdef.sed
        self.lakeT = datdef.lake
        self.faultT = datdef.fault
        print("console> data checker complete")
        """ Reading topography """
        topo = READDATA.topo_read(Ma)
        self.x = topo.x
        self.y = topo.y
        self.z = topo.z
        # Discretization variables (hydro. model specified)
        grid = int(np.sqrt(len(self.x)))
        self.nx = grid
        self.ny = grid
        griddata = GRID.XYZ(grid, self.x, self.y, self.z)
        self.X = griddata.X
        self.Y = griddata.Y
        self.top = griddata.Z
        self.top_l = self.top.copy()
        self.top2 = self.top.copy()
        self.top3 = self.top.copy()
        """ Check optional files whether they are there or not there """
        if self.sedT == True:
            print("console> reading sed data")
            sed = READDATA.sed_read(Ma)
            self.sed_x = sed.x
            self.sed_y = sed.y
            self.sed_b = sed.z
        elif self.sedT == False:
            print("console> sed data is missing for this period")
            self.sed_x = []
            self.sed_y = []
            self.sed_b = []
        if self.faultT == True:
            print("console> reading fault data")
            fault = READDATA.fault_read(Ma)
            self.fault_x = fault.x
            self.fault_y = fault.y
            self.fault_z = fault.z
            fault_rearr = READDATA.fault_rearrange(Ma)
            self.fault_x1 = fault_rearr.x
            self.fault_y1 = fault_rearr.y
            self.fault_z1 = fault_rearr.z
        elif self.faultT == False:
            print("console> fault data is missing for this period")
        if self.faultT == True and self.sedT == True:
            self.imodel = 1
        elif self.faultT == True and self.sedT == False:
            self.imodel = 2
        elif self.faultT == False and self.sedT == True:
            self.imodel = 3
        elif self.faultT == False and self.sedT == False:
            self.imodel = 4
        # Exponentially increasing dz values
        dzr = GRID.dzratio(inz)
        self.dzratio = dzr.dzratio
        bot = GRID.bot(inz, nz_fixed, grid, self.top, dat_var, idat, self.dzratio).bot
        self.bot_l = bot.copy()
        self.bot1 = bot.copy()
        self.bot2 = bot.copy()
        dzs = GRID.dzs(self.top, bot, nz, self.ny, self.nx).dzs
        self.dzs_cumsum = np.cumsum(dzs, axis = 0)
        node = GRID.node(bot, dzs, nz, self.ny, self.nx).node
        self.nod = node.copy()
        # * .9 # Initial value of water table according to the topo.


        if self.lakeT == True:
            lake = LAKE.lake(Ma, nz, self.ny, self.nx, grid, rech,\
                             dx, dy, lakeK, self.top_l, self.bot_l, lakeb)
            self.lake_x = lake.x
            self.lake_y = lake.y
            self.lakarr = lake.lakarr
            self.lake_level = lake.level
            self.ilak = 1
            self.nlakes = lake.nlakes
            self.ghb = lake.ghb
            #print(self.ghb)
        elif self.lakeT == False:
            self.ilak = 0
            print("console> lake data is missing for this period \n")
        conductivity = PROPERTIES.cond(nz, self.ny, self.nx, Kconst, hratio,\
                                       node, self.top, lamb, rho_rock1)    
        self.h_ini = conductivity.h_ini
        self.hk = conductivity.hk
        self.dens = conductivity.dens
        self.dens1 = self.dens.copy()

        self.sediment = PROPERTIES.sed(self.dzs_cumsum, self.sed_b,\
                                           self.sed_y, self.sed_x, grid,\
                                           self.hk, nz, perm_sed, self.imodel,\
                                           self.nod, self.top, self.ny, self.nx,\
                                           Kconst, lamb, rho_rock1,\
                                           self.dens1, rho_sed, self.bot1, dy, dx, dzs)
        #print(self.hk[:, 85, 105])
        if self.faultT == True:
            fault = PROPERTIES.hfb(hydchr, grid, nz, self.ny, self.nx,\
                                   self.fault_z,\
                                   self.fault_y, self.fault_x,\
                                   self.fault_x1, self.fault_y1,\
                                   self.fault_z1, node, dat)
            self.hfb_pair = fault.hfb_pair
            self.hfb_info = fault.hfb_info
            vfb = PROPERTIES.vfb(nz, self.ny, self.nx, self.fault_x1,\
                                 self.fault_y1, self.fault_z1,\
                                 node, grid, hydchr, self.dzratio)
            self.vkcb_data = vfb.vkcb_data
            self.vkcb_layers = vfb.layer
        elif self.faultT == False:
            self.hfb_pair = []
            self.vkcb_data = []
            self.hfb_info = []
            
        switch = 2
        
        laytyp = PROPERTIES.laytyp(switch, nz, inz).laytyp

        self.xx = (0.5*(grid-1))*dx
        self.xn = -(0.5*(grid-1))*dx
        self.yx = (0.5*(grid-1))*dy
        self.yn = -(0.5*(grid-1))*dy 
        self.topo = self.top.copy()
        if self.lakeT == True:
            prec = RCHRG.preci1(rech, Rflag, self.topo, self.nx, self.ny,\
                                        self.lake_x, self.lake_y, CXrain, CYrain,\
                                        self.lake_level, grid, self.xx,\
                                        self.xn, self.yx, self.yn, dx, dy,\
                                        Krain, self.ilak)
            self.rchrg = prec.rchrg
            self.prcp_t = prec.prcp_t
        elif self.lakeT == False:
            self.lake_x = []
            self.lake_y = []
            self.lake_level = []
            prec = RCHRG.preci1(rech, Rflag, self.topo, self.nx, self.ny,\
                                        self.lake_x, self.lake_y, CXrain, CYrain,\
                                        self.lake_level, grid, self.xx,\
                                        self.xn, self.yx, self.yn, dx, dy,\
                                        Krain, self.ilak)
            self.rchrg = prec.rchrg
            self.prcp_t = prec.prcp_t
            self.nlakes = 0
            self.ghb = []

        self.ibound  = PROPERTIES.boundary(bot, nz, self.ny, self.nx, grid).ibound
                                 
        self.conductance = PROPERTIES.cd(self.top3,\
                                         self.bot2, dy, dx,\
                                         self.ny, self.nx,\
                                         self.hk, self.x, self.y,\
                                         self.ibound, grid).conductance  
        self.condarr = PROPERTIES.cd(self.top3,\
                                         self.bot2, dy, dx,\
                                         self.ny, self.nx,\
                                         self.hk, self.x, self.y,\
                                         self.ibound, grid).condarr
        # Model paths from external subroutine
        self.path_info = PACKAGES.paths()
        self.fdirmodel = self.path_info.fdirmodel
        self.fnmmodel = self.path_info.fnmmodel
        # Model packages from external subroutine
        modelpackage = PACKAGES.packages(nz, self.ny, self.nx, dy, dx,\
                                         self.top, bot, self.ibound,\
                                         self.h_ini, self.hfb_pair,\
                                         self.conductance, laytyp, self.hk,\
                                         self.vkcb_data, self.rchrg, self.fdirmodel,\
                                         self.fnmmodel, self.nlakes, self.ilak\
                                         , self.ghb, self.imodel)

            
        self.mf = modelpackage.mf
        self.mfdis = modelpackage.mfdis
        self.mfbas = modelpackage.mfbas
        if self.faultT == True:
            self.mfhfb = modelpackage.mfhfb
        self.mfdrn = modelpackage.mfdrn
        self.mfupw = modelpackage.mfupw
        self.mfrch = modelpackage.mfrch
        self.mfnwt = modelpackage.mfnwt
        if self.lakeT == True:
            self.mfghb = modelpackage.mfghb
        self.output = modelpackage.output
def main(Ma, nz, nz_fixed, inz, dx, dy, dat, dat_var, idat\
         , rech, perm_sed, target_row, Kconst, hratio, hydchr,\
         iskip, ivtk, h_tol, Rflag, CXrain, CYrain, Krain,\
         lakeK, lakeb, lamb, rho_rock1, rho_sed):
    return __main__(Ma, nz, nz_fixed, inz, dx, dy, dat, dat_var, idat\
                    , rech, perm_sed, target_row, Kconst, hratio, hydchr,\
                    iskip, ivtk, h_tol, Rflag, CXrain, CYrain, Krain,\
                    lakeK, lakeb, lamb, rho_rock1, rho_sed)
#%% Model execution    
model = main(Ma, nz, nz_fixed, inz, dx, dy, dat, dat_var, idat\
             , rech, perm_sed, target, Kconst, hratio, hydchr,\
             iskip, ivtk, h_tol, Rflag, CXrain, CYrain, Krain,\
             lakeK, lakeb, lamb, rho_rock1, rho_sed)
#%% Model input checker
mf = model.mf
mf.dis.check()
mf.drn.check()
if model.lakeT == True:
    mf.ghb.check()
mf.rch.check()
mf.upw.check()
mf.write_input()
mf.run_model()
#%%
print("model recharge average: ", np.mean(model.rchrg))
print("model topography average: ", np.mean(model.top))
print("model N/K: ", np.mean(model.rchrg)/np.mean(model.hk))
#%%
# Plot the grid, head contour, and discharge vector on a cross-section
target = 75
import matplotlib.pyplot as plt
import flopy   
figheadxsect, axheadxsect = plt.subplots(figsize=(40,5))
mfxsect = PLOT.fmfxsect(mf, model.mfdis, target, axheadxsect).mfxsect
a = PLOT.head(mf, model.fdirmodel).a
headc = PLOT.headc(mfxsect, a)
if model.imodel > 4:
    p = PLOT.pat(mfxsect, model.hk, Kconst, model.vkcb_data,\
                 model.hfb_info, nz, model.ny, model.nx, model.imodel)
#p1 = PLOT.wtable(mfxsect, model.top, a, nz, model.ny, model.nx)
#patch1 = p1.plarr
#patch1 = p.bcarray
# Plotting the grid mesh
#gdplot = mfxsect.plot_grid(color='r', linewidths=0.2)
# Plotting the BC (blue = const. head; noflow = black)
BCplot = mfxsect.plot_ibound(model.ibound, color_noflow = 'black',\
                         color_ch = 'blue', head = a)
wtable = PLOT.wt_f(nz, model.ny, model.nx, a, model.top)
wtf = wtable.wtf
wtplot = mfxsect.plot_surface(wtf[0,:,:])
topplot = mfxsect.plot_surface(model.top)
cellbudget = PLOT.cbc(model.fdirmodel, mf)
times = cellbudget.times
# Plot flux vector
qx = cellbudget.qx # right face
qy = cellbudget.qy # front face
qz = cellbudget.qz # lower face
# Average flows to cell centers
avgq = PLOT.qavg(qz, qy, qx, nz, model.ny, model.nx)
qx_avg = avgq.qx_avg
qy_avg = avgq.qy_avg
qz_avg = avgq.qz_avg

y, x, z = model.mfdis.get_node_coordinates()
#diagonal
# [0,0] [1, 1], [2,2] cross-section
# x = 1000*sqrt(2)/2 increment
# z = 
zdiag = np.zeros((nz, model.nx), dtype=np.float32)
xdiag = np.zeros((nz, model.nx), dtype=np.float32)
hdiag = np.zeros((nz, model.nx), dtype=np.float32)
qxdiag = np.zeros((nz, model.nx), dtype=np.float32)
qydiag = np.zeros((nz, model.nx), dtype=np.float32)
qzdiag = np.zeros((nz, model.nx), dtype=np.float32)
tdiag = np.zeros((model.nx), dtype=np.float32)
for i in range(model.nx):
    # [40, 201]
    zdiag[:, i] = z[:, i, i]
    xdiag[:, i] = i*(dx*np.sqrt(2)/2)
    hdiag[:, i] = a[:, i, i]
    tdiag[i] = model.top[i, i]
    qxdiag[:, i] = qx_avg[:, i, i]
    qydiag[:, i] = qy_avg[:, i, i]
    qzdiag[:, i] = qz_avg[:, i, i]
if target >= 0:
    X, Z = np.meshgrid(x, z[:, target, 0])
    Z = z[:, target, :]
elif target < 0:
    X, Z = np.meshgrid(y, z[:, 0, np.abs(target)])
    Z = z[:, :, np.abs(target)]
elif target == "a":
    X = xdiag
    Z = zdiag

flowlines = PLOT.quiver(X, Z, qx_avg, qy_avg, qz_avg, target,\
                            iskip, axheadxsect).quiver
axheadxsect.set_xlabel('Model width (m)', fontsize = 12)
axheadxsect.set_ylabel('Model height (m)', fontsize = 12)
plt.show()
#%% Streamlines
import plotly.plotly as py
from t2hlib import PLOT
target=75
y, x, z = model.mfdis.get_node_coordinates()
if target >= 0:
    X, Z = np.meshgrid(x, z[:, target, 0])
#    Z = z[:, target, :]
elif target < 0:
    X, Z = np.meshgrid(y, z[:, 0, np.abs(target)])
#    Z = z[:, :, np.abs(target)]

s =  PLOT.sline( X, Z, qx_avg, qz_avg, target)
#%%
# Plot the grid, head contour, and discharge vector on a cross-section
target = 75
import matplotlib.pyplot as plt
import flopy   
figheadxsect, axheadxsect = plt.subplots(figsize=(40,5))
mfxsect = PLOT.fmfxsect(mf, model.mfdis, target, axheadxsect).mfxsect
a = PLOT.head(mf, model.fdirmodel).a
headc = PLOT.headc(mfxsect, a)
if model.imodel > 4:
    p = PLOT.pat(mfxsect, model.hk, Kconst, model.vkcb_data,\
                 model.hfb_info, nz, model.ny, model.nx, model.imodel)
#p1 = PLOT.wtable(mfxsect, model.top, a, nz, model.ny, model.nx)
#patch1 = p1.plarr
#patch1 = p.bcarray
# Plotting the grid mesh
#gdplot = mfxsect.plot_grid(color='r', linewidths=0.2)
# Plotting the BC (blue = const. head; noflow = black)
BCplot = mfxsect.plot_ibound(model.ibound, color_noflow = 'black',\
                         color_ch = 'blue', head = a)

wtable = PLOT.wt_f(nz, model.ny, model.nx, a, model.top)
wtf = wtable.wtf
wtf1 = model.top

wtplot = mfxsect.plot_surface(wtf[0,:,:])

topplot = mfxsect.plot_surface(model.top)

#%%
# HK grid plot
target = 75
import matplotlib.pyplot as plt
import flopy   
figheadxsect, axheadxsect = plt.subplots(figsize=(40,5))
mfxsect = PLOT.fmfxsect(mf, model.mfdis, target, axheadxsect).mfxsect
mfxsect.plot_array(model.hk)
#%%
cellbudget = PLOT.cbc(model.fdirmodel, mf)
times = cellbudget.times
# Plot flux vector
qx = cellbudget.qx # right face
qy = cellbudget.qy # front face
qz = cellbudget.qz # lower face
# Average flows to cell centers
avgq = PLOT.qavg(qz, qy, qx, nz, model.ny, model.nx)
qx_avg = avgq.qx_avg
qy_avg = avgq.qy_avg
qz_avg = avgq.qz_avg

y, x, z = model.mfdis.get_node_coordinates()
#diagonal
# [0,0] [1, 1], [2,2] cross-section
# x = 1000*sqrt(2)/2 increment
# z = 
zdiag = np.zeros((nz, model.nx), dtype=np.float32)
xdiag = np.zeros((nz, model.nx), dtype=np.float32)
hdiag = np.zeros((nz, model.nx), dtype=np.float32)
qxdiag = np.zeros((nz, model.nx), dtype=np.float32)
qydiag = np.zeros((nz, model.nx), dtype=np.float32)
qzdiag = np.zeros((nz, model.nx), dtype=np.float32)
tdiag = np.zeros((model.nx), dtype=np.float32)
for i in range(model.nx):
    # [40, 201]
    zdiag[:, i] = z[:, i, i]
    xdiag[:, i] = i*(dx*np.sqrt(2)/2)
    hdiag[:, i] = a[:, i, i]
    tdiag[i] = model.top[i, i]
    qxdiag[:, i] = qx_avg[:, i, i]
    qydiag[:, i] = qy_avg[:, i, i]
    qzdiag[:, i] = qz_avg[:, i, i]
if target >= 0:
    X, Z = np.meshgrid(x, z[:, target, 0])
    Z = z[:, target, :]
elif target < 0:
    X, Z = np.meshgrid(y, z[:, 0, np.abs(target)])
    Z = z[:, :, np.abs(target)]
elif target == "a":
    X = xdiag
    Z = zdiag

flowlines = PLOT.quiver(X, Z, qx_avg, qy_avg, qz_avg, target,\
                            iskip, axheadxsect).quiver
axheadxsect.set_xlabel('Model width (m)', fontsize = 12)
axheadxsect.set_ylabel('Model height (m)', fontsize = 12)
plt.show()
#%% model top
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from matplotlib import cbook
from matplotlib import cm
from matplotlib.colors import LightSource
import matplotlib.pyplot as plt
import numpy as np
y, x, z = model.mfdis.get_node_coordinates()
x, y = np.meshgrid(x, y)
z = model.top
fig, ax = plt.subplots(subplot_kw = dict(projection='3d'), figsize=(20, 10))
ls = LightSource(270, 45)
rgb = ls.shade(z, cmap=cm.gist_earth, vert_exag=10, blend_mode='soft')
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=rgb,
                       linewidth=0, antialiased=True, shade=False)

plt.show()
#%% Lake discharge
cellbudget = PLOT.cbc(model.fdirmodel, mf)
times = cellbudget.times
# Plot flux vector
qx = cellbudget.qx # right face
qy = cellbudget.qy # front face
qz = cellbudget.qz # lower face
# Average flows to cell centers
avgq = PLOT.qavg(qz, qy, qx, nz, model.ny, model.nx)
lakdis = np.zeros((model.ny, model.nx), dtype=np.float32)
for j in range (model.ny):
    for i in range (model.nx):
        if model.lakarr[0, j, i] > 0:
            if qz_avg[0, j, i] < 0:
                lakdis[j, i] = - qz_avg[0, j, i]
lakedis = np.sum(lakdis)
lakedis_avg = np.mean(lakdis)
totdis = 0
ni = 0
for j in range(model.ny):
    for i in range(model.nx):
        if -qz_avg[0, j, i] > 0:
            totdis = totdis - qz_avg[0, j, i]
            ni = ni +1
ddddd = np.sum(model.rchrg)
ddd = np.mean(model.rchrg)

#%% WTF
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from matplotlib import cbook
from matplotlib import cm
from matplotlib.colors import LightSource
import matplotlib.pyplot as plt
import numpy as np
y, x, z = model.mfdis.get_node_coordinates()
wtable = PLOT.wt_f(nz, model.ny, model.nx, a, model.top)
wtf = wtable.wtf
x, y = np.meshgrid(x, y)
z = wtf[0, :, :]
fig, ax = plt.subplots(subplot_kw = dict(projection='3d'), figsize=(20, 10))
ls = LightSource(270, 45)
rgb = ls.shade(z, cmap=cm.gist_earth, vert_exag=10, blend_mode='soft')
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=rgb,
                       linewidth=0, antialiased=False, shade=False)
# make the panes transparent
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
# make the grid lines transparent
ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.set_zlim(0, 3000)
m = cm.ScalarMappable(cmap=cm.gist_earth)
m.set_array(z)
fig.colorbar(m, shrink=0.5, aspect=5)
ax.set_title('Water table elevation %1.1f Ma' %(Ma))

plt.show()
#%% Recharge surface plot
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize = (12, 12))
xsurf = np.arange(dx/2, model.nx*dx, dx)
ysurf = np.arange(dx/2, model.ny*dx, dy)
Xs, Ys = np.meshgrid(xsurf, ysurf)
Zs = model.rchrg[1:model.ny-1, 1:model.nx-1]
im = ax.imshow(Zs, cmap = cm.RdBu)
cset = contour(Zs, np.arange(0, 2.0, 0.1))
clabel(cset, inline = True, fmt='%1.1f', fontsize = 10)
colorbar(im)
title('Recharge (m/yr) %1.1f Ma' %(Ma))
show()
#%% Lake checker
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
import matplotlib.pyplot as plt
if model.lakeT == True:
    fig, ax = plt.subplots(figsize = (12, 12))
    xsurf = np.arange(dx/2, model.nx*dx, dx)
    ysurf = np.arange(dx/2, model.ny*dx, dy)
    Xs, Ys = np.meshgrid(xsurf, ysurf)
    Zs = model.lakarr[0, 1:model.ny-1, 1:model.nx-1]
    im = ax.imshow(Zs, cmap = cm.RdBu)
    #cset = contour(Zs, np.arange(0, 1.5, 0.1))
    #clabel(cset, inline = True, fmt='%1.1f', fontsize = 10)
    colorbar(im)
    title('Lake stage (m) %1.1f Ma' %(Ma))
    show()
#%% Water table plot
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize = (15, 15))
xsurf = np.arange(dx/2, model.nx*dx, dx)
ysurf = np.arange(dx/2, model.ny*dx, dy)
Xs, Ys = np.meshgrid(xsurf, ysurf)
wtable = PLOT.wt_f(nz, model.ny, model.nx, a, model.top)
wtf = wtable.wtf + model.lakarr[0, :, :]
Zs = wtf[0, 1:model.ny-1, 1:model.nx-1]
im = ax.imshow(Zs, cmap = cm.RdBu, vmin = 0, vmax = 3E3)
cset = contour(Zs, np.arange(0, 3000, 5E2))
clabel(cset, inline = True, fmt='%1.1f', fontsize = 12)
colorbar(im)
title("Water table (m); %1.1f Ma" %(Ma))
show()
#%% DRN plot
from pylab import imshow
fig, ax = plt.subplots(figsize = (12, 12))
dny = model.ny
dnx = model.nx
imdrn = model.condarr[1:dny-1, 1:dnx-1]
cset = contour(imdrn, np.arange(np.amin(imdrn), np.amax(imdrn), 1E6))
clabel(cset, inline = True, fmt='%1.1f', fontsize = 12)
im = ax.imshow(imdrn)
colorbar(im)
title("DRN plot")
show()
#%%
dny = model.ny
dnx = model.nx
a11 = model.prcp_t * 1000 * 1000
a22 = np.zeros((dny, dnx), dtype=np.float32)
for j in range(dny):
    for i in range(dnx):
        if a[0, j, i] >= model.top[j, i]:
            if -qz_avg[0, j, i] > 0:
                a22[j, i] = -qz_avg[0, j, i]
a33 = np.sum(a22)/np.sum(a11)
#%% Discharge plot
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
fig, ax = plt.subplots(figsize = (12, 12))
dny = model.ny
dnx = model.nx
imdisc = np.zeros((dny, dnx), dtype=np.float32)
for j in range(dny):
    for i in range(dnx):
        if a[0, j, i] >= model.top[j, i]:
            if -qz_avg[0, j, i] > 0:
                imdisc[j, i] = -qz_avg[0, j, i]
im = ax.imshow(imdisc, cmap = cm.plasma, vmin = 0, vmax = 1E7)
#cset = contour(imdisc, np.arange(np.amin(imdisc), np.amax(imdisc), 1E5))
#clabel(cset, inline = True, fmt='%1.1f', fontsize = 10)
colorbar(im)
title("Discharge (flux; m/yr) %1.1f Ma" %(Ma))
#ax.set_aspect('auto')
show()
np.mean(imdisc)
np.sum(imdisc)
#%% Lake recharge flux to sediment check 
#from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
#fig, ax = plt.subplots(figsize = (12, 12))
#dny = model.ny
#dnx = model.nx
#imdisc = np.zeros((dny, dnx), dtype=np.float32)
#for j in range(dny):
#    for i in range(dnx):
#        if a[0, j, i] >= model.top[j, i]:
#            if qz_avg[0, j, i] > 0:
#                imdisc[j, i] = qz_avg[0, j, i]
#im = ax.imshow(imdisc, cmap = cm.RdBu)
##cset = contour(imdisc, np.arange(np.amin(imdisc), np.amax(imdisc), 1E5))
##clabel(cset, inline = True, fmt='%1.1f', fontsize = 10)
#colorbar(im)
#title("Lake (flux; m/yr) %1.1f Ma" %(Ma))
#show()
#%%
ax1 = plt.subplots(figsize=(12,12))
mp = flopy.plot.map.ModelMap(sr=None, ax=None, model=mf,\
                             layer=0, dis=model.mfdis, extent=None, rotation=0.0)
mpp = mp.plot_discharge(frf=qx, fff=qy, flf=-qz, head=a,\
                        istep = 5, jstep = 5, normalize=True, color='b')
plt.show()        
# water table fluctuation
wtable = PLOT.wt_f(nz, model.ny, model.nx, a, model.top)
wtf = wtable.wtf
#%%
dis_loc = PLOT.dischloc(a, model.top, model.ny, model.nx, qz_avg)
dis_loc = dis_loc.discharge
wtcontr = PLOT.wtcont(x, y, dx, dy, model.ny, model.nx, wtf)
wtcontr = wtcontr.plot
#%% Output to VTK ()
ivtk = 1
lakarr = model.lakarr
if ivtk == 1:
    PLOT.vtk_export(mf, a, qx_avg, qy_avg, qz_avg,\
                        model.hk, model.vkcb_data, model.hfb_info, wtf, lakarr)
print("console> VTK output completed \n")
#%%
plt.imshow(np.transpose(model.hk[0, :, :]), aspect='auto')
plt.colorbar()
plt.show()
#%%
figg = plt.subplots(figsize=(40,5))
figg = plt.imshow(model.hk[:, 101, :])
plt.colorbar()
plt.show()

#%%
plt.imshow()