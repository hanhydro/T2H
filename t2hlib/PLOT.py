#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 14:58:45 2018

@author: kyungdoehan
"""
# a = 3-D head array to plot
# masked_values: iterable of floats, ints
# head = 3-D array to set top of patches to the min. of the top of a layer
# or the head value. Used to create patches that conform to WT elevations

import flopy
import numpy as np
from flopy.export import vtk as fv
import matplotlib.pyplot as plt
import plotly.plotly as py
import plotly.figure_factory as ff

#class getextent:
#    def __init__(self):
#        
#def gex():
#    return getextent()
#%%
class patches:
    def __init__(self, mfxsect, hk, Kconst, vkcb_data, hfb_info, nz,\
                 ny, nx, imodel):
        # sediment patches
        self.hkk = np.zeros((nz, ny, nx), dtype=np.float32)
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    if hk[k, j, i] > Kconst:
                        self.hkk[k,j,i] = hk[k,j,i]
        #self.plotarray = mfxsect.plot_array(self.hkk, masked_values =[Kconst], head=None, color = 'y', zorder=1)
        # fault patches
        if imodel == 1 or imodel == 2:
            self.pfault = mfxsect.plot_array(vkcb_data, masked_values = [1], color='r')
        #self.phfb = mfxsect.plot_array(hfb_info, masked_values = [1], color='r')
        #self.bcarray = mfxsect.plot_bc(ftype = 'GHB', package = 'mfghb', kper=0, color = 'g')

def pat(mfxsect, hk, Kconst, vkcb_data, hfb_info, nz, ny, nx, imodel):
    return patches(mfxsect, hk, Kconst, vkcb_data, hfb_info,\
                   nz, ny, nx, imodel)

#%%
class watertable:
    def __init__(self, mfxsect, top, a, nz, ny, nx):
        self.wt = a[0, :, :] - top # for bird-eye view, 2-D
        self.wt_m = np.zeros((ny, nx), dtype=np.float32)
        self.wtarr = np.zeros((nz, ny, nx), dtype=np.float32) # for 3-, xsects
        for j in range(ny):
            for i in range(nx):
                self.wtarr[0, j, i] = self.wt[j, i] 
                if self.wt[j, i] <= 0:
                    self.wt_m[j, i] = 0           
        self.wtf = np.zeros((1, ny, nx), dtype=np.float32)
        for j in range(ny):
            for i in range(nx):
                if a[0, j, i] < 0:
                    for k in range(nz):
                        if a[k, j, i] >= 0:
                            self.wtf[0, j, i] = a[k, j, i]
                            break
                if a[0, j, i] >= 0:
                    self.wtf[0, j, i] = a[0, j, i]
        self.wtf1 = np.zeros((nz-1, ny, nx), dtype=np.float32)
        self.depth_wt = top - self.wtf
        self.wtf = np.vstack((self.wtf, self.wtf1))
#        self.plarr = mfxsect.plot_array(self.wtarr, masked_values=[0, -1], head=a)
        #self.plsurf = mfxsect.plot_surface(self.wtf[0, :, :], color='b')
        #self.plsurf = mfxsect.plot_array(self.wtf[0, :, :], color='b')
        #self.plarr1 = mfxsect.plot_fill_between(self.wtarr, masked_values=[0, -1], head=a)

        
def wtable(mfxsect, top, a, nz, ny, nx):
    return watertable(mfxsect, top, a, nz, ny, nx)
        
        
#%%
class mfxsect_data:
    def __init__(self, mf, mfdis, target, axheadxsect):
        if target >= 0:
            self.mfxsect = flopy.plot.crosssection.ModelCrossSection\
            (ax = axheadxsect, model = mf, dis = mfdis, \
             line = {'Row': target},\
             xul = 0, yul = None, rotation = 0.0)
        elif target < 0:
            self.mfxsect = flopy.plot.crosssection.ModelCrossSection\
            (ax = axheadxsect, model = mf, dis = mfdis, \
             line = {'Column': np.abs(target)},\
             xul = 0, yul = None, rotation = 0.0)            
            
        # xul = x-coordinate of upper left corner
        # yul = y-coordinate of upper left corner
def fmfxsect(mf, mfdis, target, axheadxsect):
    return mfxsect_data(mf, mfdis, target, axheadxsect)
        
class head_data:
    def __init__(self, mf, fdirmodel):
        self.mfhds = flopy.utils.HeadFile(fdirmodel + "/" + mf.name + '.hds')
        self.a = self.mfhds.get_data(kstpkper = (0, 0))
def head(mf, fdirmodel):
    return head_data(mf, fdirmodel)
# head contour plot
class head_contour:
    def __init__(self, mfxsect, a):
        self.headcontour = mfxsect.contour_array(a, linewidths=0.5, \
                                    levels=np.arange(0,np.max(a),5))
def headc(mfxsect, a):
    return head_contour(mfxsect, a)
# cell budget calculation to get flux vectors
class cellBudget:
    def __init__(self, fdirmodel, mf):
        self.mfcbc = flopy.utils.CellBudgetFile\
        (fdirmodel + "/" + mf.name + '.cbc')
        self.times = self.mfcbc.get_times()
        self.qx = self.mfcbc.get_data\
        (text='flow right face', totim = self.times[-1])[0]
        self.qy = self.mfcbc.get_data\
        (text='flow front face', totim = self.times[-1])[0]
        self.qz = self.mfcbc.get_data\
        (text='flow lower face', totim = self.times[-1])[0] # lower face
def cbc(fdirmodel, mf):
    return cellBudget(fdirmodel, mf)
# avg. flux calculations
class cellCenters:
    def __init__(self, qz, qy, qx, nz, ny, nx):
        self.qx_avg = np.empty(qx.shape, dtype=qx.dtype)
        self.qx_avg[:, :, 1:] = 0.5 * (qx[:, :, 0:nx-1] + qx[:, :, 1:nx])
        self.qx_avg[:, :, 0] = 0.5 * qx[:, :, 0]
        self.qy_avg = np.empty(qy.shape, dtype=qy.dtype)
        self.qy_avg[:, 1:, :] = 0.5 * (qy[:, 0:ny-1, :] + qy[:, 1:ny, :])
        self.qy_avg[:, 0, :] = 0.5 * qy[:, 0, :]
        self.qz_avg = np.empty(qz.shape, dtype=qz.dtype)
        self.qz_avg[1:, :, :] = 0.5 * (qz[0:nz-1, :, :] + qz[1:nz, :, :])
        self.qz_avg[0, :, :] = 0.5 * qz[0, :, :]
def qavg(qz, qy, qx, nz, ny, nx):
    return cellCenters(qz, qy, qx, nz, ny, nx)

class discharge:
    def __init__(self, qx_avg, qy_avg, qz_avg, mfxsect, a):
        self.disch = mfxsect.plot_discharge(frf=qx_avg, fff=qy_avg,\
                                            flf=qz_avg, head=a,\
                                             kstep = 5, hstep = 5,\
                                             normalize=False, color='r')
def disch(qx_avg, qy_avg, qz_avg, mfxsect, a):
    return discharge(qx_avg, qy_avg, qz_avg, mfxsect, a)

#, headwidth=.01, headlength=0.05, headaxislength=.01, width=0.002)
#%% quiver plot for cross sections
class quiv:
    def __init__(self, X, Z, qx_avg, qy_avg, qz_avg, target, iskip, axheadxsect):
        def qlog(x):
            return np.sign(x) * np.log(np.abs(x)+1)
        
        if target >= 0:
            self.HEADL = 4
            self.angle = np.arctan2(-qz_avg[::iskip, target, ::iskip],\
                                    qx_avg[::iskip, target, ::iskip])*180/np.pi
            self.quiver = axheadxsect.quiver(X[::iskip, ::iskip],\
                                             Z[::iskip, ::iskip], \
                                             qlog(qx_avg[::iskip, target,\
                                                         ::iskip]),\
                                             qlog(-qz_avg[::iskip, target,\
                                                          ::iskip]),\
                                             color='k', angles = self.angle, \
                                             headaxislength=self.HEADL,\
                                             units = 'x', scale = 0.008)
        elif target < 0:
            self.HEADL = 4
            self.angle = np.arctan2(-qz_avg[::iskip, target, ::iskip],\
                                    qy_avg[::iskip, target, ::iskip])*180/np.pi
            self.quiver = axheadxsect.quiver(X[::iskip, ::iskip],\
                                             Z[::iskip, ::iskip], \
                                             qlog(qy_avg[::iskip, target,\
                                                         ::iskip]),\
                                             qlog(-qz_avg[::iskip, target,\
                                                          ::iskip]),\
                                             color='k', angles = self.angle, \
                                             headaxislength=self.HEADL,\
                                             units = 'x', scale = 0.008)
def quiver(X, Z, qx_avg, qy_avg, qz_avg, target, iskip, axheadxsect):
    return quiv(X, Z, qx_avg, qy_avg, qz_avg, target, iskip, axheadxsect)
#%% Streamline plot for cross-sections
class streamlines:
    def __init__(self, X, Z, qx_avg, qz_avg, target):
        fig, ax = plt.subplots(figsize=(40, 5))
        self.strm = ax.streamplot(X, Z, qx_avg[:, target, :], -qz_avg[:, target, :])#, color = X, linewidth =1)
        fig.colorbar(self.strm.lines)
def sline(X, Z, qx_avg, qz_avg, target):
    return streamlines(X, Z, qx_avg, qz_avg, target)
#%% water table fluctuation plot        
class WTF:
    def __init__(self, nz, ny, nx, a, top):
        self.wtf = np.zeros((1, ny, nx), dtype=np.float32)
        for j in range(ny):
            for i in range(nx):
                if a[0, j, i] < 0:
                    for k in range(nz):
                        if a[k, j, i] >= 0:
                            self.wtf[0, j, i] = a[k, j, i]
                            break
                if a[0, j, i] >= 0:
                    self.wtf[0, j, i] = a[0, j, i]
        self.wtf1 = np.zeros((nz-1, ny, nx), dtype=np.float32)
        self.depth_wt = top - self.wtf
        self.wtf = np.vstack((self.wtf, self.wtf1))

def wt_f(nz, ny, nx, a, top):
    return WTF(nz, ny, nx, a, top)

class vtkout:
    def __init__(self, mf, a, qx_avg, qy_avg, qz_avg,\
                 hk, vkcb_data, hfb_info, wtf, lakarr):
        self.vtkout = "model.vtu"
        self.vtkobj = fv.Vtk(self.vtkout, mf)
        self.vtkobj.add_array("head", a) # Save head array from MF output
        self.vtkobj.add_array("qx_avg", qx_avg) 
        self.vtkobj.add_array("qy_avg", qy_avg)
        self.vtkobj.add_array("qz_avg", qz_avg)
        self.vtkobj.add_array("hk", hk)
        self.vtkobj.add_array("vkcb_data", vkcb_data)
        self.vtkobj.add_array("hfb_data", hfb_info)
        self.vtkobj.add_array("water table", wtf)
        self.vtkobj.add_array("lake", lakarr)
        self.vtkobj.write(ibound_filter=True)
def vtk_export(mf, a, qx_avg, qy_avg, qz_avg, hk, vkcb_data, hfb_info, wtf, lakarr):
    return vtkout(mf, a, qx_avg, qy_avg, qz_avg, hk, vkcb_data, hfb_info, wtf, lakarr)

class discharge_location:
    def __init__(self, a, top, ny, nx, qz_avg):
        self.discharge = np.zeros([ny, nx])
        for j in range(ny):
            for i in range(nx):
                if a[0, j, i] > top[j, i]:
                    self.discharge[j, i] = -qz_avg[0, j, i]
def dischloc(a, top, ny, nx, qz_avg):
    return discharge_location(a, top, ny, nx, qz_avg)

class wtcontour:
    def __init__(self, x, y, dx, dy, ny, nx, wtf):
        fig, ax = plt.subplots(figsize=(15,15))
        self.x1 = np.linspace(dx/2, nx*dx - dx/2, nx)
        self.y1 = np.linspace(ny*dy - dy/2, dy/2, ny)
        self.plot = plt.contour(self.x1, self.y1, wtf[0, :, :])
        #self.im = plt.imshow(self.x1, self.y1, wtf[0, :, :], cmap =cm.Set2)
        self.lx = plt.xlabel("x")
        self.ly = plt.ylabel("y")
        self.label = plt.clabel(self.plot, inline = 1, \
                                fmt = "%1.1f", fontsize = 10)
def wtcont(x, y, dx, dy, ny, nx, wtf):
    return wtcontour(x, y, dx, dy, ny, nx, wtf)

class elevation:
    def __init__(self, x, y, nx, dx, ny, dy, top):
        fig, ax = plt.subplots(figsize=(15,15))
        self.x1 = np.linspace(dx/2, nx*dx - dx/2, nx)
        self.y1 = np.linspace(ny*dy - dy/2, dy/2, ny)
        self.plot = plt.contour(self.x1, self.y1, top, zorder=1)
        self.im = plt.imshow(self.x1, self.y1, top, zorder = 0)
        self.lx = plt.xlabel("x")
        self.ly = plt.ylabel("y")
        self.label = plt.clabel(self.plot, inline = 1, \
                                fmt = "%1.1f", fontsize = 10)
        self.bb = plt.colorbar(self.im)
        plt.show()
def elev(x, y, nx, dx, ny, dy, top):
    return elevation(x, y, nx, dx, ny, dy, top)
