#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 20:50:47 2018

@author: kyungdoehan
"""

import os
import datetime
import flopy
#%% Model paths variables
class model_paths:
    def __init__(self):
        self.now = datetime.datetime.now()
        self.fdirmfexe = os.getcwd()
        self.fnmmfexe = '/mfnwt'
        self.fdirmodel = os.getcwd() + "/models/" + str(self.now.month)\
        + str(self.now.day) + "_" + str(self.now.hour) + str(self.now.minute)
        self.fnmmodel = 'T2H' + "_" + str(self.now.month) + str(self.now.day)\
        + "_" + str(self.now.hour) + str(self.now.minute)
def paths():
    return model_paths()
#%%
class model_packages:
    def __init__(self, nz, ny, nx, dy, dx, top, bot, ibound, h_ini, hfb_pair,\
                 conductance, laytyp, hk, vkcb_data, rech, fdirmodel, fnmmodel\
                 , nlakes, ilak, ghb, imodel):
        self.mf = flopy.modflow.Modflow\
        (modelname = fnmmodel, \
         exe_name = (paths().fdirmfexe + paths().fnmmfexe), \
         model_ws = fdirmodel, version = 'mfnwt')
        self.mfdis = flopy.modflow.ModflowDis(self.mf, \
                                     nlay = nz, nrow = ny, ncol = nx,\
                                     delr = dy, delc = dx,\
                                     top = top, botm = bot,\
                                     perlen = 1, nstp = 1,\
                                     steady = True, itmuni = 5, lenuni = 2)
        self.mfbas = flopy.modflow.ModflowBas(self.mf, ibound = ibound,\
                                              strt = h_ini, ichflg = True)
        # strt = starting heads; 
        if imodel == 1 or imodel == 2:
            self.mfhfb = flopy.modflow.ModflowHfb(self.mf, hfb_data = hfb_pair)
        self.mfdrn = flopy.modflow.ModflowDrn(self.mf, stress_period_data={0: conductance})
        # drn stage = top of the model grid (=topography)
        if imodel == 1 or imodel == 2:
            self.mfupw = flopy.modflow.ModflowUpw(self.mf, laytyp=laytyp,\
                                                  layavg = 0, chani = 0, hani = 1, \
                                                  layvka = 0.1, laywet = 0,\
                                                  ipakcb = 1, iphdry = 0, hk = hk,\
                                                  vkcb = vkcb_data)
        elif imodel == 3 or imodel == 4:
            self.mfupw = flopy.modflow.ModflowUpw(self.mf, laytyp=laytyp,\
                                                  layavg = 0, chani = 0, hani = 1, \
                                                  layvka = 0.1, laywet = 0,\
                                                  ipakcb = 1, iphdry = 0, hk = hk)            
        self.mfrch = flopy.modflow.ModflowRch(self.mf, rech = rech, nrchop = 3)
        self.mfnwt = flopy.modflow.ModflowNwt(self.mf, linmeth = 1, iprnwt = 1,\
                                 headtol = 1e-4, options = 'COMPLEX')
        if ilak == 1:
            self.mfghb = flopy.modflow.ModflowGhb(self.mf, ipakcb = 1, stress_period_data = ghb)       
        self.output = {(0,0):['save budget', 'save head']}
        self.mfoc = flopy.modflow.ModflowOc(self.mf, stress_period_data=self.output)
        print("PACKAGES console> MODFLOW-NWT input generated")
def packages(nz, ny, nx, dy, dx, top, bot, ibound, h_ini, hfb_pair,\
                 conductance, laytyp, hk, vkcb_data, rech, fdirmodel, fnmmodel\
                 , nlakes, ilak, ghb, imodel):
    return model_packages(nz, ny, nx, dy, dx, top, bot, ibound, h_ini, hfb_pair,\
                 conductance, laytyp, hk, vkcb_data, rech, fdirmodel, fnmmodel\
                 , nlakes, ilak, ghb, imodel)
