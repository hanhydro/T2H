#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 19:00:49 2018

@author: kyungdoehan
"""
import pandas as pd
import os
#%% check every input file is in the path 
class dataclasses:
    def __init__(self, Ma):
        self.dir = os.getcwd() + "/tisc_output/" 
        self.topo = os.path.isfile(self.dir + "topo_" + str(Ma) + "0Ma.txt")
        self.bdtopo = os.path.isfile(self.dir + "bdtopo_" + str(Ma) + "0Ma.txt")
        self.sed = os.path.isfile(self.dir + "sedthick_" + str(Ma) + "0Ma.txt")
        self.lake = os.path.isfile(self.dir + "lakes_" + str(Ma) + "0Ma.csv")
        self.fault = os.path.isfile(self.dir + "fault_" + str(Ma) + "0Ma.txt")
def data(Ma):
    return dataclasses(Ma)
#%% importing topography data from topo_*.txt
class Topo_data:
    def __init__(self, Ma):
        self.dir = os.getcwd()
        self.df = pd.read_csv(self.dir + "/tisc_output/topo_" + str(Ma) + "0Ma.txt",\
                 sep=" ", header=None)
        self.df.columns = ["x", "y", "elev"]
        self.x = self.df.loc[:, "x"].values.astype('int')
        self.y = self.df.loc[:, "y"].values.astype('int')
        self.z = self.df.loc[:, "elev"].values
def topo_read(Ma):
    return Topo_data(Ma)
#%% reading sediment thickness data from file
class Sed_data:
    def __init__(self, Ma):
        self.dir = os.getcwd()
        self.df_sed = pd.read_csv(self.dir + "/tisc_output/sedthick_" \
                                  + str(Ma) + "0Ma.txt",\
                                  sep=" ", header=None)
        self.df_sed.columns = ["x", "y", "thickness"]
        self.x = self.df_sed.loc[:, "x"].values.astype('int')
        self.y = self.df_sed.loc[:, "y"].values.astype('int')
        self.z = self.df_sed.loc[:, "thickness"].values
def sed_read(Ma):
    return Sed_data(Ma)
#%% 
class Fault_data:
    def __init__(self, Ma):
        self.dir = os.getcwd()
        self.df_fault = pd.read_csv(self.dir + "/tisc_output/fault_"\
                                    + str(Ma) + "0Ma.txt",\
                                    sep=" ", header=None)
        self.df_fault.columns = ["x", "y", "depth"]
        self.x = self.df_fault.loc[:, "x"].values.astype('int')
        self.y = self.df_fault.loc[:, "y"].values.astype('int')
        self.z = self.df_fault.loc[:, "depth"].values

def fault_read(Ma):
    return Fault_data(Ma)
#%%
class Fault_data_rearrange:
    def __init__(self, Ma):
        self.dir = os.getcwd()
        self.df_fault = pd.read_csv(self.dir + "/tisc_output/fault_"\
                                    + str(Ma) + "0Ma.txt",\
                                    sep=" ", header=None)
        self.df_fault.columns = ["x", "y", "depth"]
        self.fault_arr = self.df_fault.sort_values(["x", "y", "depth"], ascending=[True, True, False]) 
        self.x = self.fault_arr.loc[:, "x"].values.astype('int')
        self.y = self.fault_arr.loc[:, "y"].values.astype('int')
        self.z = self.fault_arr.loc[:, "depth"].values

def fault_rearrange(Ma):
    return Fault_data_rearrange(Ma)