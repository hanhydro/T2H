#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 14:47:47 2018

@author: kyungdoehan
"""

import flopy
import os 
class paths:
    def __init__(self):
        self.fdirmpexe = os.getcwd() + "/exe/"
        self.fnmpexe = 'mp6'
        
class particle_tracking:
    def __init__(self, fdirmpexe, )