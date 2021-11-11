# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 15:15:12 2021

@author: TEJan
"""

from pcraster import *

os.chdir("C:/Users/TEJan/Documents/Master/Master Thesis/pycatch-master/inputs")

import os
print(os.getcwd())

mergeDem=readmap("mergeDem.map")
aguila(mergeDem)

