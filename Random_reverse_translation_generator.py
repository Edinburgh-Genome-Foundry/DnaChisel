#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 09:52:54 2022

@author: weiqiyao
"""
from dnachisel import *
from dnachisel import biotools
import flametree

protein = "MAWWL*"
test = biotools.reverse_translate(protein,True,table = 'Bacterial')
print(test)
# GGTCATATTTTAAAAATGTTTCCT