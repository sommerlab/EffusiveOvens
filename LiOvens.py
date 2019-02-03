# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 18:06:41 2019

@author: ariel
"""
import numpy as np
import Delambre.MKS as u
from Oven import Oven

OHaraOven = Oven(atom_str="Li",
              diameter=0.25*u.inch,
              length=0.15*u.inch,
              T=(273+435.)*u.K,
              num_apertures=1,
              MFP_model=2)  
              
print "\nO'Hara:"
OHaraOven.report()


print "\nWeld:"
WeldOven = Oven(atom_str="Li",
              diameter=.1*u.mm,
              length=5*u.mm,
              T=(273+525.)*u.K,
              num_apertures=528*1.5,#factor of 1.5 to account for gaps
              MFP_model=1)  
WeldOven.report()    
#
#print "\nWeld2:"
#WeldOven2 = Oven(atom_str="Li",
#              diameter=.1*u.mm/2,
#              length=5*u.mm,
#              T=(273+495.)*u.K,
#              num_apertures=2*528,
#              MFP_model=2)  
#WeldOven2.report()   


print "\nCooper:"
CooperOven = Oven(atom_str="Li",
              diameter=.51*u.mm,
              length=20*u.mm,
              T=(273+370.)*u.K,#370C
              num_apertures=15,
              MFP_model=2)  
CooperOven.report()           

LiOven = Oven(atom_str="Li",
              diameter=.50*u.mm,
              length=8*u.mm,
              T=(273+435.)*u.K,
              num_apertures=np.int((.25*u.inch/(.5*u.mm))**2),
              MFP_model=2)  
              
print "Custom:"
LiOven.report()