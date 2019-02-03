# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 18:06:41 2019

@author: ariel
"""
import numpy as np
import Delambre.MKS as u
pi=np.pi
class Oven:#single species oven
    def __init__(self, atom_str,diameter,T, length=0, num_apertures=1,MFP_model=1):
        self.atom_str=atom_str
        self.diameter=diameter
        self.area = pi*(diameter/2.)**2
        self.T= T
        self.length=0*u.m if length==0 else length
        self.num_apertures=num_apertures
        self._computeClausing()
        self._setMass()
        self.vmp = (2*u.kB*T/self.mass)**0.5
        self.vavg = 2/np.sqrt(pi)*self.vmp
        self._computePressure()
        self.density = self.P/(u.kB * self.T)
        self.MFP_model = MFP_model
        self._computeMFP()
        self._computeFlux()
        
    def _setMass(self):
        if self.atom_str=="Li":
            self.mass = 9.9883414e-27 *u.kg
            
    def _computeClausing(self):
        a=self.diameter/2
        x=self.length/a
        B=3.*x/8.
        self.Clausing = 1./(1.+B)        
        
    def _computePressure(self):
        TK=self.T/u.K
        if self.atom_str=="Li":
            P_Torr = 10**(10.34540 - 8345.574/TK - 0.00008840*TK - 0.68106*np.log10(TK))
            self.P = P_Torr*u.Torr
            
    def _computeMFP(self):
        if self.MFP_model==1:
            atomDia = 3.13e-10*u.m
            self.MFP = 1./(np.sqrt(2)*np.pi*self.density*atomDia**2)
        if self.MFP_model==2:
            C6 = 1400.*u.hartree*u.a0**6
            C6cgs = C6/(u.erg*u.cm**6)
            vcgs = self.vavg/(u.cm/u.s)
            sigma=5e11*(C6cgs/vcgs)**(2./5)*u.cm**2
            self.MFP=1/(self.density*sigma)  
            
    def _computeFlux(self):
        if self.length==0:#thin aperture
            self.peakIntensity = self.num_apertures*self.density*self.vavg*self.area/(4*pi)
            self.flux = pi*self.peakIntensity
        else:#molecular regime!
            self.peakIntensity = self.num_apertures*self.density*self.vavg*self.area/(4*pi)
            self.flux = self.Clausing * pi*self.peakIntensity
        self.peakingFactor = pi*self.peakIntensity/self.flux
    
    def findLifetime(self,totalMass):
        Natoms = totalMass / self.mass
        return Natoms / self.flux
        
    def report(self):
        print "\tNum apertures:",self.num_apertures
        print "\tFlux:", self.flux
        print "\tPeak Intensity:",self.peakIntensity
        print "\tPeaking Factor:",self.peakingFactor
        print "\tMFP (mm):",self.MFP/u.mm
        print "\t 25 gram lifetime (yr):",self.findLifetime(25*u.gram)/u.year
        print "\t 10 gram lifetime (yr):",self.findLifetime(10*u.gram)/u.year
        print "\t Vapor pressure (mTorr):",self.P/u.Torr/1e-3
        
        
