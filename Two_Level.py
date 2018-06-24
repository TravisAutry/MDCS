#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 24 15:23:32 2018

@author: travis
"""

import MDCS_Diagrams as MD
import numpy as np
import matplotlib.pyplot as plt


## Define constants and time/frequency/energy axes
#global tMesh tauMesh Tdelay

#Define Energy Range Required
meV2Hz = 241.79895E9*2*np.pi;      # hbar# Define time/frequency [axes and mesh]
meV2THz = meV2Hz*1E-12
Size = 512;   # Equal size along both directions i.e. the number of points
center_Eg = 1550
EgyRange = 10;    # omega_t half-range in meV 

# Define energy axes [Used for plotting]
EgyAxis = np.linspace(-EgyRange,EgyRange,Size+1);   # Spectrum center point should exist
EgyAxis = EgyAxis[0:-1];   # omega_t [Resize]
EgyAxis2 = np.transpose(EgyAxis)    # omega_tau
# Define frequency in Hz [To obtain time axis]
FAxis = EgyAxis*meV2THz/(2*np.pi)
tStep = 1/(FAxis[-1] - FAxis[0])
tAxis =  np.arange(0,tStep*(Size),tStep)
tauAxis = np.transpose(tAxis)
tMesh = np.tile(tAxis,(Size,1))
tauMesh = np.transpose(np.tile(tauAxis,(Size,1)))
Tdelay = tStep;

EgyAxis = EgyAxis 
EgyAxis2 = EgyAxis2 


D1 = MD.Diagram(name = 'D1')
D1.updateParams( [1,0,0,0,.4,.4,.4,.8,.001,.001])
#
D2 = MD.Diagram(name = 'D2')
D2.updateParams( [1,0,0,0,.4,.4,.4,.8,.001,.001])
#
#
#
#
Spectra1 = MD.Spectra(D1,D2,name = '1')
Spectra1.FWM_S1(tMesh,tauMesh,Tdelay)
#
plt.figure();
plt.contourf(tMesh,tauMesh,np.abs(Spectra1.FWMS1))
#
#
plt.figure();
plt.contourf(EgyAxis,EgyAxis,np.abs(Spectra1.FWMS1_E)/np.amax(np.abs(Spectra1.FWMS1_E)),levels = np.linspace(0,1,21))
plt.colorbar()
#
#
#
