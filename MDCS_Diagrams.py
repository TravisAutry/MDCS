# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 16:42:25 2016

@author: Travis Autry

"""


import numpy as np
import pandas as pd

''' This program generates 2D spectra 

# Functions that evaluate FWM signal for different diagrams are called
S1_Diagram = Rephasing Pulse Sequence
S2_Diagram = Non-Rephasing Pulse Sequence
S3_Diagram = Double Quantum Pulse Sequence diagrams

All simulations are done in a rotating frame.  So if your signal is at 1500meV. Set your
center_freq = 1500 meV.


'''


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


class Spectra():
        def __init__(self,*args, name = []):
            print ('args: ', args)
            self.args = args
            self.name = name
            self.FWMS1 = np.complex(0,1)*0
            self.FWMS2 =  np.complex(0,1)*0
            self.FWMS3 =  np.complex(0,1)*0
            
            self.FWMS1_E = np.complex(0,1)*0
            self.FWMS2_E =  np.complex(0,1)*0
            self.FWMS3_E =  np.complex(0,1)*0
            
            self.dlist = list()
            for i in range(np.shape(self.args)[0]):
                self.dlist.append(self.args[i].name)
            
            self.diag_dict = dict(zip(self.dlist,self.args))
            
        def FWM_S1(self,tMesh,tauMesh,Tdelay):
            I = np.complex(0,1)

            for i in self.diag_dict.keys():
                FWM = self.diag_dict[i].S1_Diagram(tMesh,tauMesh,Tdelay)
                self.FWMS1 = self.FWMS1 + FWM
            
            self.FWMS1 = self.FWMS1*(-1*I)
            self.FWMS1_E = self.FFT2D(self.FWMS1)
                
        def FWM_S2(self,tMesh,tauMesh,Tdelay):
            I = np.complex(0,1)
       
            for i in self.diag_dict.keys():
                FWM = self.diag_dict[i].S2_Diagram(tMesh,tauMesh,Tdelay)
                self.FWMS2 = self.FWMS2 + FWM
           
            self.FWMS2 = self.FWMS2*(-1*I)
            self.FWMS2_E = self.FFT2D(self.FWMS2)

        def FWM_S3(self,tMesh,tauMesh,Tdelay):
            I = np.complex(0,1)
 
            for i in self.diag_dict.keys():
                FWM = self.diag_dict[i].S3_Diagram(tMesh,tauMesh,Tdelay)
                self.FWMS3 = self.FWMS3 + FWM
           
            self.FWMS3 = self.FWMS3*(-1*I)
            self.FWMS3_E = self.FFT2D(self.FWMS3)


        def FFT2D(self,A):
            
            B = np.fft.fft2(A)
            B = np.fft.fftshift(np.fft.fftshift(B,axes=0),axes=1)
            
            return B


class Diagram():
     

    def __init__(self,name = []):
       self.meV2Hz = 241.79895E9*2*np.pi;      # hbar# Define time/frequency [axes and mesh]
       self.meV2THz = self.meV2Hz*1E-12  
         
       self.name = name
       self.dipole = 0
       self.w_tau = 0
       self.w_T = 0
       self.w_t = 0
       self.g_tau = 0
       self.g_T = 0
       self.g_t = 0
       self.sigma = 0
       self.delta = 0
       self.delta_g = 0
       self.Diag_Param_meV2Hz = []
       self.Diag_Param = pd.DataFrame(columns = ['meV'],data = [self.dipole,self.w_tau,self.w_T,self.w_t,self.g_tau,self.g_T,self.g_t,self.sigma,self.delta,self.delta_g],index = ['Dipole^4','W_tau','W_T','W_t','g_tau','g_T','g_t','Sigma','Delta','Delta_g'])
       
       for i in range(self.Diag_Param.shape[0]-1):
           self.Diag_Param['THz'] = self.Diag_Param['meV'][i+1]*self.meV2THz

    def updateParams(self,values = []):
        
       #if no updates values then update defaults to the current values
       if values == []:
           values = self.Diag_Param['meV'].values
    
        #if values given then update 'meV' and then convert to THz
       self.Diag_Param['meV'] = values 
       for i in range(self.Diag_Param.shape[0]):
           if i != 0:
               self.Diag_Param['THz'].values[i] = self.Diag_Param['meV'].values[i]*self.meV2THz
           else:
               self.Diag_Param['THz'].values[i] = self.Diag_Param['meV'].values[i]
       print('*****   {0} is {1} *****'.format('Name:', self.name))
       print(self.Diag_Param)
       
     

    def S1_Diagram(self,tMesh,tauMesh,Tdelay):
        
        '''This does is the diagram for a rephasing pulse sequence i.e. photon echo '''
        Amp = self.Diag_Param['THz'][ 'Dipole^4']
        w1 = self.Diag_Param['THz']['W_tau']#; %freq during tau
        w2 = self.Diag_Param['THz']['W_T']##; %freq during T
        w3 = self.Diag_Param['THz']['W_t']#]#;  % freq during t.
        g1 = self.Diag_Param['THz']['g_tau']##; %decay during tau
        g2 = self.Diag_Param['THz']['g_T']#Diag_Param['THz']['W_tau']##; %decay during T
        g3 = self.Diag_Param['THz']['g_t']#; %decay during t
        Delta = self.Diag_Param['THz']['Delta']#
        Delta_g = self.Diag_Param['THz']['Delta_g']#
        
        w3 = w3 + Delta
        g3 = g3 + Delta_g
        
        sigma1 = self.Diag_Param['THz']['Sigma'];
        
        I = complex(0,1);
        exp = np.exp
        FWM_t = Amp*I*exp(- tauMesh* (g1+I*w1))*exp(- tMesh* (g3-I*w3))* exp((-1/2)*((sigma1*tMesh)**2 + (1*sigma1*tauMesh)**2- (2*tauMesh*tMesh*sigma1**2)))*exp(-g2*Tdelay);
        
        
        return FWM_t
    
    def S2_Diagram(self,tMesh,tauMesh):
        
        Amp = self.Diag_Param['THz'][ 'Dipole^4']
        w1 = self.Diag_Param['THz']['W_tau']#; %freq during tau
        w2 = self.Diag_Param['THz']['W_T']##; %freq during T
        w3 = self.Diag_Param['THz']['W_t']#]#;  % freq during t.
        g1 = self.Diag_Param['THz']['g_tau']##; %decay during tau
        g2 = self.Diag_Param['THz']['g_T']#Diag_Param['THz']['W_tau']##; %decay during T
        g3 = self.Diag_Param['THz']['g_t']#; %decay during t
        Delta = self.Diag_Param['THz']['Delta']#
        Delta_g = self.Diag_Param['THz']['Delta_g']#
        
        w3 = w3 + Delta
        g3 = g3 + Delta_g
        
        sigma1 = self.Diag_Param['THz']['Sigma'];
        
        
        I = complex(0,1);
        exp = np.exp
    
        FWM_t = Amp*I*exp(- tauMesh* (g1-I*w1))*exp(- tMesh* (g3-I*w3))* exp((-1/2)*((sigma1*tMesh)**2 + (1*sigma1*tauMesh)**2+ (2*tauMesh*tMesh*sigma1**2)))*exp(-g2*Tdelay);
        
        
        return FWM_t 
    
    def S3_Diagram(self,tMesh,tauMesh):
        
        Amp = self.Diag_Param['THz'][ 'Dipole^4']
        w1 = self.Diag_Param['THz']['W_tau']#; %freq during tau
        w2 = self.Diag_Param['THz']['W_T']##; %freq during T
        w3 = self.Diag_Param['THz']['W_t']#]#;  % freq during t.
        g1 = self.Diag_Param['THz']['g_tau']##; %decay during tau
        g2 = self.Diag_Param['THz']['g_T']#Diag_Param['THz']['W_tau']##; %decay during T
        g3 = self.Diag_Param['THz']['g_t']#; %decay during t
        Delta = self.Diag_Param['THz']['Delta']#
        Delta_g = self.Diag_Param['THz']['Delta_g']#
        
        w3 = w3 + Delta
        g3 = g3 + Delta_g
        
        sigma1 = self.Diag_Param['THz']['Sigma'];
        
        
        I = complex(0,1);
        exp = np.exp
    
        FWM_t = Amp*I*exp(- tauMesh* (g1-I*w1))*exp(- tMesh* (g3-I*w3))* exp((-1/2)*((sigma1*tMesh)**2 + (1*sigma1*tauMesh)**2+ (2*tauMesh*tMesh*sigma1**2)))*exp(- Tdelay* (g2-I*w2));
        
        
        return FWM_t




#
#                        
#                       
#
#D1 = Diagram(name = 'D1')
#D1.updateParams( [1,0,0,0,.4,.4,.4,.8,.001,.001])
#
#D2 = Diagram(name = 'D2')
#D2.updateParams( [1,0,0,0,.4,.4,.4,.8,.001,.001])
#
#
#
#
#Spectra1 = Spectra(D1,D2,name = '1')
#Spectra1.FWM_S1(tMesh,tauMesh,Tdelay)
#
#plt.figure();
#plt.contourf(tMesh,tauMesh,np.abs(Spectra1.FWMS1))
#
#
#plt.figure();
#plt.contourf(EgyAxis,EgyAxis,np.abs(Spectra1.FWMS1_E)/np.amax(np.abs(Spectra1.FWMS1_E)),levels = np.linspace(0,1,21))
#plt.colorbar()
#
#
#
