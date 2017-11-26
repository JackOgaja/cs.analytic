#!/usr/bin/env python
#
# -*- coding: utf-8 -*-
#
# file:     hos_utils.py
# created:  2014/10/10  Jack Ogaja <jack.ogaja'at'gmail.com>
# Copyright (C) 2014 Jack Ogaja
# See LICENSE for details
#
# purpose:  Plotting accuracy analyses of different Central Difference
#           advection schemes implemented in a NWP and RCM model
# usage:    "python schemes.py" from shell or via
#           "./hosanalytics.py" after making this script executable
##
# Attributes: hosFunctions
# (a) Phase errors:
#            omegaSeond      =>Symmetric Second order (QCHOS)
#            omegaSecondTrad =>Traditional Second order (HOS)
#            omegaFourth     =>Symmetric Fourth order (QCHOS)
#            omegaFourthTrad =>Traditional Fourth order (HOS)
#            omegaSixth      =>Symmetric Sixth order (QCHOS)
#            omegaSixthTrad  =>Traditional Sixth order (HOS)
# (b) Group velocities
#            groupVSeond      =>Symmetric Second order (QCHOS)
#            groupVSecondTrad =>Traditional Second order (HOS)
#            groupVFourth     =>Symmetric Fourth order (QCHOS)
#            groupVFourthTrad =>Traditional Fourth order (HOS)
#            groupVSixth      =>Symmetric Sixth order (QCHOS)
#            groupVSixthTrad  =>Traditional Sixth order (HOS)
#--------------------------------------------------------------#
# >> modules
from __future__ import absolute_import

import numpy as np
import math as mt

#--- The following can be bundled into a config file ---#
import os, sys, inspect
libPath = os.path.abspath('hosUtilities')
sys.path.append(libPath)
#------------------------#

import hosFunctions as hf
import hosGraphics as hg
#----------------------------------------------------------------------------80
# >> defintions

#---- Wavenumber ----#
c=10.
dx=2.
k=np.linspace(0.00095,np.pi/dx,50)
#--alias Error settings---#
#m=np.pi/dx
#m=8.*np.pi/9./dx
#m=3.*np.pi/4./dx
m=np.pi/2./dx
#m=np.pi/20./dx
k2=np.linspace(-np.pi/dx,np.pi/dx,100)
k1=np.linspace(-np.pi/dx,np.pi/dx,100)
#+++k1=m-k2
k3=k1+k2
kmax=np.pi/dx

if m >= 0.0:
    b=k1-kmax
    #b=m-kmax
else:
    b=k1+kmax
    b=m+kmax

print(' K3[4] = {0}, K3[83] = {1}'.format(k3[4],k3[83]))
kt=[]
for i in range(len(k3)):
    if k3[i] > kmax:
        kt.append(-2.0*kmax + k3[i])
    if k3[i] < -kmax:
        kt.append(2.0*kmax + k3[i])

print('LEN(KT) = {0}'.format(len(kt)))
print(kt)
#--------------------------#
#+==============+#
# Error objects +#
#+==============+#
hos     = hf.hosSymbolic()
hos.dx  = dx
hos.c   = c
hos1    = hos.PhaseError(k)
hos2    = hos.GroupVelocity(k)
hos3    = hos.aliasError(k1,k2)
#+--------------+#

(oSec,oSecTrad,oFourth,oFourthTrad,oSixth,oSixthTrad) = hos1
(vSec,vSecTrad,vFourth,vFourthTrad,vSixth,vSixthTrad) = hos2
(aSec,aSecTrad,aThird,aFou,aFouTrad,aFifth,aSix,aSixTrad)\
 =[abs(scheme) for scheme in hos3]

#hos4\
# = [hos.aliasError(k1[i],k2)
#    for i in range(len(k1))]
#(aSec2d,aSecTrad2d,aThird2d,aFou2d,aFouTrad2d,
# aFifth2d,aSix2d,aSixTrad2d)\
#  =[[hos4[m][n] for m in range(len(k1))]
#                      for n in range(len(hos4[0]))]

hos4\
 = [hos.aliasError(kt[i],k2)
    for i in range(len(kt))]
(aSec2d,aSecTrad2d,aThird2d,aFou2d,aFouTrad2d,
 aFifth2d,aSix2d,aSixTrad2d)\
  =[[hos4[m][n] for m in range(len(kt))]
                      for n in range(len(hos4[0]))]

bp=[[b[i][j] for i in range(len(kt))] for j in range(len(k2))]

#+=================================+#
#+ Plots of the analytic functions +#
#+=================================+#
#-- Plot settings/Attributes---#
plt = hg.hosPlot()  # call the plot utility
plt.title='Some plot'
plt.styles  = ['b-.','r--','k+'] #['b-.','r--','k+','kx','g-','k:']
#plt.lwidth=[2,2,2,2,2]
plt.labels  = ['2nd order']
plt.savePlot= True
plt.plotName= 'MyPlot'
plt.plotDir = 'Dietz'         # Folder in the cwd
#plt.legLoc  = 'lower center'  # Legend location
plt.lblsize = 14
plt.ttlsize = 16
plt.minorticks=True
plt.lgframe=True
plt.lgFontSize=12
#plt.ymax=0.14
#plt.xmax=0.14
plt.ylabel='Effecitve wavenumber $k_{eff}$'
plt.xlabel='wave number $k$'
#plt.vline=True
#plt.hline=True
#plt.plotLine(k,oSec,oFourth,oFourthTrad,oSixth,oSixthTrad,refLine=(k,'Toto'))

plt2=plt
plt2.title='You Know'
#plt2.ymax=13.
#plt2.ymin=-23.
plt2.plotName='MyPlot2'
plt2.lwidth=[2,2,2,2,2]
#plt2.plotLine(k,vSec,vFourth,vFourthTrad,vSixth,vSixthTrad,hLine=(c,'Babe'))

plt3=plt
plt3.savePlot=False
plt3.lwidth=[]
#plt3.xmin=-np.pi/dx
#plt3.xmax=np.pi/dx

#plt3.plotLine(k2,aSec,aSecTrad,aThird,aFou,aFouTrad,aFifth,
#              aSix,aSixTrad,hLine=0.0,refLine=(abs(k2),'refline'))

#print(aSec2d)
#print(aThird2d)
plt4=plt
plt4.plotContour(k2,kt,aSec2d)
plt4.plotContour(k2,kt,bp)

#==== Some test =====
#plt5=plt
#plt5.testContour()
