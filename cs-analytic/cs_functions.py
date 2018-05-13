# -*- coding: utf-8 -*-
#
# module:    hos_functions
# created:   Jack Ogaja <jack.ogaja'at'gmail.com>
#
# purpose:   module for analytic functions of
#            Conservative HOS
#
#-----------------------------------------------------------------------------80

import numpy as np
import math as mt
#from abc import ABCMeta, abstractmethod
#-----------------------------------------------------#

class hosSymbolic(object):
    """
    Attributes:
    :c       phase velocit
    :dx      space interval
    :preC4a  first coefficient
    :preC4b  second coefficient
    :preS4a  first coefficient
    :preS4b  second coefficient
    :preS6a  third coefficient
    :preS6b  fourth coefficient
    :preS6c  fifth coefficient

    :omegaSeond      Symmetric Second order (QCHOS)
    :omegaSecondTrad Traditional Second order (HOS)
    :omegaFourth     Symmetric Fourth order (QCHOS)
    :omegaFourthTrad Traditional Fourth order (HOS)
    :omegaSixth      Symmetric Sixth order (QCHOS)
    :omegaSixthTrad  Traditional Sixth order (HOS)

    :groupVSeond      Symmetric Second order (QCHOS)
    :groupVSecondTrad Traditional Second order (HOS)
    :groupVFourth     Symmetric Fourth order (QCHOS)
    :groupVFourthTrad Traditional Fourth order (HOS)
    :groupVSixth      Symmetric Sixth order (QCHOS)
    :groupVSixthTrad  Traditional Sixth order (HOS)

    """

    c  = 10.
    dx = 10.
    preC4a = 4./3.
    preC4b = 1./6.
    preS4a = 9./8.
    preS4b = 1./24.
    preC6a = 3./2.
    preC6b = 3./10.
    preC6c = 1./30.
    preS6a = 150./128.
    preS6b = 25./(128.*3.)
    preS6c = 3./(128.*5.)

    def __calcconst(self):
        coeff= self.c/self.dx
        return coeff

    def __calcsin(self):
        a = np.sin(self.k*self.dx)
        b = np.sin(2.*self.k*self.dx)
        c = np.sin(3.*self.k*self.dx)
        d = np.sin(5.*self.k*self.dx)
        return (a,b,c,d)

    def __calccos(self):
        a = np.cos(self.k*self.dx)
        b = np.cos(2.*self.k*self.dx)
        c = np.cos(3.*self.k*self.dx)
        d = np.cos(5.*self.k*self.dx)
        return (a,b,c,d)

    def PhaseError(self,arg):
        """
        Calculate phase error

        """
        self.k  = arg
        cf = self.__calcconst()
        cs = self.__calcsin()
        #--------------------------------------#
        # Calculate Symbolic Phase errors:     #
        #--------------------------------------#
        #- Second order
        omegaSecond = cf*cs[0]/self.dx
        omegaSecondTrad = cf*cs[0]/self.dx
        #- Fourth order
        omegaFourth = cf*(self.preS4a*cs[0]-self.preS4b*cs[2])/self.dx
        omegaFourthTrad = cf*(self.preC4a*cs[0]-self.preC4b*cs[1])/self.dx
        #- Sixth order
        omegaSixth = cf*(self.preS6a*cs[0]-self.preS6b*cs[2]+\
                     self.preS6c*cs[3])/self.dx
        omegaSixthTrad = cf*(self.preC6a*cs[0]-self.preC6b*cs[1]+\
                         self.preC6c*cs[2])/self.dx
        return (omegaSecond,omegaSecondTrad,omegaFourth,omegaFourthTrad,
                omegaSixth,omegaSixthTrad)

    def GroupVelocity(self,arg):
        """
        Compute Group velocity
        """

        self.k  = arg; cc = self.__calccos()
        #--------
        # Calculate Symbolic group velocities: 
        #--------
        #- Second order
        groupVSecond = self.c*cc[0]
        groupVSecondTrad = self.c*cc[0]
        #- Fourth order
        groupVFourth = self.c*(self.preS4a*cc[0]-self.preS4b*3.*cc[2])
        groupVFourthTrad = self.c*(self.preC4a*cc[0]-\
                           self.preC4b*2.*cc[1])
        #- Sixth order
        groupVSixth = self.c*(self.preS6a*cc[0]-self.preS6b*3.*cc[2]+\
                      self.preS6c*5.*cc[3])
        groupVSixthTrad = self.c*(self.preC6a*cc[0]-self.preC6b*2.*cc[1]+\
                          self.preC6c*3.*cc[2])
        return (groupVSecond,groupVSecondTrad,groupVFourth,
                groupVFourthTrad,groupVSixth,groupVSixthTrad)

    def aliasError(self,k1,k2):
        """
        alias error
        """

        self.k1=k1; self.k2=k2
        #--------------------------------------#
        # Calculate Interaction coeffients:    #
        #--------------------------------------#
        #- Second order
        aSec = 1./(2.*self.dx)*(np.sin((self.k1+self.k2)*self.dx)\
                             +np.sin(self.k2*self.dx)\
                             -np.sin(self.k1*self.dx))
        aSecTrad = 1./self.dx*(np.sin(self.k2*self.dx))
        #- Third order
        aThird = 1./self.dx*(4./3.*np.sin(self.k2*self.dx)\
                            -1./6.*np.sin(2.*self.k2*self.dx)\
                            -1j/3.*(1.-np.cos(self.k2*self.dx))**2)
        #- Fourth order
        aFou = 1./(384.*self.dx)*(243.*np.sin((self.k1+self.k2)*self.dx)\
                             +243.*np.sin(self.k2*self.dx)\
                             -27.*np.sin((2.*self.k1+self.k2)*self.dx)\
                             +27.*np.sin((self.k1-self.k2)*self.dx)\
                             +36.*np.sin(2.*self.k1*self.dx)\
                             -263.*np.sin(self.k1*self.dx)\
                             -9.*np.sin((2.*self.k1+3.*self.k2)*self.dx)\
                             -9.*np.sin((self.k1+3.*self.k2)*self.dx)\
                             +np.sin((3.*self.k1+3.*self.k2)*self.dx)\
                             +np.sin(3.*self.k2*self.dx)\
                             -np.sin(3.*self.k1*self.dx))
        aFouTrad = 1./self.dx*(4./3.*np.sin(self.k2*self.dx)\
                             -1./6.*np.sin(2.*self.k2*self.dx))
        #- Fifth order
        aFifth=1./self.dx*(3./2.*np.sin(self.k2*self.dx)\
                          -3./10.*np.sin(2.*self.k2*self.dx)\
                          +1./30.*np.sin(3.*self.k2*self.dx)\
                          +1j*2./15.*(1.-np.cos(self.k2*self.dx))**3)
        #- Sixth order
        aSix=150./(2.*2.*self.dx*128.*128.)*(300.*np.sin((self.k1+self.k2)\
                                   *self.dx)\
                                   +300.*np.sin(self.k2*self.dx)\
                                   -50.*np.sin((2.*self.k1+self.k2)*self.dx)\
                                   +50.*np.sin((self.k1-self.k2)*self.dx)\
                                   +6.*np.sin((3.*self.k1+self.k2)*self.dx)\
                                   -6.*np.sin((2.*self.k1-self.k2)*self.dx)\
                                   -350.*np.sin(self.k1*self.dx)\
                                   +56.*np.sin(2.*self.k1*self.dx)\
                                   -6.*np.sin(3.*self.k1*self.dx))\
                                   -25./(2.*6.*self.dx*128.*128.)\
                                   *(300.*np.sin((2.*self.k1+3.*self.k2)*self.dx)\
                                   +300.*np.sin((self.k1+3.*self.k2)*self.dx)\
                                   -50.*np.sin((3.*self.k1+3.*self.k2)*self.dx)\
                                   -50.*np.sin(3.*self.k2*self.dx)\
                                   +6.*np.sin((4.*self.k1+3.*self.k2)*self.dx)\
                                   -6.*np.sin((self.k1-3.*self.k2)*self.dx)\
                                   -300.*np.sin(2.*self.k1*self.dx)\
                                   -294.*np.sin(self.k1*self.dx)\
                                   +50.*np.sin(3.*self.k1*self.dx)\
                                   -6.*np.sin(4.*self.k1*self.dx))\
                                   +3./(2.*10.*self.dx*128.*128.)\
                                   *(300.*np.sin((3.*self.k1+5.*self.k2)*self.dx)\
                                   +300.*np.sin((2.*self.k1+5.*self.k2)*self.dx)\
                                   -50.*np.sin((4.*self.k1+5.*self.k2)*self.dx)\
                                   -50.*np.sin((self.k1+5.*self.k2)*self.dx)\
                                   +6.*np.sin((5.*self.k1+5.*self.k2)*self.dx)\
                                   +6.*np.sin(5.*self.k2*self.dx)\
                                   -300.*np.sin(3.*self.k1*self.dx)\
                                   -300.*np.sin(2.*self.k1*self.dx)\
                                   +50.*np.sin(4.*self.k1*self.dx)\
                                   +50.*np.sin(self.k1*self.dx)\
                                   -6.*np.sin(5.*self.k1*self.dx))
        aSixTrad=1./(30.*self.dx)*(45.*np.sin(self.k2*self.dx)\
                                 -9.*np.sin(2.*self.k2*self.dx)\
                                 +np.sin(3.*self.k2*self.dx))
        return (aSec,aSecTrad,aThird,aFou,aFouTrad,aFifth,aSix,aSixTrad)
        #return (abs(aSec),abs(aSecTrad),abs(aThird),abs(aFou),
        #        abs(aFouTrad),abs(aFifth),abs(aSix),abs(aSixTrad))


    #----- (C) ------------#
    #----- stability ------#
    #----- (C) ------------#
    #----- effectiveRes ---#

 def main():
     """
     main function
     """
     pass

    #------------------------------------------------------------------#

#*** Make the module executable as a script******#
if __name__ == "__main__":
    print ('Running as a script')


