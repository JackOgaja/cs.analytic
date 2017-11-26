# -*- coding: utf-8 -*-
#
"""
module for plotting analytical analysis
of Higher Order Schemes for Geophysical models
"""
# module:    hos_graphics
# created:   Jack Ogaja <jack.ogaja'at'gmail.com>
#
# purpose:   module for analytic functions of
#            Conservative HOS
#
#-----------------------------------------------------------------------------#

import matplotlib.pyplot as pl                   # more control than pylab
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator, MaxNLocator
import os, sys, inspect
#-----------------------------------------------------------------------------#

class hosPlot:
    """
    Attributes ->
              labels:    legend labels
              title:     plot title
              styles:    plot styles
              lwidth:    plot line width
              legLoc:    legend location
              savePlot:  option to save a plot to a file/dir
              plotName:  plot name
              plotDir:   plot directory/folder
              outFormat: format of the output file |pdf/eps/png|
              lblsize:   label fornt size
              ttlsize:   title font size
              ymin:      minimum y-axis limit
              ymax:      maximum y-axis limit
              xmin:      minimum x-axis limit
              xmax:      maximum x-axis limit
              ylabel:    y-axis label
              xlabel:    x-axis label
              lgframe:   legend box 'on' or 'off'
              lgFontSize:legend font size
    """

    def __init__(self):
        self.labels = []
        self.title  = ''
        self.styles = []
        self.lwidth = []
        self.legLoc = 'best'
        self.lblsize = 10
        self.ttlweight ='normal'
        self.ttlsize = 10
        self.savePlot =None
        self.minorticks =None
        self.plotName ='plot'
        self.plotDir  ='out'
        self.outFormat ='png'
        self.ymin =None
        self.ymax =None
        self.xmin =None
        self.xmax =None
        self.ylabel ='y'
        self.xlabel ='x'
        self.lgframe =None
        self.lgFontSize =10

    def __plotSettings(self):
        rc('text',usetex=True)                      # LaTeX support
        rc('font',family='serif',serif='cmr10')     # LaTeX Computer Modern font
                                                    # (save: raw strings r'...')
        rc('legend',**{'fontsize':self.lgFontSize}) # global: legend's fontsize

        if self.minorticks: pl.minorticks_on()
        if self.ymin: pl.ylim(ymin=self.ymin)
        if self.ymax: pl.ylim(ymax=self.ymax)
        if self.xmin: pl.xlim(xmin=self.xmin)
        if self.xmax: pl.xlim(xmax=self.xmax)

    def __consistency(self):
        if len(self.styles) > len(self.y):
            errmsg1 = 'There are more plot styles than the plots!'
            self.__printError(errmsg1)
            self.__standarderror()
            return 1
        elif len(self.styles) < len(self.y):
            warningmsg1 = 'There are less plot styles than the plots...'
            warningmsg2 = 'Plotting will continue with default settings...'
            self.__printWarning(warningmsg1,warningmsg2)
            diff = len(self.y)-len(self.styles)
            for d in range(diff):
                self.styles.append('')
        if len(self.labels) > len(self.y):
            errmsg1 = 'There are more plot labels than the plots!'
            self.__printError(errmsg1)
            self.__standarderror()
            return 1
        elif len(self.labels) < len(self.y):
            warningmsg1 = 'There are less plot labels than the plots...'
            warningmsg2 = 'Plotting will continue without labels...'
            self.__printWarning(warningmsg1,warningmsg2)
            diff = len(self.y)-len(self.labels)
            for d in range(diff):
                self.labels.append('plot'+str(len(self.labels)+1))

    def __createDir(self):
        scriptDir = os.getcwd()
        pltDir    = os.path.join(scriptDir, self.plotDir)
        if not os.path.isdir(pltDir):
            os.makedirs(pltDir)
        return pltDir

    #---- Line plot method ------#
    def plotLine(self, k, *daten, **kwargs):
        "Plot multiple curves"
        self.y=daten

        #--- Update the settings -------#
        self.__plotSettings()
        (refLine,hLine,vLine)=(None,None,None)
        lw = []

        #fig1 = pl.figure(1)
        #ax11 = fig1.add_subplot(111)
        #ax11 = fig1.subplots(111)
        fig, ax11 = pl.subplots()

        #---- Add titles ----#
        ax11.set_title(self.title, fontsize=self.ttlsize,\
                        fontweight=self.ttlweight)
        ax11.set_xlabel(self.xlabel, fontsize=self.lblsize)
        ax11.set_ylabel(self.ylabel, fontsize=self.lblsize)

        #---- Evaluations ---#
        cons=self.__consistency()
        if cons: sys.exit(cons)

        if kwargs:
            for key, value in kwargs.items():
                if key == 'refLine':
                    refLine=True
                    if value.__class__ == tuple:
                    #if hasattr(value, '__iter__'):
                        refs=self.__evaluatereferenceLines(key,value)
                        (refValue, reflabel)=refs
                    else:
                        (refValue,reflabel)=(value, None)
                        print('Key is:{0}, Value is:{1}'.format(key,value))
                elif key == 'hLine':
                    hLine=True
                    if value.__class__ == tuple:
                        refs=self.__evaluatereferenceLines(key,value)
                        print('REFS = {0}'.format(refs))
                        (hValue, hlabel)=refs
                    else:
                        (hValue,hlabel)=(value, None)
                        print('Key is:{0}, Value is:{1}'.format(key,value))
                elif key == 'vLine':
                    vLine=True
                    if value.__class__ == tuple:
                        refs=self._evaluatereferenceLines(key,value)
                        (vValue,vlabel)=refs
                    else:
                        (vValue,vlabel)=(value, None)
                        print('Key is:{0}, Value is:{1}'.format(key,value))
                else:
                    warningmsg=('{0} is not an attribute '
                                'for the plot instance'.format(key))
                    self._printWarning(warningmsg)

        for i in range(len(self.y)):
            if len(self.lwidth)==0: lw.append(1)
            else: lw.append(self.lwidth[i])
            ax11.plot(k, self.y[i], self.styles[i],linewidth=lw[i],
                      label=self.labels[i])

        if refLine:
            ax11.plot(k,refValue,ls=':',linewidth=2,
                      color='k', label=reflabel)
        if hLine:
            ax11.axhline(hValue, linewidth=2, ls=':',
                         color='k', label=hlabel)
        if vLine:
            ax11.axvline(vValue, linewidth=2, ls=':',
                         color='k', label=vlabel)

        ax11.legend(loc=self.legLoc,ncol=1,handlelength=2.5,
                    frameon=self.lgframe)

        pl.draw()

        #---- Save the plot ---#
        pdir=self.__createDir()
        fileN = self.plotName+'.'+self.outFormat
        if self.savePlot:
            fig1.savefig(os.path.join(pdir, fileN), dpi=200)
        else:
            pl.show()

        pl.close(1)

    #---- Contour plot method ------#
    def plotContour(self, k1, k2, daten, **kwargs):
        "Contour plot"
        pl.figure()
        CS = pl.contour(k1, k2, daten,
                        colors='k',  # negative contours dashed by default
                       )
        pl.clabel(CS, fontsize=9, inline=1)
        pl.title('Single color - negative contours dashed')
        pl.ylim(-1.5, 0.)
        pl.show()

    def testContour(self):
        " A test for contour plots "
        import numpy as np
        X = np.linspace(0, 1, 10)
        Y = np.linspace(0, 1, 10)
        x,y = np.meshgrid(X,Y)
        f1 = np.cos(x*y)
        f2 = x-y
        pl.contour(x,y,f2,colors='red')
        pl.contour(x,y,f1,colors='blue')
        pl.show()

    def __evaluatereferenceLines(self,lkey,lvalue):
        err1=0; pdata=lvalue[0]
        #exec('try:\n\tint(lvalue[1])\nexcept ValueError:\n\terr1=42',{"err1:42})
        try: int(lvalue[1])
        except ValueError: err1=42
        if err1==42:
            if len(lvalue) > 2:
                errmsg=('Error occured in values of \'{0}\''
                        .format(lkey))
                self.__printError(errmsg)
                self.__standarderror()
                sys.exit(err1)
            else:
                return(pdata, lvalue[1])

    def __standarderror(self):
        errst=inspect.stack()[1][3]
        print >> sys.stderr, ('Script quits in \'{0}\' function '
                              'in \'{1}\' module\n'
                              .format(errst,__name__))

    def __printWarning(self,*args):
        print('*'*10)
        for arg in args: print('WARNING: '+arg)
        print('*'*10)

    def __printError(self,*args):
        print('*'*10)
        for arg in args: print('ERROR: '+arg)
        print('STOP : The script stops here')
        print('*'*10)

#-----------------------------------------------------------------------------#
#*** Make the module executable as a script******#
if __name__ == "__main__":
    fname=__file__.split('/')
    print('\nThis is a test run for the \'{0}\' module.\n'
          'The test run is the simplest use of the module\n'
          .format(fname[-1]))
    import numpy as np
    x = np.linspace(0, 2 * np.pi, 50, endpoint=True)
    y1 = 3 * np.sin(x)
    y2 = np.sin(2*x)
    y3 = 0.3 * np.sin(x)
    y4 = np.cos(x)
    hp=hosPlot()
    hp.styles=['b-','r--','g:','m-.']
    hp.lwidth=[2.5,1.5,2,2]
    hp.plotLine(x,y1,y2,y3,y4,hLine=(1,'hLine'))
