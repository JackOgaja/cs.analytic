#!/usr/bin/env python
"""
------------------------------
This is a test script for the hosFunctions and hosGraphics module
by: Jack, 2015
------------------------------
"""

# >> modules <<<<
from hosFunctions import * 
from hosGraphics import * 
#import analyticFns 
#------------------

k  = np.linspace(0.00095,np.pi/205,5)
dx = 5.
c  = 5.
#task="dummy"

# Function definition is here
def printme( str ):
   "This prints a passed string into this function"
   print str
   return;

# Now you can call printme function
printme("I'm first call to user defined function!")
printme("Again second call to the same function")

numContent=dir(np)
loc=locals() #--Returned in the form of a dictionary!--
ky=loc.keys()
print ""
print "=================="
#print numContent
print "=================="
print ""
print ".............."
#print loc
print ""
#print ky
#-------
#test=CalcAnalytic("phase") 
#t   =CalcAnalytic("notNow")
#CalcAnalytic.task="phase"
test     = CalcAnalytic() # Analytic functions object 
t        = CalcAnalytic()
t.new    = np.pi # dynamic attribute creation!
test.task="dummy"
test.k   = np.linspace(0.0010,np.pi/505,5)
test.c   = 500.
test.dx  = 9999.
d        = test.__dict__
print ""
print "=================="
print "----MODULE TEST---" 
print "=================="
print " test = ", test
print " testOmegaSecond = ", test.omegaSecond(test.k, dx, c)
print " testOmegaFourth = ", test.omegaFourth(k, test.dx, test.c)
print " testOmegaSixth  = ", test.omegaSixth(k, dx, c)
print ".............."
print ""
print " NEW PIE = ", t.new
print " DICTIONARY = ", d
print "=================="
print ""

plotSave=False
addRefline=False
daten=test.omegaSecond(test.k, dx, c) 
#+p = hosPlot(plotSave) # Analytic functions object
#p.legs="anAriekKoa90"
#+p.legs=p.legs.append("anAriekKoa90")
#+plt=p.plotPhase(daten,addRefline)

print ""
print "=================="
print " ANALYTICPLOT = ", p
#print " DICTIONARY = ", d
print "=================="
print ""


