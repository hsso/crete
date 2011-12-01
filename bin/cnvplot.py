#! /usr/bin/python
#
# Plot amc convergence history.
# Argument 1: name of input file (i.e., amc output file)
# Argument 2: transition to be plotted (in order listed in molfile)
#
# Floris van der Tak & Michiel Hogerheijde
# first version: 15 Dec 2004
# this  version: 22 Feb 2005
#
# This script uses the ppgplot library -- the Python version of pgplot.
# See http://efault.net/npat/hacks/ppgplot/ for details & downloads.
#
import os
import sys
import math
from ppgplot import *
from Numeric import array,ones,log
from numarray import where

# command line input: generic file name & line to be plotted
if (len(sys.argv) < 3):
    print "Usage: "+sys.argv[0]+" [file name] [line number]"
    sys.exit()
fname    = str(sys.argv[-2])
plotline = int(sys.argv[-1])
llabel   = 'File: '+fname+'   Line: '+str(sys.argv[-1])

# Minimum population which amc checks for convergence.
# Default is 1e-6 ; change if necessary 
minpop = 1.0e-8
tiny   = 1.0e-20

# fundamentals
hp = 6.6256e-27
kb = 1.3805e-16

# find the input files
flist = []
for file in os.listdir('.'):
    if (file[0:len(fname)] == fname) and (file[-3:] == 'his'):
        flist.append(file)
# make sure that real output file comes last
for file in os.listdir('.'):
    if (file == fname):
        flist.append(file)

# minimum 1, maximum 10 files (for clarity; only got 10 colours)
nfl = len(flist)
red = nfl / 10 + 1
if (nfl < 1):
    print "Error: No files found"

# get line freqs & stat weights from molecular data file
file = open(flist[0],'r')
line = ''
while (line[0:7] != 'molfile'):
    line = file.readline()
file.close()
molfile = line[8:-1]

file = open(molfile,'r')
for idum in range(5):
    file.readline()
nlev = int(file.readline())
file.readline()
gstat = ones(nlev,'d')

for ilev in range(nlev):
    line  = file.readline()
    words = line.split()
    gstat[ilev] = float(words[2])

file.readline()
nlin = int(file.readline())
file.readline()
frq = ones(nlin,'d')
gup = ones(nlin,'d')
glo = ones(nlin,'d')
jup = ones(nlin,'i')
jlo = ones(nlin,'i')

for iline in range(nlin):
    line  = file.readline()
    words = line.split()
    jup[iline] = int(words[1])
    jlo[iline] = int(words[2])
    gup[iline] = gstat[jup[iline]-1]
    glo[iline] = gstat[jlo[iline]-1]
    frq[iline] = float(words[4])*1.0e9

def exct(iline,upper,lower):
    above = -1.0 * hp * frq[iline] / kb
    below = upper / lower * glo[iline] / gup[iline]
    return above / log(below)

# load populations & calculate T_ex
file = open(flist[0],'r')
line = ''
while (line[0:5] != 'ncell'):
    line = file.readline()
ncell = int(line[6:-1])
txc   = ones((ncell,nfl),'d')

ifile = 0
flag  = 0 
for ff in flist:
    file = open(ff,'r')
    while (line != '@\n'):
        line = file.readline()
    for icell in range(ncell):
        line  = file.readline()
        words = line.split()
        popup = float(words[9+jup[plotline]-1])
        poplo = float(words[9+jlo[plotline]-1])
        #set T_ex=-999.99 if pops are untested
        if ((popup < minpop) or (poplo < minpop)):
            txc[icell,ifile] = -999.99
            # but don't flag for vacuum case
            if ((popup > tiny) and (poplo > tiny)):
                flag = 1
        else:
            txc[icell,ifile] = exct(plotline,popup,poplo)
    file.close()
    ifile += 1
    
if (flag == 1):
    print
    print "Warning: Some populations were below limit of ",minpop
    print

# Representation 1: x-y plot of Txc vs radius; colour code time curves

# plot limits -- take care to exclude 'blanking' value of -999 
x0 = 0
x1 = ncell*1.2
y0 = -999.99
icell = -1
while (y0 == -999.99):
    icell += 1
    y0 = txc[icell,0]
for icell in range(ncell):
    if ((txc[icell,0] > -999) and (txc[icell,0] < y0)):
        y0 = txc[icell,0]
y1 = max(txc[:,0])

for ifl in range(nfl):
    for icell in range(ncell):
        if ((txc[icell,0] > -999) and (txc[icell,0] < y0)):
            y0 = txc[icell,0]
    y1 = max(y1,max(txc[:,ifl]))
y0 /= 1.2
y1 *= 1.2

print "min/max T_ex:",y0,y1

# draw box 
pgbeg('?')
pgslw(4)
pgsch(1.4)
pgscf(2)
pgenv(x0,x1,y0,y1,0,0)
pglab('Cell Number','Excitation Temperature (K)',llabel)
pgsch(1.2)

# draw curves (only one transition for now)
xpnts = range(ncell)
icc   = 0
for ifl in range(1,nfl,red):
    ypnts = txc[:,ifl]
    icc  += 1
    pgsci(icc)
    pgline(array(xpnts),array(ypnts))
    ltext = 'Iteration '+str(ifl)
    pgtext(ncell-0.7,y1-10.0-ifl/red/10.0*(y1-y0),ltext)

pgsci(1)
pgpage()

# Representation 2: contour plot 
nclev = 9
clevs = ones(nclev,"d")
# linear contours if dynamic range < 10, otherwise logarithmic ones.
if (y1/y0 < 10):
    clevs = array(range(nclev))*(y1-y0)/nclev + y0
    print "Using linear contour spacing"
else:
    print "Using logarithmic contour spacing"
    step = math.log10(y1/y0)/nclev
    for ilevs in range(1,nclev+1):
        clevs[ilevs-1] = y0 * 10**(step*ilevs)
print "Contour levels: ",clevs

#space and time axes
pgenv(0,nfl*1.2,x0,x1)
pglab('Iteration Number','Cell Number',llabel)
pgcont_s(txc,nclev,clevs)
pgconl_s(txc,clevs[1],str(round(clevs[1],2)),10,5)
pgconl_s(txc,clevs[4],str(round(clevs[4],2)),10,5)
pgconl_s(txc,clevs[7],str(round(clevs[7],2)),10,5)

pgpage()
pgend
