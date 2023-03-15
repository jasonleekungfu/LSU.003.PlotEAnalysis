# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 14:45:36 2020

@author: jason
"""

import numpy as np
import re
from scipy.special import lpn
import scipy.interpolate as interpolate
from scipy.optimize import curve_fit


#-----------------------------------------------------------------------
#                           Constants
#-----------------------------------------------------------------------
    
Pi = 3.141592653589793238462643             # Pi
eV = 1/27.211396                            # eV
nm = 1/0.05291772192                        # nm   
Ang = 1/0.5291772192                        # Angstrom
c  = 137.035999173                          # Speed of light
me = 1                                      # Electron mass
asec = 1/24.18884326505                     # as
fsec = 1000/24.18884326505                  # fs
Wpcm2 = 1/(3.51e16)                         # W/cm^2
deg = Pi/180                                # Degree
small = 1e-7                                # Infinitasimal


#-----------------------------------------------------------------------
#                           Functions
#-----------------------------------------------------------------------

lnum = { "s": 0, "p": 1, "d": 2, "f": 3, "g": 4}
llet = {v: k for k, v in lnum.items()}

def read2DData(fp, emit=0):
    """
    Read 2D data file.
    """
    data = []
    with open(fp, "r") as f:
        for i in (range(emit)):   # Emit given number of lines
            f.readline()
        for line in f:
            dataline = list(map(float, line.split()))
            data.append(dataline)
    return(np.array(data).transpose())

def readPopulation(intensities, folder, n=2, l="p"):
    """
    Read, proccess, and export population data.
    """
    populations = np.zeros(intensities.size, dtype=float)
    with open(folder+"populations."+str(n)+l+".txt", "w") as f:
    
        # Read and process
        for i in range(intensities.size):
    
            # Read populations
            data = read2DData(folder + str(intensities[i]) + \
                              "/populations."+l+".out")[:,-1]
            populations[i] = data[n-lnum[l]]
                
        # Export processed data
        for i in range(intensities.size):
            f.write(str(intensities[i]) + "\t" + str(populations[i]) + "\n")

def readPhase(intensities, folder, n=2, l="p"):
    """
    Read, proccess, and export phase data.
    """

    lnum = { "s": 0, "p": 1, "d": 2, "f": 3 }
    phases = np.zeros(intensities.size, dtype=float)
    with open(folder+"phases."+str(n)+l+".txt", "w") as f:
    
        # Read and process
        for i in range(intensities.size):
    
            # Read phases and connect 2Pi jumps
            data = read2DData(folder + str(intensities[i]) + \
                              "/phases."+l+".out")[:,-1]
            phases[i] += data[n-lnum[l]]
            if (i>2 and (phases[i]-phases[i-1])*(phases[i-1]-phases[i-2])<0 \
                    and (phases[i-1]-phases[i-2])*(phases[i-2]-phases[i-3])<0 \
                    and abs(phases[i-1]-phases[i-2]) > np.pi \
                or i==2 and (phases[2]-phases[1])*(phases[1]-phases[0])<0 \
                        and abs(phases[1]-phases[0]) > np.pi):
                if (phases[i]-phases[i-1]>0):
                    phases[i-1:] += 2*np.pi
                elif (phases[i]-phases[i-1]<0):
                    phases[i-1:] -= 2*np.pi
                
        # Export processed data
        phases -= phases[0]
        for i in range(intensities.size):
            f.write("{:f}\t{:f}\n".format(intensities[i], -phases[i]))

def readAmpLost(intensities, folder):
    """
    Read, proccess, and export amplitude lost data.
    """

    amplost = np.zeros(intensities.size, dtype=float)
    p = re.compile('total lost amplitude')
    with open(folder+"amplost.txt", "w") as f:
    
        # Read and process
        for i in range(intensities.size):
            
            # Open file
            with open(folder + str(intensities[i]) + "/TDSE.o") as f2:
                
                # Read all lines
                for line in f2:
                    
                    # Proceed if found keyword "total lost amplitude"
                    if (p.findall(line)):
                        amplost[i] = float(line.split()[-1])
                
        # Export processed data
        for i in range(intensities.size):
            f.write("{:f}\t{:f}\n".format(intensities[i], amplost[i]))
            
def print1DData(ary1, ary2, filename):
    with open(filename, 'w') as f:
        #np.savetxt(f, ary1, ary2)
        for i in range(ary1.size):
            f.write(str(ary1[i])+'\t'+str(ary2[i])+'\n')
            #f.write("{:f}\t{:.8f}\n".format(ary1[i], ary2[i]))

def getYield(filename, angle=Pi):
    """
    Parameters
    ----------
    filename : string
        Data file path and name.
    angle : float (0, Pi), optional
        Collection angle. The default is Pi (integrate over all angles).

    Returns
    -------
    None.
    """
    
    # Integral from opening angle to Pi
    x = np.cos(angle)
    yieldAry = read2DData(filename,1)
    res = np.zeros(yieldAry[0].size)
    legendres = lpn(yieldAry[:,0].size-1,x)
    for i in range(yieldAry[:,0].size-1):
        if (i==0):
            res = res + yieldAry[1]*(1-x)
        else:
            res = res + yieldAry[i+1]*(1-x**2)/i/(i+1)*legendres[1][i]
    
    # Integral over all angles
    return(res)

def getP(yieldAry, EAry, order1, order2=-1, delta=0.4, ACShift=0):
    """
    Parameters
    ----------
    yieldAry : array (float)
        Array of yield.
    EAry : array (float)
        Array of energies (in the units of harmonic orders).
    order1, order2 : integer
        Particular orders of two sidebands (ascendant) to retrieve the yield.
        If order2 is not given (default -1), then yield of order 1 is returned.
    delta : float, optional
        Integrating limit, plus/minus delta. The default is 0.4.
    ACShift : float, optional
        AC Stark shift correction amount. The default is 0.

    Returns
    -------
    P parameter between two given orders.
    """
    
    tck = interpolate.splrep(EAry, yieldAry)
    if (order2 > 0):
        S1 = interpolate.splint(order1+ACShift-delta, order1+ACShift+delta, tck)
        S2 = interpolate.splint(order2+ACShift-delta, order2+ACShift+delta, tck)
        return((S2-S1)/(S2+S1))
    else:
        return(interpolate.splint(order1+ACShift-delta, order1+ACShift+delta, tck))

def getRho(P1, P2):
    
    return((np.average(P1*P2)-np.average(P1)*np.average(P2)) / \
           np.sqrt(np.average(P1**2)-np.average(P1)**2) / \
           np.sqrt(np.average(P2**2)-np.average(P2)**2))

def mycos(x, omega, dphi, A, B, C):
    return(A*np.cos(omega*x+dphi) + B*x + C)

def fitcos(x, y, params=[6*Pi, 0, 1, 0, 0]):
    params = curve_fit(mycos, x, y, p0=params)[0]
    if (params[2] < 0):
        params[2] = -params[2]
        params[1] = params[1] + Pi
    if (params[1] < 0):
        params[1] = params[1] + 2*Pi
    return(curve_fit(mycos, x, y, p0=params)[0])

def getPhaseDif(x,y):
    dphi = x-y-Pi
    while abs(dphi)>Pi:
        if (dphi>0):
            dphi = dphi - 2*Pi
        else:
            dphi = dphi + 2*Pi
    return(dphi)
        
    
