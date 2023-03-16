from lib import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.special import lpn


#-----------------------------------------------------------------------
#                        User-defined ariables
#-----------------------------------------------------------------------

atom = "he"         # "he", "ne"
scale = 'log'       # 'log' or 'linear'
ndelay = 95         # Max delay index
lambdaIR = 800      # Integer [nm]
iPath = "/work/jasonli3/Program/2.EAnalysis/Output"         # Input path
oPath = "/work/jasonli3/Program/3.PlotEAnalysis/Output"     # Output path
subfolder = 'test'  # Subfolder name

#-----------------------------------------------------------------------
#                            Derived variables
#-----------------------------------------------------------------------

lambdaIR = lambdaIR*nm
TIR = lambdaIR/c
omegaIR = c/lambdaIR*2*Pi
filename = iPath + "/" + subfolder + "/B/"

delayAry = np.zeros(ndelay+1, dtype=float) 
EAry = (read2DData(filename+"0",1)[0]/(omegaIR/eV)+Ip[atom]/omegaIR)/3
    # Instead of eV, use Harmonics
nE = EAry.size
yieldAry = np.zeros([ndelay+1,nE], dtype=float)
      
for i in range(ndelay+1):
    yieldAry[i] = getYield(filename+str(i), angle=Pi/2)
    with open(filename+str(i), "r") as f:
        delayAry[i] = float(list(map(float, f.readline().split()))[-1])


#-----------------------------------------------------------------------
#                               Plot
#-----------------------------------------------------------------------

fig, ax = plt.subplots(figsize=[7,10],\
            gridspec_kw = { 'left':         0.1,\
                            'right':        0.9,\
                            'top':          0.9,\
                            'bottom':       0.15})  

ax.set_title('')    
ax.set_xlabel('XUV(APT)-IR time delay [$T_{IR}$]')
ax.set_ylabel('Harmonic order')
ax.set_yticks([6,7,8,9,10,11,12])
if (scale == 'log'):
    ax.pcolormesh(delayAry, EAry, yieldAry.transpose(), edgecolors='face',\
                  norm=colors.LogNorm(vmin=yieldAry.max()*1e-2, \
                                      vmax=yieldAry.max()))
else:
    ax.pcolormesh(delayAry, EAry, yieldAry.transpose(), edgecolors='face',\
                  vmin=yieldAry.min(), vmax=yieldAry.max())
    
plt.savefig(oPath + "/"+subfolder+".png", dpi=100)  

