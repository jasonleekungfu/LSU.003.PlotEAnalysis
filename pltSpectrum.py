from lib import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.special import lpn


#-----------------------------------------------------------------------
#                               Plot
#-----------------------------------------------------------------------

scale = 'log'       # 'log' or 'linear'
ndelay = 95         # Max delay index
IIR = '3.0'         # In W/cm^2
lambdaIR = 800      # In nm
iPath = "/work/jasonli3/Program/2.EAnalysis/Output"         # Input path
oPath = "/work/jasonli3/Program/3.PlotEAnalysis/Output"     # Output path
subfolder = 'test'  # Subfolder name
filename = iPath + "/" + subfolder + "/B/"


lambdaIR = lambdaIR*nm
TIR = lambdaIR/c
omegaIR = c/lambdaIR*2*Pi

delayAry = np.linspace(-1,1,ndelay+1) #* (TIR/fsec)
#EAry = read2DData(filename+"0",1)[0]/(omegaIR/eV)+0.793253117866589/omegaIR #Ne
EAry = read2DData(filename+"0",1)[0]/(omegaIR/eV)+0.903247020452234/omegaIR #He
nE = EAry.size
yieldAry = np.zeros([ndelay+1,nE], dtype=float)
      
for i in range(ndelay+1):
    yieldAry[i] = getYield(filename+str(i), angle=Pi/2)

#-----------------------------------------------------------------------

fig, ax = plt.subplots(figsize=[7,10],\
            gridspec_kw = { 'left':         0.1,\
                            'right':        0.9,\
                            'top':          0.9,\
                            'bottom':       0.15})  

ax.set_title('IR intensity: '+IIR+' W/cm$^2$')    
ax.set_xlabel('XUV(APT)-IR time delay [$T_{IR}$]')
ax.set_ylabel('Harmonic order')
ax.set_yticks([6,7,8,9,10,11,12])  # Ne
#ax.set_yticks([6,7,8,9,10,11])  #He
ax.set_ylim(5,12)
if (scale == 'log'):
    print(yieldAry.min(),yieldAry.max())
    ax.pcolormesh(delayAry, EAry/3, yieldAry.transpose(), edgecolors='face',\
#                  norm=colors.LogNorm(vmin=4e-5, \
#                                      vmax=3e-3))
#                  norm=colors.LogNorm(vmin=yieldAry.min(), \
                  norm=colors.LogNorm(vmin=yieldAry.max()*1e-2, \
                                      vmax=yieldAry.max()))
else:
    print(yieldAry.min(),yieldAry.max())
    ax.pcolormesh(delayAry, EAry/3, yieldAry.transpose(), edgecolors='face',\
                  vmin=yieldAry.min(), vmax=yieldAry.max())
ax2 = ax.twinx()
#ax2.set_ylim(5,24) # Ne
ax2.set_ylim(0,42) # He  0,40
ax2.set_ylabel('Photoelectron energy [eV]')
    
plt.savefig(oPath + "/"+subfolder+".png", dpi=100)  

# Print data
with open(oPath + '/data.txt', "w") as f:
    for i in range(0,yieldAry[0].size,8):
        f.write(str(EAry[i])+'\t'+str(yieldAry[25,i])+'\n')
