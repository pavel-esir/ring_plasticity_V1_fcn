# -*- coding: utf-8 -*-
"""
Created on Wed May 13 12:00:56 2015

@author: Pavel Esir
"""
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
from matplotlib.pylab import *
from numpy import *

SimTime = 2000
C = 40.0
dNread = 20
freq = 4
T = 0.05
N = 200
Mm = 0.5

Nreads = arange(dNread, N + 1, dNread, dtype='int')
Urange = arange(0.05, 0.951, 0.05)
#Ierange = load('U_Iex_SimTime_20.0_h_0.0020_D_2.0_N_200_eps_0.010_m_{:.1f}.npy'.format(Mm))

#folderName = 'res_h_0.0020_D_2.0_freq_{:.1f}_T_{:.2f}/'.format(freq, T)
folderName = 'res_h_0.0020_D_2.0_freq_{:.1f}_T_{:.2f}_m_{:.1f}/'.format(freq, T, Mm)
name = folderName + 'U_{:.2f}_C_{:.1f}_N_{:d}_SimTime_{:d}.npy'
nameExactR = folderName + 'U_{:.2f}_C_{:.1f}_N_{:d}_SimTime_{:d}_exactPV.npy'

errR = np.load(name.format(Urange[0], C, N, SimTime))
arrShape = len(Urange), shape(errR)[0], shape(errR)[1], shape(errR)[2]
errR = zeros(arrShape)

arrShape = len(Urange), shape(errR)[1], shape(errR)[3]
errExactR = zeros(arrShape)

for idx, u in enumerate(Urange):
    er = load(name.format(u, C, N, SimTime))
    errR[idx] = er
    er = load(nameExactR.format(u, C, N, SimTime))
    errExactR[idx] = er  
#%%

meanErrR = mean(errR, axis=3)
minLags = argmin(meanErrR, axis=1)
minLagErr = amin(meanErrR, axis=1)

meanErrExactR = mean(errExactR, axis=2)
minLagsExact = argmin(meanErrExactR, axis=1)
minLagErrExact = amin(meanErrExactR, axis=1)


minErrNread = zeros(len(Nreads))
NreadsForPlot = [20, 40, 80, 200]
for idx, Nread in enumerate(Nreads):
    errs = minLagErr[:, int(Nread/dNread - 1)]*360/(4*np.pi)

    minErrNread[idx] = Urange[argmin(errs)]
    if Nread in NreadsForPlot:
        figure(1, figsize=(4*2.5, 3*2.5))
        plot(Urange, errs, '-o', label=r'$N_{{read}}={}$'.format(Nread))
        
        figure(2, figsize=(4*2.5, 3*2.5))
        plot(Urange, minLags[:, int(Nread/dNread - 1)]*2., '-o', label=r'$N_{{read}}={}$'.format(Nread))

figure(1, figsize=(4*2.5, 3*2.5))
plot(Urange, minLagErrExact*360/(4*pi), '--o', label='exact readout')
title('$C={}$'.format(C))
legend(loc='upper right', fontsize=14.0)
xlabel('U')
ylabel('$Error[degree \; ^{\circ}$]')
savefig("Freq_{}_C_{}_Tsim_{}_mean_{}_N_{}.png".format(freq, C, SimTime, Mm, N), dpi=260.)

figure(2, figsize=(4*2.5, 3*2.5))
plot(Urange, minLagsExact*2, '--o', label='exact readout')
title('$C={}$'.format(C))
legend(loc='upper right', fontsize=14.0)
xlabel('U')
ylabel('$lag[ms]$')
savefig("lag_Freq_{}_C_{}_Tsim_{}_mean_{}_N_{}.png".format(freq, C, SimTime, Mm, N), dpi=260.)

#%%
figure(3, figsize=(4*2.5, 3*2.5))
#title(r"$C = {}$".format(C))
plot(Nreads, minErrNread, '-o', label='$C= {}$'.format(C))
xlabel(r'$N_{read}$')
ylabel(r'U for minimal error')
savefig("U_for_min_Freq_{}_C_{}_Tsim_{}_mean_{}_N_{}.png".format(freq, C, SimTime, Mm, N), dpi=260.)
