# -*- coding: utf-8 -*-
"""
Created on Wed May 13 12:00:56 2015

@author: Pavel Esir
"""
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
from matplotlib.pylab import *
from numpy import *

SimTime = 20
C = 10.0
dNread = 20
freq = 8
N = 200
Mm = 0.5

Nreads = arange(dNread, N + 1, dNread, dtype='int')
Urange = arange(0.05, 1.0, 0.05)
Ierange = load('U_Iex_SimTime_20.0_h_0.0020_D_2.0_N_200_eps_0.010_m_{:.1f}.npy'.format(Mm))

name = 'U_{:.2f}_C_{:.1f}_N_{:d}_SimTime_{:n}.npy'

errR = np.load(name.format(Urange[0], C, N, SimTime))
arrShape = len(Urange), shape(errR)[0], shape(errR)[1], shape(errR)[2]
errR = zeros(arrShape)

for idx, u in enumerate(Urange):
    er = load(name.format(u, C, N, SimTime))
    errR[idx] = er
#%%

meanErrR = mean(errR, axis=3)
minLags = argmin(mean(errR, axis=3), axis=1)
minLagErr = amin(meanErrR, axis=1)

minErrNread = zeros(len(Nreads))

figure()
title('$C={}$'.format(C))
for idx, Nread in enumerate(Nreads):
    errs = minLagErr[:, int(Nread/dNread - 1)]*360/(4*np.pi)

    minErrNread[idx] = Urange[argmin(errs)]
    plot(Urange, errs, '-o', label=r'$N_{{read}}={}$'.format(Nread))

legend(loc='upper right', fontsize=14.0)
xlabel('U')
ylabel('$Error[degree \; ^{\circ}$]')
savefig("Freq_{}_C_{}_Tsim_{}_mean_{}_N_{}.png".format(freq, C, SimTime, Mm, N), dpi=300.)
#%%
figure()
title(r"$C = {}$".format(C))
plot(Nreads, minErrNread, '-o')
xlabel(r'$N_{read}$')
ylabel(r'U for minimal error')
savefig("U_for_min_Freq_{}_C_{}_Tsim_{}_mean_{}_N_{}.png".format(freq, C, SimTime, Mm, N), dpi=300.)
