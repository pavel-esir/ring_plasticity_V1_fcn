# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 11:55:48 2015

@author: Pavel Esir
"""
from __future__ import division
from matplotlib.gridspec import GridSpec
from matplotlib.pylab import *
from numpy import *
seed(0)

# if CalcMode is equal to 0 then a single realisation is calculated and then printed

# if CalcMode is equal to 1 then realisation for different U is calculated and saved into files

# if CalcMode is equal to 2 then U is set from argv and then single realisation is calculated
# can be used for parallel calculation with GNU parallel

# if CalcMode is equal to 3 then I0 vs U curve is calculated by bisection method
CalcMode = 0


SimTime = 20.0        # seconds
h = 0.002             # seconds
pltSampl = 0.02       # variable save interval in s
Tsim = int(SimTime/h)

# to load empirical calculated U vs I0 dependence
Mm = 0.5
Urange = arange(0.05, 1.0, 0.05)

Ierange = [-1.170, -0.917, -0.669, -0.485, -0.485, -0.542, -0.711,
           -0.934, -1.154, -1.398, -1.690, -1.950, -2.256, -2.389,
           -2.543, -2.717, -2.855, -3.033, -3.150]
#Ierange = load('U_Iex_SimTime_20.0_h_0.0020_D_2.0_N_200_eps_0.010_m_{:.1f}.npy'.format(Mm))
#%%
U = 0.25
I0 = Ierange[int(U/0.05) - 1]
#U = 0.0
#I0 = Ierange[0]

D = 2.0
J0 = -12
J1 = 30

# duration of events in sec
T = 1.00
# amplitude of events
#C = 20.0
C = 10.0
# input poisson events rate, Hz
freq = 100

tau_r = 0.01

tau = 0.01
tau_n = 0.1
tau_rec = 0.8
N = 200

folderName = 'res_h_0.0020_D_2.0_freq_{:.1f}_T_{:.2f}/'.format(freq, T)

stime = arange(0, SimTime, h)
m = zeros(N)

ActM = zeros((int(SimTime/pltSampl), N), dtype='float32')
ActX = zeros((int(SimTime/pltSampl), N), dtype='float32')
x = ones(N)

Inoise = zeros(N)
Iex = zeros(Tsim)
ThetaEs = zeros(Tsim)

th = linspace(-pi, pi, N, endpoint=False)
alpha, beta = meshgrid(th, th)
W = (J0 + J1*cos(alpha - beta))/N
del alpha, beta
#%%
# number of input events, nearly equals poisson rate*time
Nev = int(freq*SimTime)

# times of Poisson events have exponential distribution
inpTimes = (exponential(1/freq, Nev) + 2*T).cumsum()/h
inpTimes = array(inpTimes[inpTimes < Tsim - 2*T/h], dtype='int')

# cutoff events which occur later than simulation time
Nev = len(inpTimes)

# uniformly distributed angles of input events
inpTheta = uniform(-pi, pi, Nev)

for t, theta in zip(inpTimes, inpTheta):
    Iex[t:int(t + T/h)] = C
    ThetaEs[t:int(t + T/h)] = theta

Nreads = linspace(int(N/10), N, 10, dtype='int')
R = zeros((Tsim, len(Nreads)), dtype='complex')

exactR = zeros(Tsim, dtype='complex')

choises = [[]]*len(Nreads)
for idx, n in enumerate(Nreads):
    choises[idx] = choice(range(N), n, replace=None)
#%%
def integrate():
    global Inoise, x, m
    for t in xrange(0, Tsim - 1):
        Inoise = Inoise + (-Inoise*h/tau_n + D*sqrt(2*h/tau_n)*randn(N))

        x = x + ((1 - x)/tau_rec - U*x*m)*h
        m = m + (- m + log(1 + exp((dot(W, U*x*m) + I0 +
                Iex[t]*cos(th - ThetaEs[t]) + Inoise))))*(h/tau)

        if t % int(pltSampl/h) == 0:
            ActM[int(t/(pltSampl/h))] = m
            ActX[int(t/(pltSampl/h))] = x
        exactR[t+1] = sum(exp(1j*th)*m)/N
        if CalcMode != 3:
            for idx, Nread in enumerate(Nreads):
                R[t+1, idx] = R[t, idx] + (-R[t, idx] + sum(exp(1j*th[choises[idx]])*poisson(m[choises[idx]]*h, Nread)))*h/tau_r

# error estimation function
def estimErrDiffNread(lag=0):
    errEsR = zeros((len(Nreads), Nev))

    for i, (t, theta) in enumerate(zip(inpTimes, inpTheta)):
        for idx, Nread in enumerate(Nreads):
            nev = abs(angle(R[t + lag:int(t + lag + T/h), idx]) - ThetaEs[t:int(t + T/h)])
            anglDiff = amin([nev, 2*pi - nev], axis=0)
            errEsR[idx, i] = mean(anglDiff)

    return errEsR

def estimErrExactR(lag=0):
    errEsR = zeros(Nev)
    for i, (t, theta) in enumerate(zip(inpTimes, inpTheta)):
        nev = abs(angle(exactR[t + lag:int(t + lag + T/h)]) - ThetaEs[t:int(t + T/h)])
        anglDiff = amin([nev, 2*pi - nev], axis=0)
        errEsR[i] = mean(anglDiff)
    return errEsR

def calcErrDiffNread():
    lags = arange(0, int(T/h), 1, dtype='int')

    errR = zeros((len(lags), len(Nreads), Nev))
    errExactR = zeros((len(lags), Nev))

    for j, lag in enumerate(lags):
        errR[j] = estimErrDiffNread(lag)
        errExactR[j] = estimErrExactR(lag)
    return errR, errExactR

if CalcMode == 0:
    integrate()
#%%
    # plotting results
    figure(figsize=(4*4, 3*4))
    gs = GridSpec(4, 2, height_ratios=[1, 1, 1, 1./3], width_ratios=[30, 1])
    gs.update(wspace=0.1, hspace=0.1)
    axM = subplot(gs[0, 0])
    axSpec = subplot(gs[1, 0], sharex=axM)
    axAngle = subplot(gs[2, 0], sharex=axM)
    axEx = subplot(gs[3, 0], sharex=axM)
    axCbar = subplot(gs[:, 1])

    egg = axM.pcolormesh(arange(0, SimTime, pltSampl), th*360/(4*pi), ActM.T)
    spam =colorbar(egg, cax=axCbar)
    axM.set_ylim([-90, 90])
    axM.set_ylabel(r"$\theta$")
    setp(axM.get_xticklabels(), visible=False)
    axM.set_title('$U={}\quad I_0={:.2f}$'.format(U, I0))

    axSpec.plot(stime, abs(exactR))
#    axSpec.plot(stime, abs(R[:, -1]))
#    axSpec.plot(arange(0, SimTime, pltSampl), mean(ActM, axis=1))
    setp(axSpec.get_xticklabels(), visible=False)
    axSpec.set_ylabel(r"$|R|$")

    axAngle.plot(stime, angle(exactR)*360/(4*pi))
#    axAngle.plot(stime, angle(R[:, -1])*360/(4*pi))
    axAngle.hlines(inpTheta*360/(4*pi), inpTimes*h, inpTimes*h + T, 'r', linewidth=3)
    setp(axAngle.get_xticklabels(), visible=False)
    axAngle.set_ylabel(r"$angle(R)$")

    axEx.plot(stime, Iex, 'r')
    axEx.set_ylim([0, C + 1])
    axEx.set_xlim([0, 50])
    axEx.locator_params(axis = 'y', nbins=3)
    axEx.set_xlabel(r"$Time [s]$")
#%%
    errR, errExactR = calcErrDiffNread()

    # errR shape lags, Nreads, Number of stimulus
    minLags = argmin(mean(errR, axis=2), axis=0)*h*1000
    minLagErr = amin(mean(errR, axis=2), axis=0)*360/(4*np.pi)

    print("Lags for minimal error {} ms".format(minLags))
    print("Minimal errors {}".format(minLagErr))

    minLag = argmin(mean(errExactR, axis=-1), axis=0)*h*1000
    minErr = amin(mean(errExactR, axis=-1), axis=0)*360/(4*np.pi)

    print("Lag for minimal error {} ms (exact PV)".format(minLag))
    print("Minimal error {} (exact PV)".format(minErr))

    df = mean(ma.array(abs(R[:, -1]), mask=Iex))
#    print("Mean value between intervals {}".format(mean(df)))

    df = mean(ma.array(abs(R[:, -1]), mask=~array(Iex, dtype='bool')))
#    print("Mean value when stimulus apply {}".format(mean(df)))
#%%
elif CalcMode == 1:
    for U, I0 in zip(Urange, Ierange):
        integrate()
        print("Calculating for U: {}".format(U))

        errR, errExactR = calcErrDiffNread()
        np.save(folderName + 'U_{:.2f}_C_{:.1f}_N_{:n}_SimTime_{:n}.npy'.format(U, C, N, SimTime), errR)
elif CalcMode == 2:
    import sys
    U = float(sys.argv[1])
    I0 = Ierange[int(U/0.05) - 1]
#    U = 0.0
#    I0 = Ierange[0]
    integrate()

    errR, errExactR = calcErrDiffNread()
    np.save(folderName + 'U_{:.2f}_C_{:.1f}_N_{:n}_SimTime_{:n}.npy'.format(U, C, N, SimTime), errR)
    np.save(folderName + 'U_{:.2f}_C_{:.1f}_N_{:n}_SimTime_{:n}_exactPV.npy'.format(U, C, N, SimTime), errExactR)
elif CalcMode == 3:
#    Iex[:] = 0
    fname = 'U_Iex_SimTime_{:.1f}_h_{:.4f}_D_{:.1f}_N_{:n}_eps_{:.3f}_m_{:.1f}.npy'

    Eps = .01
    meanAct = 0.5

    Imax = 15.0
    Imin = -15.0
    Ilow, Ihigh = Imin, Imax

    Ierange = zeros_like(Urange) + np.nan

    for idx, U in enumerate(Urange):
        for j in xrange(100):
            m = zeros(N)
            x = ones(N)
            Inoise = zeros(N)

            seed(1)

            integrate()
            cmeanAct = mean(ActM)

            if abs(meanAct - cmeanAct) < Eps:
                Ilow, Ihigh = Imin, Imax
                Ierange[idx] = I0
                break
            if (cmeanAct - meanAct) < 0:
                Ilow = I0
            else:
                Ihigh = I0
            I0 = (Ihigh + Ilow)*0.5
        print U, I0, mean(ActM)
    save(fname.format(SimTime, h, D, N, Eps, meanAct), Ierange)
