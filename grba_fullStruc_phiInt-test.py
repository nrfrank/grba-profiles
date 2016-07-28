import sys
import numpy as np
import numpy.ma as ma
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import matplotlib.colors as c
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from grba_profiles import intG, thetaPrime, r0_max
from time import time

def rootFuncR(R, R0, kap, sig, thv, phi):
    thP = thetaPrime(R, thv, phi)
    lExp = -np.power(np.divide(thP, sig), 2.0*kap)
    left = (np.tan(thv)**2.0 + 2.0*R*np.tan(thv)*np.cos(phi) + R**2.0)*np.exp2(lExp)
    thP0 = thetaPrime(R0, thv, 0.0)
    rExp = -np.power(np.divide(thP0, sig), 2.0*kap)
    right = (np.tan(thv) + R0)**2.0*np.exp2(rExp)
    return(left - right)

def fluxG_fullStr(r, y, kap, sig, thv, gA = 1.0, k = 0.0, p = 2.2):
    Gk = (4.0 - k)*gA**2.0
    thP0 = thetaPrime(r, thv, 0.0)
    exp0 = np.power(np.divide(thP0, sig), 2.0*kap)
    chiVal = np.divide(y - Gk*np.exp2(-exp0)*(np.tan(thv) + r)**2.0, np.power(y, 5.0 - k))
    return r*intG(y, chiVal)

def surfacePlot_r0Phi(num, yVal, kap, sig, thv, gA = 1.0, k = 0.0, p = 2.2):
    R0_MAX = r0_max(yVal, kap, sig, thv, gA, k, p)
    rVals = np.linspace(0.0, R0_MAX, num)
    pVals = np.linspace(0.0, 2.0*np.pi, num)
    P, R = np.meshgrid(pVals, rVals)
    F_ = np.zeros((num, num))
    for i in range(num):
        # flux = []
        guess = R[i][0]
        for j in range(num):
            rP = fsolve(rootFuncR, guess, args = (R[i][j], kap, sig, thv, P[i][j]))
            # print R[i][j], P[i][j], rP
            if rP == 0.0 and R[i][0] == 0.0:
                f = 1.0
            else:
                f = (rP/R[i][j])**2.0
            guess = rP
            # F_[i][j] = f*fluxG_fullStr(R[i][j], yVal, kap, sig, thv, gA, k, p)
            F_[i][j] = f
            # flux.append(fluxG_fullStr(R[i][j], yVal, kap, sig, thv)*f)
        # flux = [f/max(flux) for f in flux]
    
    plt.figure()
    F_ = ma.masked_invalid(F_)
    plt.pcolormesh(P, R, F_, cmap = plt.cm.viridis, norm = c.LogNorm())
    plt.colorbar()
    # plt.axhline(y = np.pi/2.0, color = 'r')
    # plt.axhline(y = np.pi*3.0/2.0, color = 'r')
    plt.plot([np.pi, np.pi], [0.0, R0_MAX], 'r-')
    plt.axis([min(pVals), max(pVals), min(rVals), max(rVals)])
    plt.show()
    plt.clf()
    plt.close()
    
    return F_
    
def surfacePlot_r0Phi_3D(num, yVal, kap, sig, thv, gA = 1.0, k = 0.0, p = 2.2):
    R0_MAX = r0_max(yVal, kap, sig, thv, gA, k, p)
    rVals = np.linspace(0.0, R0_MAX, num)
    pVals = np.linspace(0.0, 2.0*np.pi, num)
    P, R = np.meshgrid(pVals, rVals)
    F_ = np.zeros((num, num))
    for i in range(num):
        # flux = []
        guess = R[i][0]
        for j in range(num):
            rP = fsolve(rootFuncR, guess, args = (R[i][j], kap, sig, thv, P[i][j]))
            # print R[i][j], P[i][j], rP
            if rP == 0.0 and R[i][0] == 0.0:
                f = 1.0
            else:
                f = (rP/R[i][j])**2.0
            guess = rP
            F_[i][j] = f*fluxG_fullStr(R[i][j], yVal, kap, sig, thv, gA, k, p)
            # flux.append(fluxG_fullStr(R[i][j], yVal, kap, sig, thv)*f)
        # flux = [f/max(flux) for f in flux]
    
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    F_ = ma.masked_invalid(F_)
    surf = ax.plot_surface(P, R, F_, cmap = plt.cm.viridis, norm = c.LogNorm())
    fig.colorbar(surf)
    plt.show()
    plt.clf()
    plt.close()
    
    return F_

def main():
    N = 1000
    TINY = np.power(10.0, -3.0)
    KAPPA = 1.0
    SIGMA = 2.0
    THETA_V = np.radians(6.0)
    Y_VAL = 0.1
    R_MAX = r0_max(Y_VAL, KAPPA, SIGMA, THETA_V, gA = 1.0, k = 0.0, p = 2.2)
    F_PHI = surfacePlot_r0Phi(N, Y_VAL, KAPPA, SIGMA, THETA_V)
    # F_PHI = surfacePlot_r0Phi_3D(N, Y_VAL, KAPPA, SIGMA, THETA_V)

if __name__ == "__main__":
    sys.exit(int(main() or 0))