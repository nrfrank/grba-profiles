import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import seaborn as sns
from grba_profiles import intG, thetaPrime, fluxG_fullStr, r0_max

def rootFuncR(R, R0, kap, sig, thv, phi):
    thP = thetaPrime(R, thv, phi)
    lExp = -np.power(np.divide(thP, sig), 2.0*kap)
    left = (np.tan(thv)**2.0 + 2.0*R*np.tan(thv)*np.abs(np.cos(phi)) + R**2.0)*np.exp2(lExp)
    thP0 = thetaPrime(R0, thv, 0.0)
    rExp = -np.power(np.divide(thP0, sig), 2.0*kap)
    right = (np.tan(thv) + R0)**2.0*np.exp2(rExp)
    return(left - right)

N = 100
KAPPA = 1.0
SIGMA = 2.0
THETA_V = np.radians(0.0)
Y_VAL = 0.7
R_MAX = r0_max(Y_VAL, KAPPA, SIGMA, THETA_V, gA = 1.0, k = 0.0, p = 2.2)

print R_MAX

rVals = np.linspace(0.0, R_MAX, N)
pVals = np.linspace(0.0, 2.0*np.pi, N)
P, R = np.meshgrid(pVals, rVals)
RP = np.zeros((N, N))

start = time.time()
for i in range(N):
    guess = R[i][0]
    for j in range(N):
        rP = fsolve(rootFuncR, guess, args = (R[i][j], KAPPA, SIGMA, THETA_V, P[i][j]))
        # print R[i][j], P[i][j], rP
        if rP == 0.0 and R[i][0] == 0.0:
            RP[i][j] = 1.0
        else:
            RP[i][j] = (rP/R[i][j])**2.0
        guess = rP

stop = time.time()
runtime = stop - start
print "runtime = ", runtime, ", final R' value = ", rP

plt.figure()
# RP = ma.masked_invalid(RP)
plt.pcolormesh(P, R, RP, cmap = plt.cm.viridis)
plt.colorbar()
# plt.axhline(y = np.pi/2.0, color = 'r')
# plt.axhline(y = np.pi*3.0/2.0, color = 'r')
plt.axis([min(pVals), max(pVals), min(rVals), max(rVals)])
plt.show()