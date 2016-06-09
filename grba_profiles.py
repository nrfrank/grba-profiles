import sys
import time
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import nquad
from scipy.special import lambertw

def intG(y, chi, k = 0.0, p = 2.2):
    bG = (1.0 - p)/2.0
    ys = np.power(y, 0.5*(bG*(4.0 - k) + 4.0 - 3.0*k))
    chis = np.power(chi, np.divide(7.0*k - 23.0 + bG*(13.0 + k), 6.0*(4.0 - k)))
    factor = np.power((7.0 - 2.0*k)*chi*np.power(y, 4.0 - k) + 1.0, bG - 2.0)
    return ys*chis*factor

def fluxG(y, chi, k = 0.0, p = 2.2):
    Ck = (4.0 - k)*np.power(5.0 - k, np.divide(k - 5.0, 4.0 - k))
    cov = np.divide(np.power(y, 5.0 - k), 2.0*Ck)
    return 2.0*np.pi*cov*intG(y, chi, k, p)

def bounds_Y(chi, k = 0.0):
    return [0.0, np.power(chi, np.divide(1.0, k - 4.0))]

def bounds_Chi():
    return [1.0, np.inf]

def fluxG_oaStr(y, chi, kap, sig, gA = 1.0, k = 0.0, p = 2.2):
    Gk = (4.0 - k)*gA**2.0
    C_const = np.divide(np.power(np.log(2.0), 1.0/kap), (4.0 - k)*gA**2.0*sig**2.0)
    C_yChi = np.divide(y - chi*np.power(y, 5.0 - k), y**2.0)
    C = C_const*C_yChi
    #C = (np.power(np.log(2.0), 1.0/kap)/((4.0 - k)*np.power(gA*sig, 2.0)))*((y - chi*np.power(y, 5.0 - k))*np.power(y, -2.0))
    lamb = np.real(lambertw(-np.power(C, kap)*kap))
    denom = (1.0 + lamb)*np.exp(lamb/kap)
    cov = np.divide(np.power(y, 5.0 - k), 2.0*Gk)
    return 2.0*np.pi*cov*intG(y, chi)/denom

def bounds_y(chi, t0, t1, k = 0.0):
    # return [0.0, 1.0]
    return [0.0, np.power(chi, np.divide(1.0, k - 4.0))]

def bounds_chi(t0, t1):
    #return [1.0, np.power(y, k - 4.0)]
    return [1.0, np.inf]

def thetaPrime(r, thv, phi):
    top = r*(np.cos(thv)**2.0 - 0.25*np.sin(2.0*thv)**2.0*np.cos(phi)**2.0)**2.0
    bot = 1.0 + 0.5*r*np.sin(2.0*thv)*np.cos(phi)
    return np.divide(top, bot)

def r0_max(y, kap, sig, thv, gA = 1.0, k = 0.0, p = 2.2):
    Gk = (4.0 - k)*gA**2.0
    def root(rm):
        thP0 = thetaPrime(rm, thv, 0.0)
        rExp = -np.power(np.divide(thP0, sig), 2.0*kap)
        lhs = np.divide(y - np.power(y, 5.0 - k), Gk)
        rhs = (np.tan(thv) + rm)**2.0*np.exp2(rExp)
        return rhs - lhs

    rootVal = fsolve(root, 0.4)[0]
    return rootVal

def fluxG_fullStr(r, y, kap, sig, thv, gA = 1.0, k = 0.0, p = 2.2):
    Gk = (4.0 - k)*gA**2.0
    thP0 = thetaPrime(r, thv, 0.0)
    exp0 = np.power(np.divide(thP0, sig), 2.0*kap)
    chiVal = np.divide(y - Gk*np.exp2(-exp0)*(np.tan(thv) + r)**2.0, np.power(y, 5.0 - k))
    return 2.0*np.pi*r*intG(y, chiVal)

def bounds_yr(kap, sig, thv):
    return [0.0, 1.0]

def bounds_ry(y, kap, sig, thv):
    return [0.0, r0_max(y, kap, sig, thv)]

def main():
    # base_int = nquad(fluxG, [bounds_Y, bounds_Chi])
    tiny = np.power(10.0, -3.0)
    SIGMA = 2.0
    #KAPPA = tiny
    for kap in range(10):
        KAPPA = np.power(10.0, -float(kap + 1))
        #KAPPA = float(kap) + 0.01
        #oa_str_int = nquad(fluxG_oaStr, [bounds_chi, bounds_y], args = (KAPPA, SIGMA))[0]
        oa_str_int = nquad(fluxG_oaStr, [bounds_y, bounds_chi], args = (KAPPA, SIGMA))[0]
        str_int = nquad(fluxG_fullStr, [bounds_ry, bounds_yr], args = (KAPPA, SIGMA, 0.0))[0]
        print oa_str_int, str_int, np.abs(oa_str_int - str_int)/str_int*100.0


if __name__ == "__main__":
    sys.exit(int(main() or 0))