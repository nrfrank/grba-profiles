import sys
import time
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import nquad

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

def bounds_y(chi, k = 0.0):
    return [0.0, np.power(chi, np.divide(1.0, k - 4.0))]

def bounds_chi():
    return [1.0, np.inf]

def main():
    start = time.time()
    integral = nquad(fluxG, [bounds_y, bounds_chi])
    stop = time.time()
    print integral, stop - start

if __name__ == "__main__":
    sys.exit(int(main() or 0))