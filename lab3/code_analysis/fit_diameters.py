"""
Script to fit things
Astro 121, Spring 2014
Aaron Tran

Compute fringe modulating function
Brute-force least squares fit to obtain
angular diameters of sun and moon.
"""


import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.integrate

def main():
    ffreqs = np.linspace(0,10,100)
    mfs = np.array([modulating_funct(freq) for freq in ffreqs])

    plt.plot(ffreqs, mfs)
    plt.show()


def modulating_funct(ffreq):
    """Fringe modulating function"""
    R = 1
    x = np.linspace(-R,R,1000)  # Ranges over the width of object, finely
    fx = 1./R * np.sqrt(R**2 - x**2) * np.cos(ffreq*x)
    return sp.integrate.simps(fx, x)


if __name__ == '__main__':
    main()