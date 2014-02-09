"""
Astro 121 Radio Lab
Spring 2014
Aaron Tran

Obtain many sums of n random samples, from some distribution
Show that sums are normally distributed by central limit theorem
Show that distribution of sums has variance ~ 1/n
Compare to a theoretical expectation
"""

import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    Unif. dist.: mean of means = 0.5, stdev. of means = 1/sqrt(12) / sqrt(n)
    Poisson dist.: mean of means = 0.5, stdev of means = 1/sqrt(0.5) / sqrt(n)
    If plotting sums, mean of sums is m*0.5.  Stdev of sums = m/sqrt(?) /sqrt(n)
    """
    m = 10 # number of random samples to sum
    n = 1000000 # number of sums to calculate
    bins = 100 # number of bins for histogram
    
    data_means = mean_unif_samp(m,n)
    std = std_of_mean_unif(m)
    mu = 0.5
    
    x = np.linspace(mu-5*std, mu+5*std, 100)
    y = n*(10*std)/bins * gaussian(x, mu, std) # Scale by n * (bin width)
    
    plt.plot(x, n*0.1*std*y, '-k', linewidth=2.0)
    plt.hold(True)
    plt.hist(data_means, bins, (mu-5*std,mu+5*std), color='c')
    plt.savefig('test.png')

def makeNormalPlot(m,n,bins):
    print('hi')

def gaussian(x, mean, std):
    # Take input vector x, output normalized gaussian(x)
    power = np.power(x-mean,2) / (2.*(std**2))
    return 1./(std*np.sqrt(2*np.pi)) * np.exp(-power)


def mean_poisson_samp(m,n):
    # Returns n sums of m samples from Poisson dist., \lambda = 0.5
    x = np.random.poisson(0.5, (m,n))
    return np.mean(x, 0)

def std_of_mean_poisson(m):
    # Standard deviation of m-sample mean from Poisson dist., \lambda = 0.5
    poissonstd = 1/sqrt(0.5)
    return poissonstd/np.sqrt(m)

def mean_unif_samp(m,n):
    # Returns n sums of m samples from uniform dist. [0,1)
    x = np.random.uniform(0, 1, (m,n))
    return np.mean(x, 0)

def std_of_mean_unif(m):
    # Standard deviation of m-sample mean from uniform dist. on [0,1)
    unifstd = 1/np.sqrt(12)
    return unifstd/np.sqrt(m)

# Boilerplate to call main()
if __name__ == '__main__':
    main()
