"""
Astro 121 Radio Lab
Spring 2014, Prof. A Parsons
Aaron Tran

Obtain n sums of m random samples, from a random distribution
Show that sums are normally distributed by central limit theorem
Show that distribution of sums has variance ~ 1/m
Compare to theoretical expectation

Some magic numbers/etc floating around, but works well enough
"""

import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    REFERENCE
    m = number of random samples to sum
    n = number of sums to calculate
    bins = number of bins for histogram
    """
    
    b = 100 # Number of histogram bins
    
    # Fix m = 10.  Vary n to show that means /are/ normally distributed
    # This shows that getting more samples (increasing n),
    # better approximates the underlying distribution
    m = 10
    nvals = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1e3, 2e3, 5e3, 1e4, 1e5, 1e6]
    make_gaussian_conv_plots(m, nvals)
    
    # Fix n = 1e5.  Vary m to show that standard deviation tightens
    # as you increase m (number of samples summed/averaged)
    # Expect: (std of underlying distribution) / sqrt(m)
    n = 1e5
    mvals = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
    make_stdev_plot(mvals, n, b)
    
    # Make some representative plots!
    make_normal_plot(10, 1e1, 100)
    make_normal_plot(10, 1e3, 100)
    make_normal_plot(10, 1e5, 100)
    make_normal_plot(10, 1e5, 100)
    make_normal_plot(100, 1e5, 100)
    make_normal_plot(1000, 1e5, 100)

def make_gaussian_conv_plots(m, nvals):
    """
    Compute mean of means, std of means for many values of n
    Plot against expected mean, std of mean as function of n
    """
    mean_of_means = np.vectorize(lambda x: np.mean(mean_unif_samp(m,x)))
    std_of_means = np.vectorize(lambda x: np.std(mean_unif_samp(m,x)))
    
    mns = mean_of_means(nvals)
    stds = std_of_means(nvals)
    
    mu = 0.5
    std_expect = std_of_mean_unif(m)
    
    plt.plot(nvals,mns, 'ko')
    plt.axhline(mu)
    plt.xscale('log')
    plt.xlabel('Number of %d-sample sums' % (m))
    plt.ylabel('Mean of means')
    plt.show()
    
    plt.plot(nvals,stds, 'ko')
    plt.axhline(std_expect)
    plt.xscale('log')
    plt.xlabel('Number of %d-sample sums' % (m))
    plt.ylabel('Standard deviations of means')
    plt.show()


def make_stdev_plot(mvals, n, b):
    """
    Compute std of means for various values of m
    Plot against expected standard deviation as a function of m
    """
    std_of_means = np.vectorize(lambda x: np.std(mean_unif_samp(x, n)))
    
    stds = std_of_means(mvals)
    
    m_expect = np.logspace(0,3,200) # Expected stdev. scaling
    std_expect = std_of_mean_unif(m_expect)
    
    plt.plot(mvals, stds, 'ko')
    plt.plot(m_expect, std_expect, '--k')
    plt.xscale('log')
    plt.xlabel('Number of samples in computed mean')
    plt.ylabel(r'Standard deviation of mean with $n=%d$ means' % (n))
    plt.savefig('allan_variance.pdf')
    plt.show()


def make_poisson_plot(m,n,bins):
    """Samples from Poisson dist. to make hist., predicted Gaussian"""
    data_means = mean_poisson_samp(m,n)
    std = std_of_mean_poisson(m) # std, mu are /expected/ values
    mu = 0.5
    
    x = np.linspace(mu-5*std, mu+5*std, 100)
    y = n*(10*std)/bins * gaussian(x, mu, std) # Scale by n * (bin width)
    
    plt.plot(x, y, '-k', linewidth=2.0)
    plt.hist(data_means, bins, (mu-5*std,mu+5*std), color='c')
    plt.xlabel('Mean of means')
    plt.ylabel('Number of occurences')
    plt.savefig('hist_poisson_m=%d_n=%d.pdf' % (m,n))
    plt.show()


def make_normal_plot(m,n,bins):
    """Samples from unif. dist. to make hist., theoretical prediction"""
    data_means = mean_unif_samp(m,n)
    std = std_of_mean_unif(m) # std, mu are /expected/ values
    mu = 0.5
    
    x = np.linspace(mu-5*std, mu+5*std, 100)
    y = n*(10*std)/bins * gaussian(x, mu, std) # Scale by n * (bin width)
    
    plt.plot(x, y, '-k', linewidth=2.0)
    plt.hist(data_means, bins, (mu-5*std,mu+5*std), color='c')
    plt.xlabel('Mean of means')
    plt.ylabel('Number of occurences')
    plt.savefig('hist_normal_m=%d_n=%d.pdf' % (m,n))
    plt.show()


def gaussian(x, mean, std):
    # Take input vector x, output normalized gaussian(x)
    power = np.power(x-mean,2) / (2.*(std**2))
    return 1./(std*np.sqrt(2*np.pi)) * np.exp(-power)

def mean_poisson_samp(m,n):
    # Returns n sums of m samples from Poisson dist., \lambda = 0.5
    x = np.random.poisson(0.5, (m,n))
    return np.mean(x, 0)

def std_of_mean_poisson(m):
    # Expected standard deviation of m-sample mean
    # from Poisson dist., \lambda = 0.5
    poissonstd = np.sqrt(0.5)
    return poissonstd/np.sqrt(m)

def mean_unif_samp(m,n):
    # Returns n sums of m samples from uniform dist. [0,1)
    x = np.random.uniform(0, 1, (m,n))
    return np.mean(x, 0)

def std_of_mean_unif(m):
    # Expected standard deviation of m-sample mean
    # from uniform dist. on [0,1)
    unifstd = 1/np.sqrt(12)
    return unifstd/np.sqrt(m)

# Boilerplate to call main()
if __name__ == '__main__':
    main()
