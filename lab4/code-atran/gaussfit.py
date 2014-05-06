"""
Module for performing SINGLE Gaussian fit
I wanted to throw something for multiple Gaussians
but ain't got time and not useful
Utility module for Astro 121, Spring 2014
Aaron Tran
"""

import numpy as np

import least_squares as lsq

def main():
    """Test of stuff"""
    # gaussian(x, x0, a std)
    x = np.linspace(-2, 2, 100)
    y = gaussian(x, x0=3, a=1, sigma=3) + 0.5 * (np.random.rand(x.size) - 0.5)
    gaussfit(x, y, verbose=True)


def gaussfit(x, y, rconvf=1e-7, verbose=False):
    """Non-linear single Gaussian fit, for one Gaussian (with offset = 0)

    y(x) = A exp[ -(x-x0)^2 / (2 sigma^2) ]

    Input:
        x (np.array): x values to be fitted
        y (np.array): y values to be fitted
        rconvf (float): convergence factor for sum sq. residuals
        verbose (bool): babble about the fitting procedure
    Output:
        x0, a, std best fit values
    """
    # Initial guesses from linear lsq fits in log space
    # ln(y) = -(1/(2 sigma^2)) x^2 + (x0 / sigma^2) x + ln(A) - (x_0)^2
    # a_2 = -1/(2 sigma^2), a_1 = x0 / sigma^2
    a0, a1, a2 = lsq.poly_lsq(x, np.log(y), n=2)[0]
    std_init = np.sqrt(-1./(2*a2))
    x0_init = a1 * std_init**2
    a_init = np.exp(a0 + x0_init**2)

    # Set up interation variables
    x0_fit, x0_fit_p = x0_init, 0
    a_fit, a_fit_p = a_init, 0
    std_fit, std_fit_p = std_init, 0
    s2, s2_p = 1, 0

    # Iterate until stuff converges
    while abs(s2 - s2_p) > s2 * rconvf:
        # Linearization, Gauss-Newton essentially (see lsq_lite)
        y_fit = gaussian(x, x0_fit, a_fit, std_fit)  # previous fit
        dy = y - y_fit
        da_coeff = y_fit / a_fit
        dx0_coeff = y_fit * (x-x0_fit) / (std_fit**2)
        dstd_coeff = y_fit * np.power(x-x0_fit,2) / (std_fit**3)
        # Construct design matrix for linear least squares
        xmat = np.vstack((da_coeff, dx0_coeff, dstd_coeff)).T
        # Solve for changes in fit parameters
        da, dx0, dstd = lsq.lsq_solve(xmat, dy)[0]

        # Update values of std, x0, a, following linearization
        std_fit_p = std_fit
        x0_fit_p = x0_fit
        a_fit_p = a_fit
        std_fit = std_fit_p + dstd
        x0_fit = x0_fit_p + dx0
        a_fit = a_fit_p + da

        # Check fit residual
        y_fit = gaussian(x, x0_fit, a_fit, std_fit)
        s2_p = s2
        s2 = lsq.sum_sq_residuals(y, y_fit, 3)
        if verbose:
            print 'Sum sq. residuals is %.8f' % s2

    if verbose:
        print 'Best fit x0 = %g, a = %g, std = %g' % (x0_fit, a_fit, std_fit)

    return x0_fit, a_fit, std_fit


def gaussian(x, x0, a, sigma):
    """Computes Gaussian for vector x of values"""
    return a * np.exp( -np.power(x-x0,2) / (2*sigma**2))


if __name__ == '__main__':
    main()
