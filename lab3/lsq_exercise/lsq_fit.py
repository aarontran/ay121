"""
Least squares fitting exercise
Astro 121
March 18, 2014
"""

import numpy as np
import matplotlib.pyplot as plt

def main():
    """Import and attempt to model Aaron Parson's data
    First run linear fit
    Then fit residuals of linfit to sinusoid
    Then, run a simpler combined fit

    Spits out plots and sample variances (scaled by DoFs)
    Cannot give chi^2 or reduced chi^2 without uncertainties
    """
    # 1000 pts from 0 to 10, evenly spaced
    # Need 1001 pts to have spacing 0.01 exactly
    data = np.load('ay121_lsq_data.npz')
    x, y = (data['x'], data['y'])
    n = np.size(x)
    dx = float(x[-1]) / n  # Sample spacing
    
    # Linear fit to data
    a1, b1 = lin_lsq(x, y)  # Fit coefficients
    y_fit = a1 * x + b1
    y_res = y - y_fit  # Fit residuals
    s2 = 1./(n-2) * np.sum(np.power(y_res,2))  # 2 param. fit s^2
    print 'Linfit sample variance: %g' % s2

    # Plot linear fit to data
    plt.plot(x, y, '.-k')
    plt.plot(x, y_fit, '-r', linewidth=2)
    plt.show()

    # Find signal frequency of residuals
    res_freqs = np.fft.fftfreq(n, dx)
    res_pow_spec = np.abs( np.fft.fft(y_res) )**2
    f = res_freqs[ np.argmax(res_pow_spec[:n/2]) ]
    print 'Main signal frequency at: %g' % f
    print 'Frequency resolution: %g' % res_freqs[1]
    k = 2*np.pi*f  # Convert to a wavenumber

    # Fit residuals to linear combination of sine, cosine
    a2, b2 = sin_cos_lsq(x, y_res, k)
    print a2, b2
    y_res_fit = a2 * np.cos(k*x) + b2 * np.sin(k*x)
    s2_res_fit = 1./(n-2) * np.sum(np.power(y_res - y_res_fit, 2))
    print 'Residual fit sample variance: %g' % s2_res_fit

    # Plot sin/cos fit to residuals
    plt.plot(x, y_res, '.-k')
    plt.plot(x, y_res_fit, '-r', linewidth=2)
    plt.show()

    # Simultaneously fit everything to a single model
    a3, b3, m, c = combined_lsq(x,y,k)
    y_all_fit = a3*np.cos(k*x) + b3*np.sin(k*x) + m*x + c
    plt.plot(x, y, '.k')
    plt.plot(x, y_all_fit, '-r', linewidth=2)
    plt.show()
    s2_all = 1./(n-4) * np.sum(np.power(y - y_all_fit,2))
    print 'Combined fit sample variance: %g' % s2_all


def combined_lsq(x, y, k):
    """Returns [a, b, m, c] fitting a*cos(k*x) + b*sin(k*x) + m*x + c"""
    cos_x = np.cos(k*x)
    sin_x = np.sin(k*x)
    x_mat = np.vstack((cos_x, sin_x, x, np.ones(x.size))).T
    return lsq_solve(x_mat, y)


def sin_cos_lsq(x, y, k):
    """Returns [a, b] s.t. a*cos(k*x) + b*sin(k*x) is best fit"""
    cos_x = np.cos(k*x)
    sin_x = np.sin(k*x)
    x_mat = np.vstack((cos_x,sin_x)).T
    return lsq_solve(x_mat, y)


def lin_lsq(x, y):
    """Returns [a,b] s.t. y = a*x + b is a least squares best fit"""
    x_mat = np.ones((x.size, 2))
    x_mat[:,0] = x  # Assign first column to measured x values
    return lsq_solve(x_mat, y)
        

def lsq_solve(x_mat, y):
    """Least squares solution to x_mat * a = y, a = coefficients
    x_mat = design matrix
    Solve X.T * X * a = X.T * y
    Returns a
    """
    A_mat = np.dot(x_mat.T, x_mat)
    b_vec = np.dot(x_mat.T, y)
    return np.linalg.solve(A_mat, b_vec)


if __name__ == '__main__':
    main()
