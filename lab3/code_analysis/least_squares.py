"""
Homebrewed least squares fitting for Astro 121
Astro 121, Spring 2014
Aaron Tran
"""

import numpy as np


def main():
    pass


def sin_cos_lsq(x, y, k):
    """
    Returns [a, b] s.t. a*cos(k*x) + b*sin(k*x) is best fit
    
    Wrapper method for lsq_solve
    
    Input:
        x (np.array): abscissa
        y (np.array): ordinate
        (yes I know this is maximally unhelpful)
        k (float): angular frequency of fitted sinusoid
    Output:
        (see lsq_solve output)
        a: coefficients, 1-D np.array of length x_mat.shape[1] (num. columns)
        s2: sum of squared residuals
        a_vars: variances in coefficients
    """
    cos_x = np.cos(k*x)
    sin_x = np.sin(k*x)
    x_mat = np.vstack((cos_x,sin_x)).T
    
    return lsq_solve(x_mat, y)


def poly(cffs):
    """Cute one-liner to get polynomial function
    Coeffs are ordered in increasing degree (c_0, c_1, c_2, ..., c_max)
    """
    return np.vectorize(lambda x: sum([c*(x**n) for c,n in zip(cffs,xrange(len(cffs)))]))


def poly_lsq(x, y, n=2):
    """n-degree polynomial fit to data
    
    Wrapper method for lsq_solve
    
    Input:
        x (np.array): abscissa
        y (np.array): ordinate
        (yes I know this is maximally unhelpful)
        n (int): degree of polynomial fit, must be at least 0
    Output:
        (see lsq_solve output)
        a: coefficients, 1-D np.array of length x_mat.shape[1] (num. columns)
        s2: sum of squared residuals
        a_vars: variances in coefficients
    """
    x_mat = np.ones((x.size, n+1), dtype='float')
    for i in xrange(n+1):
        x_mat[:,i] = x**i
    return lsq_solve(x_mat, y)


def lin_lsq(x, y):
    """Returns [a,b] s.t. y = a*x + b is a least squares best fit
    
    Wrapper method for lsq_solve
    
    Input:
        x (np.array): abscissa
        y (np.array): ordinate
        (yes I know this is maximally unhelpful)
        n (int): degree of polynomial fit, must be at least 0
    Output:
        (see lsq_solve output)
        a: coefficients, 1-D np.array of length x_mat.shape[1] (num. columns)
        s2: sum of squared residuals
        a_vars: variances in coefficients
    """
    #x_mat = np.ones((x.size, 2))
    #x_mat[:,0] = x  # Assign first column to measured x values
    #return lsq_solve(x_mat, y)
    return poly_lsq(x, y, n=1)


def lsq_solve(x_mat, y):
    """Least squares solution to x_mat * a = y, a = coefficients
    
    Define alpha = X.T * X, beta = X.T * y (X = X_mat)
    Solve alpha * a = beta for a
    
    Input:
        x_mat (2-D np.array): design matrix
        y (1-D np.array): y values to be fitted
    Output:
        a: coefficients, 1-D np.array of length x_mat.shape[1] (num. columns)
        s2: sum of squared residuals
        a_vars: variances in coefficients
    """
    alpha = np.dot(x_mat.T, x_mat)
    beta = np.dot(x_mat.T, y)
    
    a = np.linalg.solve(alpha, beta)
    cov = np.linalg.inv(alpha)  # covariance matrix = inv(X.T * X), not normalized
    
    y_fit = np.dot(x_mat, a)
    s2 = sum_sq_residuals(y, y_fit, y.size - 2)
    a_vars = s2 * np.diag(cov)
    
    return a, s2, a_vars


def sum_sq_residuals(y, y_model, dof):
    """Sum of squared residuals, divided by # degrees of freedom"""
    y_res = y - y_model
    return 1./(dof) * np.sum(np.power(y_res,2))


if __name__ == '__main__':
    main()