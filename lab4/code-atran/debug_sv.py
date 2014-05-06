"""
Short code (from iPython sessions) used to plot output data, after reduction
and conversion to velocities...
sanity checking my computations of mean velocities, etc..

Astro 121, Spring 2014
Aaron Tran
May 1, 2014 (i.e., while I was working on spec_veloc.py)
"""

import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt

with open('data_l_359.50_b_028.00_grid_014_066_reduce.pkl', 'rb') as f:
    a = pickle.load(f)

#with open('data_l_210.00_b_004.00_grid_002_000_reduce.pkl', 'rb') as f:
#    a = pickle.load(f)

"""
plt.plot(a['Frequency'], a['Spectrum'])
plt.axis('tight')
plt.xlabel('Frequency (Hz))')
plt.ylabel('Temperature (K)')
plt.ticklabel_format(useOffset=False, axis='x')
plt.title(a['fname'])
plt.show()

plt.plot(a['Velocity'], a['Spectrum'])
plt.axis('tight')
plt.xlabel('Velocity (km/s)')
plt.ylabel('Temperature (K)')
plt.title(a['fname'])
plt.show()
"""

a['Science']

# Idea -- show how integral varies as we move along velocity axis...
# Makes it very clear that baseline is killing us
intspec = [np.trapz(a['Spectrum'][::-1][:i] * a['Velocity'][::-1][:i], a['Velocity'][::-1][:i]) / np.trapz(a['Spectrum'][::-1], a['Velocity'][::-1]) for i in xrange(len(a['Spectrum']))]
