import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def main():
    rlc_filter_gain_plot()

def rlc_filter_gain_plot():
    R, L, C = 33, 1e-6, 22e-9

    x = np.logspace(5,7,200) # centered on f ~ 10^6 Hz
    w = 2*np.pi*x
    zlc = (1j*w*L) / (1 - w**2 * L * C)
    y = zlc / (R + zlc)
    # y = 1/(1+ R*(L-(2*np.pi*x)**2 * L * C)/((2*np.pi*x)*L*1j))
    
    plt.axvspan(1.045e6-0.1e6, 1.045e6+0.1e6, hold=True, alpha=0.3,
            facecolor='g')
    plt.plot(x, abs(y), '-k')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel(r'Gain $V_{out}/V_{in}$ (-)')
    plt.show()

if __name__ == '__main__':
    main()

