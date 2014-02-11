import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def main():
    rlc_filter_gain_plot()
    diode_rc_lowpass_gain()

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
    apply_gain_plot_settings()
    plt.savefig('rlc_filter_gain.pdf')
    plt.show()

def diode_rc_lowpass_gain():
    R, C = 15, 100e-9
    
    x = np.logspace(0,7,200) # Cut-off is ~106 kHz
    w = 2*np.pi*x
    y = 1./(1 + 1j*w*R*C)
    
    plt.axvspan(20, 20000, hold=True, alpha=0.3, facecolor='g') # Audible range
    plt.plot(x, abs(y), '-k')
    apply_gain_plot_settings()
    plt.savefig('diode_rc_lowpass_gain.pdf')
    plt.show()

def apply_gain_plot_settings():
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Frequency (Hz)', fontsize=16)
    plt.ylabel(r'Gain $\left|V_{out}/V_{in}\right|$ (-)', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout() # Prevent clipping labels


if __name__ == '__main__':
    main()

