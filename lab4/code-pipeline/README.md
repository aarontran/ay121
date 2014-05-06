CIA 21cm Leuschner Radio Pipeline
================================================

Authors: Caleb Levy, Isaac Domagalski, Aaron Tran (CIA team)

This repository contains the observing and reduction pipeline used to control
the 21cm dish at Leuschner Observatory in Lafayette, CA. This pipeline consists
of two main parts. The observing and tracking scripts, written by Aaron Tran and
Isaac Domagalski, can be found in the obs-ctrl directory. The data reduction
scripts, written by Caleb Levy and Isaac Domagalski, are located in this current
directory. See the README.md in the obs-crtl directory for a description of the
observing and tracking software used for controlling the dish and recording the
raw data that later gets reduced to extract the velocities, dispersions, and
column densities of hydrogen gas in the Milky Way galaxy.

Spectrum Calibration
----------------------------

Coming soon....

Velocity Reduction
-------------------------------------

Calculating the velocity and column density from the calibrated spectra is very
straightforward. First, the line profile is extracted out of the total signal.
This is done by fitting the background to a 5th order polynomial. The background
region used in the fit was the region outside of +/- 0.5 MHz from the maximum
frequency in the spectrum. Once the background was subtracted, the cut used to
extract the line was determined by finding where the line crossed zero. Once the
line region was determined, the velocity and dispersion were computed by taking
the weighted mean and variance of velocity, which was computed using Doppler
shift, with using the height of the spectrum as the weights. The sum of all of
the weights was also used to determine the column density of hydrogen.
Additionally, the number of peaks and their locations was computed. This was
done by smoothing the spectrum, then determining the relative maximum. In order
for a relative maximum to be a peak, the points surrounding it in the smoothed
spectra had to vary greater than twice the noise level of the spectrum.
