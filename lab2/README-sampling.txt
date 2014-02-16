Karto's set-up for
1. controlling SRS function generators (local oscillators)
2. obtaining samples from computer 'pulsar'

On a lab computer, in Python (computer should have numpy):
import DFEC
Do NOT run remotely (ssh), be here in the lab to run the sampler!!!!!

DFEC.set_srs
------------
Controls the SRS local oscillators in the lab!
Remember that SRS 2 is stacked on top of SRS 1.

First time around: set offset, phase, dbm = 0 everything
(in case someone else had edited the settings)

DFEC.sampler 
------------
Lets us get data from sampler in 'pulsar'
For sampling: 2 channel, use two BNCs on left.  1 channel, right BNC.

timeWarn -- it may start slowing / struggling if you collect too much
