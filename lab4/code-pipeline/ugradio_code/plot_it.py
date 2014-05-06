#! /usr/bin/env python

import readspec_mod
import pylab as p, sys, numpy as n

for filename in sys.argv[1:]:
    print 'REading', filename
    d = readspec_mod.readSpec(filename)
    print d.shape
    d = n.average(d, axis=-1)
    p.plot(d)
p.show()
