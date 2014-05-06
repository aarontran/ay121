
README
======

Thrown together very quickly (May 6, 2014), Aaron Tran.

Scripts may not work.  If they are missing ugradio imports, you may need to add:

    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                 'ugradio_code'))

Main scripts
------------

    spec\_veloc.py
    spec\_average.py
    single\_pt\_leuschner.py

spec\_veloc.py is meant to fit in the data pipeline after Caleb's data
processing scripts, and relies on the directory structure implemented on
heiles (the remote computer used to control Leuschner).  This has been
superseded by Isaac's much more complete get\_img\_data.py and other routines.
More information / directory structure, etc. are on Isaac's pipeline at
[bitbucket](https://bitbucket.org/domagalski/cia-leuschner-pipeline/src).

single\_pt\_leuschner.py is meant built for quickly taking data on Leuschner
at a single pointing.  No capability for continuous tracking, best for debugging
and testing of Leuschner code.  Does perform checking to see if target of interest (in l,b) is visible.

Supporting files
----------------

For astronomical coordinates, least squares fitting, pickle management

    catalog.py
    coord\_conv.py
    gaussfit.py
    least\_squares.py
    nps\_ephem.py
    pklio.py

All just utility things.