def precess(ra, dec, equinox1, equinox2, FK4 = None,radian=False):
    """
    NAME:
       PRECESS
    PURPOSE:
       Precess coordinates from EQUINOX1 to EQUINOX2.  
    EXPLANATION:
       For interactive display, one can use the procedure ASTRO which calls 
       PRECESS or use the /PRINT keyword.   The default (RA,DEC) system is 
       FK5 based on epoch J2000.0 but FK4 based on B1950.0 is available via 
       the /FK4 keyword.

    CALLING SEQUENCE:
       PRECESS, ra, dec, [ equinox1, equinox2, /PRINT, /FK4, /RADIAN ]

    INPUT - OUTPUT:
       RA - Input right ascension (scalar or vector) in DEGREES, unless the 
               /RADIAN keyword is set
       DEC - Input declination in DEGREES (scalar or vector), unless the 
               /RADIAN keyword is set
               
       The input RA and DEC are modified by PRECESS to give the 
       values after precession.

    OPTIONAL INPUTS:
       EQUINOX1 - Original equinox of coordinates, numeric scalar.  If 
               omitted, then PRECESS will query for EQUINOX1 and EQUINOX2.
       EQUINOX2 - Equinox of precessed coordinates.

    OPTIONAL INPUT KEYWORDS:
       FK4   - If this keyword is set, the FK4 (B1950.0) system
               will be used otherwise FK5 (J2000.0) will be used instead.
       RADIAN - If this keyword is set and non-zero, then the input and 
               output RA and DEC vectors are in radians rather than degrees

    RESTRICTIONS:
       Accuracy of precession decreases for declination values near 90 
       degrees.  PRECESS should not be used more than 2.5 centuries from
       2000 on the FK5 system (1950.0 on the FK4 system).

    EXAMPLES:
       (1) The Pole Star has J2000.0 coordinates (2h, 31m, 46.3s, 
               89d 15' 50.6") compute its coordinates at J1985.0
       In [1]: ra, dec = precess(ten(2,31,46.3)*15, ten(89,15,50.6), 2000, 1985)

               ====> 2h 16m 22.73s, 89d 11' 47.3"

       (2) Precess the B1950 coordinates of Eps Ind (RA = 21h 59m,33.053s,
       DEC = (-56d, 59', 33.053") to equinox B1975.

       In [2]: ra = ten(21, 59, 33.053)*15
       In [3]: dec = ten(-56, 59, 33.053)
       In [4]: ra,dec = precess(ra, dec ,1950, 1975, /fk4)

    PROCEDURE:
       Algorithm from Computational Spherical Astronomy by Taff (1983), 
       p. 24. (FK4). FK5 constants from "Astronomical Almanac Explanatory
       Supplement 1992, page 104 Table 3.211.1.

    PROCEDURE CALLED:
       Function PREMAT - computes precession matrix 

    REVISION HISTORY
       Written, Wayne Landsman, STI Corporation  August 1986
       Correct negative output RA values   February 1989
       Added /PRINT keyword      W. Landsman   November, 1991
       Provided FK5 (J2000.0)  I. Freedman   January 1994
       Precession Matrix computation now in PREMAT   W. Landsman June 1994
       Added /RADIAN keyword                         W. Landsman June 1997
       Converted to IDL V5.0   W. Landsman   September 1997
       Converted to Python                   April 2014
    """
    import numpy as np
    from premat import premat 
    deg_to_rad = np.pi/180.

    #Is RA a vector or scalar?
    try:
        npts = np.min( [len(ra), len(dec)] )
    except:
        npts = 1


    if not radian:
          ra_rad = ra*deg_to_rad    
          dec_rad = dec*deg_to_rad 
    else:
        ra_rad= float(ra) ; dec_rad = float(dec)
    
    sec_to_rad = deg_to_rad/3600

    a = np.cos(dec_rad)  
    x = np.zeros((3,npts))
    x[0,:] = a*np.cos(ra_rad)
    x[1,:] = a*np.sin(ra_rad)
    x[2,:] = np.sin(dec_rad)

    # Use PREMAT function to get precession matrix from Equinox1 to Equinox2
    r = premat(equinox1, equinox2, fk4 = FK4)
    x2 = np.dot(r,x)      #rotate to get output direction cosines

    ra_rad = np.arctan2(x2[1,:],x2[0,:])
    dec_rad = np.arcsin(x2[2,:])

    if not radian:
        ra = ra_rad/deg_to_rad
        ra = ra + (ra < 0.)*360.            #RA between 0 and 360 degrees
        dec = dec_rad/deg_to_rad
    else:
        ra = ra_rad ; dec = dec_rad

    
    return (ra, dec)

