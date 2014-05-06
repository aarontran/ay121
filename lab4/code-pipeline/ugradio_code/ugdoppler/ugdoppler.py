import numpy as np
import radiolab as rlab
from baryvel import baryvel
from premat import premat
from precess import precess

def ugdoppler(ra, dec, julday, nlat=None, wlong=None, light=False, obspos_deg=None,lst_mean=None):
    """

    NAME: ugdoppler

    PURPOSE:
       computes the projected velocity of the telescope wrt
       four coordinate systems: geo, helio, bary, lsr.
       negative velocities mean approach

       the standard LSR is defined as follows: the sun moves at 20.0 km/s
       toward ra=18.0h, dec=30.0 deg in 1900 epoch coords

    CALLING SEQUENCE:
       vel = ugdoppler( ra, dec, julday, $
       nlat=nlat, wlong=wlong, $
       path=path, light=light, $
       obspos_deg=obspos_deg, lst_mean=lst_mean)

    INPUTS: fully vectorized...ALL THREE INPUTS MUST HAVE SAME DIMENSIONS!!
       ra[n] - the source ra in DECIMAL HOURS, equinox 2000
       dec[n] - the source dec in decimal degrees, equinox 2000
       julday[n] - the full (unmodified) julian day JD. MJD = JD - 2400000.5

    KEYWORD PARAMETERS
       nlat, wlong - specify nlat and wlong of obs in
       degrees.  if you set one, you must set the other
       also. For Leuschner, nlat=37.8732, wlong=+122.2573

       light - returns the velocity as a fraction of c

    OUTPUTS:
       program returns the velocity in km/s, or as a faction of c if
       the keyword /light is specified. the result is a 4-element
       vector whose elements are [geo, helio, bary, lsr]. quick
       comparison with phil's C doppler routines gives agreement to
       better than 100 m/s one arbitrary case.

    OPTIONAL OUTPUTS:
       obspos_deg: observatory [lat, wlong] in degrees that was used
       in the calculation. This is set by either (nlat and wlong) or
       path default is Arecibo.

       lst_mean: the lst at the observatory for the specified JD

    REVISION HISTORY: carlh 29oct04.
       from idoppler_ch changed calculation epoch to 2000
       19nov04: correct bad earth spin calculation
       7 jun 2005: vectorize to make faster for quantity calculations.
       20 Mar 2007: CH updated documentation for chdoppler and
       created this version, ugdoppler, which uses the locally-derived lst
       (from ilst.pro).
       5apr2011: updated documentation, tested with tst.ugdopp.idl and
       tst1.ugdopp.ilprc
    """

    dtor = np.pi/180.

    #------------------ORBITAL SECTION-------------------------
    try:
        nin = len(ra)
    except:
        ra, dec, julday = np.array([ra]), np.array([dec]), np.array([julday])
        nin = 1

    #GET THE COMPONENTS OF RA AND DEC, 2000u EPOCH
    rasource=ra*15.*dtor
    decsource=dec*dtor

    xxsource = np.zeros((3, nin))
    xxsource[0, :] = np.cos(decsource) * np.cos(rasource)
    xxsource[1, :] = np.cos(decsource) * np.sin(rasource)
    xxsource[2, :] = np.sin(decsource)
    pvorbit_helio= np.zeros( nin)
    pvorbit_bary= np.zeros( nin)
    pvlsr= np.zeros( nin)


    # GET THE EARTH VELOCITY WRT THE SUN CENTER
    # THEN MULTIPLY BY SOURCE TO GET PROJECTED VELOCITY
    # OF EARTH CENTER WRT SUN TO THE SOURCE
    for NR in range(nin):
        vvorbit, velb = baryvel(julday[NR], 2000.)
        pvorbit_helio[NR]= np.sum(vvorbit* xxsource[:,NR])
        pvorbit_bary[NR]= np.sum(velb* xxsource[:,NR])


    #-----------------------LSR SECTION-------------------------
    # THE STANDARD LSR IS DEFINED AS FOLLOWS: THE SUN MOVES AT 20.0 KM/S
    # TOWARD RA=18.0H, DEC=30.0 DEG IN 1900 EPOCH COORDS
    # using PRECESS, this works out to ra=18.063955 dec=30.004661 in 2000 coords.
    ralsr_rad= 2.*np.pi*18./24.
    declsr_rad= dtor*30.
    ralsr_rad, declsr_rad = precess(ralsr_rad, declsr_rad, 1900., 2000.,radian=True)

    #FIND THE COMPONENTS OF THE VELOCITY OF THE SUN WRT THE LSR FRAME
    xxlsr = np.zeros(3)
    xxlsr[0] = np.cos(declsr_rad) * np.cos(np.pi+ralsr_rad) #additional pi because Python and IDL just...
    xxlsr[1] = np.cos(declsr_rad) * np.sin(ralsr_rad)
    xxlsr[2] = np.sin(declsr_rad)
    vvlsr = 20.*xxlsr

    #PROJECTED VELOCITY OF THE SUN WRT LSR TO THE SOURCE
    for NR in range(nin): pvlsr[NR]=np.sum(vvlsr*xxsource[:, NR])


    #---------------------EARTH SPIN SECTION------------------------
    #NOTE: THE ORIGINAL VERSION WAS FLAWED. WE comment out those bad statements...

    # LEUSCHNER LAB COORDINATES...
    northlat = 37.91934
    westlong = 122.15385
    obspos_deg= [ northlat, westlong]


    # COORDS FROM NLAT, WLONG INPUT...
    if nlat and wlong:
        obspos_deg= [nlat, wlong]

    # GET THE LATITUDE...
    lat= obspos_deg[0]

    # GET THE LST
    lst_mean= rlab.getLST(juldate=julday, lon=-obspos_deg[1])

    # MODIFIED EARTH SPIN FROM GREEN PAGE 270
    pvspin= -0.465* np.cos(dtor* lat) * np.cos( decsource) * np.sin(( lst_mean- ra)* 15.*dtor)


    #---------------------NOW PUT IT ALL TOGETHER------------------

    vtotal= np.zeros((4, nin))
    vtotal[ 0,:]= -pvspin
    vtotal[ 1,:]= -pvspin- pvorbit_helio
    vtotal[ 2,:]= -pvspin- pvorbit_bary
    vtotal[ 3,:]= -pvspin- pvorbit_bary- pvlsr

    if light: vtotal=vtotal/(2.99792458e5)
    return vtotal

