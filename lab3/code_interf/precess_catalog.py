"""
Short script to precess RA, dec values for least-squares fits etc
Astro 121, Spring 2014
Aaron Tran

"""

import ephem

import catalog
import ay121_ephem

def main():
    """Precess coordinates for point sources of interest

    Specifically, M17, Orion, and Virgo A
    """
    obs = ay121_ephem.obs_berkeley()
    obs_day = '2014/04/07 00:00'  # PDT
    init_time = ephem.Date(ephem.Date(obs_day) + 7*ephem.hour)
    obs.date = init_time

    targets = [catalog.orion(), catalog.M17(), catalog.virgoA()]
    
    for src in targets:
        src.compute(obs)
        print src.name
        print 'J2000 ra, dec: ', src._ra, ',', src._dec
        print 'Current ra, dec: ', src.ra, ',', src.dec, '\n'


if __name__ == "__main__":
    main()

