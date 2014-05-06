"""
Some coordinate conversions
Placed here for laziness
Aaron Tran
Astro 121, Spring 2014
"""

import numpy as np

def rad2deg(pt):
    """Convert radians to degrees"""
    return pt * 180./np.pi

def deg2rad(pt):
    """Convert degrees to radians"""
    return pt * np.pi/180

def hr2rad(pt):
    """Convert hours to radians"""
    return pt * np.pi/12

def rad2hr(pt):
    """Convert radians to hours"""
    return pt * 12/np.pi

def dublinJD2JD(t):
    """Dublin JD to unmodified JD (days, both)"""
    return t + 2415020

def JD2dublinJD(t):
    """Unmodified JD to Dublin JD (days, both)"""
    return t - 2415020


if __name__ == '__main__':
    pass
