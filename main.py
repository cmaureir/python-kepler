#!/usr/bin/env python

from KeplerianOrbit import *

def get_pos_vel_from_orb_elem(a, e, m):

    mu = G * m
    if e < 1:
        r = a*(1+e)
        v = np.sqrt((mu/a) * (1-e)/(1+e))
    else:
        r = a*(1-e)
        v = np.sqrt((mu/a) * (1+e)/(1-e))


    # central object
    r0 = np.array([0,0,0])
    v0 = np.array([0,0,0])
    m0 = m

    # star
    r1 = np.array([r,0,0])
    v1 = np.array([0,v,0])
    m1 = m/1000

    return r1-r0, v1-v0

if __name__ == "__main__":

    # Integration time
    MAX_TIME = 1
    # Integration step
    step = 0.001
    # Initialization
    time = 0

    ## Example case for hyperbola
    a = -1.0
    e = 1.5

    ## Example case for ellipse
    #a = 1
    #e = 0.5

    m = 1000
    r,v = get_pos_vel_from_orb_elem(a, e, m)
    print r,v, m

    kp = KeplerianOrbit(r,v,m)

    while time < MAX_TIME:
        kp.calculate_orbital_elements()
        kp.print_oe()
        raw_input()
        if e < 1:
            kp.get_elliptical_pos_vel(time)
            print kp.r
        else:
            kp.get_hyperbolic_pos_vel(time)
        kp.print_log(time)
        time += step
