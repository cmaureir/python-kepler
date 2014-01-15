#!/usr/bin/env python

import numpy as np
from const import *

class KeplerianOrbit:
    def __init__(self, r, v, m):

        self.r = r
        self.v = v
        self.m = float(m)

        # orbital elements
        self.j          = np.array([0,0,0])
        self.e          = np.array([0,0,0])
        self.a_vec      = np.array([0,0,0])
        self.b_vec      = np.array([0,0,0])
        self.ecc        = 0
        self.a          = 0
        self.b          = 0
        self.w          = 0
        self.m_anomaly0 = 0
        self.e_anomaly0 = 0
        self.e_anomaly  = 0
        self.m_anomaly  = 0
        self.b_mag      = 0
        self.r_mag      = 0

    def print_log(self,time):
        print time, self.r[0], self.r[1], self.r[2],\
                    self.v[0], self.v[1], self.v[2],\
                    self.e_anomaly, self.m_anomaly

    def print_oe(self):
        print "j_vec=(%s, %s, %s) e_vec=(%s, %s, %s) e=%s a=%s b=%s w=%s"\
            % (self.j[0], self.j[1], self.j[2], self.e[0], self.e[1], self.e[2],\
            self.ecc, self.a, self.b, self.w)

    def calculate_orbital_elements(self):

        mu = G * self.m

        # r magnitude
        self.r_mag = np.linalg.norm(self.r)

        # Angular momentum
        # j = r x v
        self.j = np.cross(self.r, self.v)

        # Runge-Lenz-vector
        # e = { (v x j) / (G * m) }  - { r / |r| }
        self.e = np.cross(self.v, self.j)/mu - self.r/self.r_mag

        # Eccentricity
        self.ecc = np.linalg.norm(self.e)

        # Semi-major axis
        # a = ( j * j ) / (G * m * | 1 - ecc^2 | )
        self.a = self.j.dot(self.j)/(mu * np.fabs(1-self.ecc**2))

        # Frequency of the orbit
        # w = np.sqrt(( G * m )/ a^3 )
        self.w = np.sqrt((G * self.m)/ self.a**3)

        # Semi-major vector
        # a[] = a * e[]/ecc
        self.a_vec = self.a * self.e/self.ecc

        # Semi-minor axis
        self.b = self.a * np.sqrt(np.fabs(1 - self.ecc**2))

        # Semi-minor vector
        # b[] = (j x e) * b/b_mag
        jxe = np.cross(self.j,self.e)
        self.b_vec = jxe * self.b/np.linalg.norm(jxe)


    def get_elliptical_pos_vel(self, dt):

        if dt == 0:
            # Calculate cosine of Eccentric Anomaly
            self.e_anomaly0 = (self.a - self.r_mag) / (self.ecc * self.a)

            # Fixing cosine argument
            if self.e_anomaly0 >= 1.0:
                self.e_anomaly0 = 0.0
            elif self.e_anomaly0 <= -1.0:
                self.e_anomaly0 = np.pi
            else:
                self.e_anomaly0 = np.arccos(self.e_anomaly0)

            # Calculate Mean Anomaly
            self.m_anomaly0  = (self.e_anomaly0 - self.ecc * np.sin(self.e_anomaly0))

        self.m_anomaly = self.m_anomaly0 + dt * self.w

        # Adjusting M anomaly to be < 2 * np.pi
        if self.m_anomaly >= 2 * np.pi:
            self.m_anomaly = np.fmod(self.m_anomaly, 2 * np.pi)

        # Solving the Kepler equation for elliptical orbits
        e_anomaly_new = 0
        if self.ecc > 0.8:
            e_anomaly_new = np.pi
        else:
            e_anomaly_new = self.m_anomaly

        d = 1e4
        iteration = 0
        while np.fabs(d) > DEL_E:
            d = e_anomaly_new - self.ecc * np.sin(e_anomaly_new) - self.m_anomaly
            if iteration - 1 >= KEPLER_ITE:
                break
            e_anomaly_new -= d / (1.0 - self.ecc * np.cos(e_anomaly_new))
            iteration += 1

        self.e_anomaly = e_anomaly_new

        cos_e = np.cos(self.e_anomaly)
        sin_e = np.sin(self.e_anomaly)

        r_const = cos_e - self.ecc
        v_const = self.w / (1.0 - self.ecc * cos_e)

        # Based on Loeckmann and Baumgardt simulations criteria
        # This work better with 0.99 < e < 1 and |E| < 1e-3
        if self.ecc > 0.99:
            if self.e_anomaly > 2.0 * np.pi - 1e-3:
                e_tmp -= 2.0 * np.pi
            else:
                e_tmp = self.e_anomaly

            if e_tmp < 1e-3:
                e_tmp *= e_tmp
                j_mag = self.j.dot(self.j)
                ecc_const = j_mag/(self.m0 * self.a * (1 + self.ecc))
                cos_const = -0.5 * e_tmp * (1 - e_tmp / 12.0 * (1 - e_tmp / 30.0))

                r_const = ecc_const + cos_const
                v_const = self.w / (ecc_const - self.ecc * cos_const)

        # Update position and velocity
        self.r = self.a_vec * r_const + self.b_vec * sin_e
        self.v = (-self.a_vec * sin_e + self.b_vec * cos_e) * v_const


    def get_hyperbolic_pos_vel(self, dt):

        if dt == 0:
            rdotv = self.r.dot(self.v)
            # calculate eccentric anomaly e at t+dt
            self.e_anomaly0 = (self.a + self.r_mag) / (self.ecc * self.a)

            if self.e_anomaly0 < 1.0:
                self.e_anomaly0 = 0.0
            elif rdotv < 0:
                self.e_anomaly0 = -np.arccosh(self.e_anomaly0)
            else:
                self.e_anomaly0 = np.arccosh(self.e_anomaly0)

            self.m_anomaly0 = self.ecc * np.sinh(self.e_anomaly0) - self.e_anomaly0

        self.m_anomaly = self.m_anomaly0 + dt * self.w

        #if self.m_anomaly >= 2 * np.pi:
        #    self.m_anomaly = np.fmod(self.m_anomaly, 2 * np.pi)

        # Solving Kepler's equation for Hyperbolic/Parabolic orbits
        e_anomaly_new = self.m_anomaly / np.fabs(1. - self.ecc)

        is_negative = False
        if self.m_anomaly < 0:
            self.m_anomaly = -self.m_anomaly
            is_negative = True

        #thresh = DEL_E_HYP * self.m_anomaly
        thresh = DEL_E_HYP * np.fabs(1 - self.ecc)

        # cubic term is dominant
        if e_anomaly_new * e_anomaly_new > 6. * np.fabs(1. - self.ecc):
            if self.m_anomaly < np.pi:
                e_anomaly_new = np.exp(np.log(6. * self.m_anomaly)/3.0)
            # hyperbolic w/ 5th & higher-order terms predominant
            else:
                e_anomaly_new = np.arcsinh(self.m_anomaly / self.ecc)


        e_anomaly_new = np.log(2.*self.m_anomaly/self.ecc+1.8) # taken from Burkardt & Danby, CeMec 31 (1983), 317-328
        e_anomaly_new_abs = np.fabs(e_anomaly_new)
        d = ((self.ecc * np.sinh(e_anomaly_new)) - e_anomaly_new) - self.m_anomaly

        iteration = 0
        while np.fabs(d) > thresh:
            iteration += 1
            if e_anomaly_new_abs < .72 and self.ecc < 1.1:
                # [e * np.sinh(E) - E] / E << 1, and/or e * np.cosh(E) - 1 << 1
                # so don't calculate it directly
                e_anomaly_new2 = e_anomaly_new * e_anomaly_new

                # relative error when omitting nth order term needs to be smaller than resolution 1.e-15:
                # .5 * E^2 > 1e15 * E^n/n!, i.e. E < (n!/2e15)^(1/(n-2))
                # n = 16: E < .72, n = 10: E < .08

                if e_anomaly_new_abs > .08:
                    e_anomaly_new -= d / ((self.ecc - 1) * np.cosh(e_anomaly_new) + \
                            (((((((_1_16_15 * e_anomaly_new2 + 1.) \
                            * _1_14_13 * e_anomaly_new2 + 1.) \
                            * _1_12_11 * e_anomaly_new2 + 1.) \
                            * _1_10_9 * e_anomaly_new2 + 1.) \
                            * _1_8_7 * e_anomaly_new2 + 1.) \
                            * _1_6_5 * e_anomaly_new2 + 1.) \
                            * _1_4_3 * e_anomaly_new2 + 1.) \
                            * .5 *e_anomaly_new2)
                else:
                    e_anomaly_new -= d / ((self.ecc - 1) * np.cosh(e_anomaly_new) + \
                            ((((_1_10_9 * e_anomaly_new2 + 1.) \
                            * _1_8_7 * e_anomaly_new2 + 1.) \
                            * _1_6_5 * e_anomaly_new2 + 1.) \
                            * _1_4_3 * e_anomaly_new2 + 1.) \
                            * .5 *e_anomaly_new2)

                e_anomaly_new2 = e_anomaly_new * e_anomaly_new
                e_anomaly_new_abs = np.fabs(e_anomaly_new)

                if e_anomaly_new_abs > .08:
                    d = ((ecc - 1) * np.sinh(e_anomaly_new) +
                            (((((((_1_17_16 * e_anomaly_new2 + 1.)
                            * _1_15_14 * e_anomaly_new2 + 1.)
                            * _1_13_12 * e_anomaly_new2 + 1.)
                            * _1_11_10 * e_anomaly_new2 + 1.)
                            * _1_9_8 * e_anomaly_new2 + 1.)
                            * _1_7_6 * e_anomaly_new2 + 1.)
                            * .05 * e_anomaly_new2 + 1.)
                            * _1_3_2 * e_anomaly_new2 * e_anomaly_new) - self.m_anomaly
                else:
                    d = ((self.ecc - 1) * np.sinh(e_anomaly_new) +
                            ((((_1_11_10 * e_anomaly_new2 + 1.)
                            * _1_9_8 * e_anomaly_new2 + 1.)
                            * _1_7_6 * e_anomaly_new2 + 1.)
                            * .05 * e_anomaly_new2 + 1.)
                            * _1_3_2 * e_anomaly_new2 * e_anomaly_new) - self.m_anomaly
            else:
                e_anomaly_new -= d / (self.ecc * np.cosh(e_anomaly_new) - 1.0)
                d = ((self.ecc * np.sinh(e_anomaly_new)) - e_anomaly_new) - self.m_anomaly

            if iteration - 1 >= KEPLER_ITE:
                break

        self.e_anomaly = e_anomaly_new
        if is_negative:
            self.e_anomaly *= -1

        cos_e = np.cosh(self.e_anomaly)
        sin_e = np.sinh(self.e_anomaly)
        v_const = self.w / (self.ecc * cos_e - 1.)

        # New position and velocity
        self.r = self.a_vec * (self.ecc - cos_e) + self.b_vec * sin_e
        self.v = (-self.a_vec * sin_e + self.b_vec * cos_e) * v_const

