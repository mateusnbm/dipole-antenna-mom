#
# implementation.py
#


import math
import numpy as np
import scipy.constants
import matplotlib.pyplot as plot


'''
Constants.
'''

PI = scipy.constants.pi
lightspeed = scipy.constants.c
electric_permittivity_vacuum = scipy.constants.epsilon_0
magnetic_permeability_vacuum = scipy.constants.mu_0


'''
Solve the linear system [Z][I]=[V] to find the currents along the z axis.
'''

def solve_dipole_antenna_mom(
        frequency,
        antenna_length,
        antenna_radius,
        segments_count,
        powered_segment
    ):

    '''
    Inferred parameters.
    '''

    angular_frequency = 2*PI*frequency

    wavelength = float(lightspeed)/frequency
    wavenumber = (2*PI)/wavelength

    half_base_width = float(antenna_length)/(segments_count+1)

    electric_permittivity = electric_permittivity_vacuum
    magnetic_permeability = magnetic_permeability_vacuum

    '''
    Auxiliary function.
    '''

    def phi(m, n):

        if m == n:

            a = (1.0/(2*PI*half_base_width))
            b = math.log(half_base_width/antenna_radius)
            c = a * b

            d = (1/(4*PI))
            e = np.complex(0, wavenumber)
            f = d * e

            return c - f

        else:

            z_m = (m*half_base_width) - (antenna_length/2)
            z_n = (n*half_base_width) - (antenna_length/2)

            a = math.pow((z_m - z_n), 2)
            b = math.pow(antenna_radius, 2)
            c = math.sqrt(a + b)

            d = np.complex(0, -(wavenumber*c))
            e = np.exp(d)

            f = 4*PI*c

            return e / f

    '''
    Initialize matrices.
    '''

    Z = np.zeros((segments_count, segments_count), dtype=np.complex)
    V = np.zeros((segments_count, 1), dtype=np.complex)

    '''
    Calculate the voltages column matrix (V).
    '''

    V[powered_segment-1][0] = np.complex(0, -(angular_frequency*electric_permittivity))

    '''
    Calculate the impedance matrix.
    '''

    for m in range(segments_count):

        for n in range(segments_count):

            a = phi((m-0.5), (n-0.5))
            b = phi((m-0.5), (n+0.5))
            c = phi((m+0.5), (n-0.5))
            d = phi((m+0.5), (n+0.5))

            k = math.pow(wavenumber, 2)
            A_mn = math.pow(half_base_width, 2) * phi(m, n)
            O_mn = a - b - c + d

            Z[m][n] = (k*A_mn) - O_mn

    '''
    Solve the linear system [Z][I]=[V] (find currents).
    '''

    I = np.linalg.solve(Z, V)

    return I
