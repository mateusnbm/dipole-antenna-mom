#
# main.py
#


import math
import numpy as np
import matplotlib.pyplot as plot


'''
Constants.
'''

PI = 3.14
lightspeed = 3*math.pow(10, 8)
electric_permittivity_vacuum = 8.85*math.pow(10, -12)
magnetic_permeability_vacuum = 4*PI*math.pow(10, -7)


'''
Parameters.
'''

frequency  = 300*math.pow(10, 6)
angular_frequency = 2*PI*frequency

wavelength = float(lightspeed)/frequency
wavenumber = (2*PI)/wavelength

antenna_length = wavelength#*(0.5)
antenna_radius = wavelength*math.pow(10, -5)

segments_count  = 29 #19  #7 #3 #19
powered_segment = 15 #10  #4 #2 #10
half_base_width = float(antenna_length)/(segments_count+1)

electric_permittivity = electric_permittivity_vacuum
magnetic_permeability = magnetic_permeability_vacuum


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
Solve linear system (find currents).
'''

I = np.linalg.solve(Z, V)


'''
Plot the current values along the z axis.
'''

I_absolute = np.absolute(I)
I_mV = [(i*1000) for i in I_absolute]

#x = [-0.25, -0.15, 0, 0.15, 0.25] #[-0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25]

I_mV = [0] + I_mV + [0]
segments_count += 2

x_neg = [ ( ((n*half_base_width)-(antenna_length/2)) / wavelength)  for n in range((segments_count)//2)]
x_pos = [-x for x in x_neg]
x_rev = [x for x in reversed(x_pos)]

x = x_neg + [0] + x_rev

for i in range(len(x)):
    print(str(x[i]) + " " + str(I_mV[i]))

plot.figure()
plot.plot(x, I_mV, "s-")
plot.show()


'''
Plot antenna input impedance.
'''

# ...
