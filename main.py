#
# main.py
#


import os
import math
import implementation
import numpy as np
import scipy.constants
import matplotlib.pyplot as plot


frequency  = 300*math.pow(10, 6)
wavelength = float(scipy.constants.c)/frequency


'''
Plot the current distribution of a half-wavelength dipole antenna
along the z axis using 3, 7 and 19 triangular base functions. Also,
print the input impedance of the configurations.
'''

plots = []
impedances = []

segments_count = [3, 7, 19]
antenna_length = wavelength*(0.5)
antenna_radius = wavelength*math.pow(10, -4)

for count in segments_count:

    powered_segment = (count//2)+1
    half_base_width = float(antenna_length)/(count+1)

    I = implementation.solve_dipole_antenna_mom(
        frequency=frequency,
        antenna_length=antenna_length,
        antenna_radius=antenna_radius,
        segments_count=count,
        powered_segment=powered_segment
    )

    I_in = I[powered_segment-1]
    Z_in = 1/I_in

    I_absolute = np.absolute(I)
    I_mV = [(i*1000) for i in I_absolute]
    I_mV = [0] + I_mV + [0]

    x_neg = [ ( ((n*half_base_width)-(antenna_length/2)) / wavelength)  for n in range((count+2)//2)]
    x_pos = [-x for x in x_neg]
    x_rev = [x for x in reversed(x_pos)]
    x = x_neg + [0] + x_rev

    plots.append([x, I_mV])
    impedances.append(Z_in)

plot.figure()
plot.plot(plots[0][0], plots[0][1], "g-+")
plot.plot(plots[1][0], plots[1][1], "b-x")
plot.plot(plots[2][0], plots[2][1], "r-")
plot.legend(["N=3", "N=7", "N=19"], loc="upper right")
plot.title("Current Distribution (L = λ/2)")
plot.xlabel("Z/λ")
plot.ylabel("|I| (mA)")
plot.savefig("graphs/half-wavelength-currents.pdf", bbox_inches="tight")

print("")
print("Half-wavelength dipole antenna input impedance:")
for i in range(len(segments_count)):
    print(format(segments_count[i], '02d') + " " + str(impedances[i]))

'''
Plot the inpute impedance active and reactive components for a
half-wavelength dipole antenna using 3..49 triangular base functions.
'''

actives = []
reactives = []

segments_count = [x for x in range(3, 50) if (x%2) != 0]
antenna_length = wavelength*(0.5)
antenna_radius = wavelength*math.pow(10, -4)

for count in segments_count:

    powered_segment = (count//2)+1
    half_base_width = float(antenna_length)/(count+1)

    I = implementation.solve_dipole_antenna_mom(
        frequency=frequency,
        antenna_length=antenna_length,
        antenna_radius=antenna_radius,
        segments_count=count,
        powered_segment=powered_segment
    )

    I_in = I[powered_segment-1]
    Z_in = 1/I_in

    actives.append([count, Z_in.real])
    reactives.append([count, Z_in.imag])

plot.figure()
plot.plot([x[0] for x in actives], [x[1] for x in actives], "g-")
plot.plot([x[0] for x in reactives], [x[1] for x in reactives], "b-")
plot.legend(["R (Ω)", "X (Ω)"], loc="center right")
plot.title("Input Impedance (L = λ/2)")
plot.xlabel("Z/λ")
plot.savefig("graphs/half-wavelength-impedances.pdf", bbox_inches="tight")


'''
Plot the current distribution of a wavelength dipole antenna along
the z axis using 3, 7 and 19 triangular base functions. Also, print
the input impedance of the configurations.
'''

plots = []
impedances = []

segments_count = [3, 7, 19]
antenna_length = wavelength
antenna_radius = wavelength*math.pow(10, -4)

for count in segments_count:

    powered_segment = (count//2)+1
    half_base_width = float(antenna_length)/(count+1)

    I = implementation.solve_dipole_antenna_mom(
        frequency=frequency,
        antenna_length=antenna_length,
        antenna_radius=antenna_radius,
        segments_count=count,
        powered_segment=powered_segment
    )

    I_in = I[powered_segment-1]
    Z_in = 1/I_in

    I_absolute = np.absolute(I)
    I_mV = [(i*1000) for i in I_absolute]
    I_mV = [0] + I_mV + [0]

    x_neg = [ ( ((n*half_base_width)-(antenna_length/2)) / wavelength)  for n in range((count+2)//2)]
    x_pos = [-x for x in x_neg]
    x_rev = [x for x in reversed(x_pos)]
    x = x_neg + [0] + x_rev

    plots.append([x, I_mV])
    impedances.append(Z_in)

plot.figure()
plot.plot(plots[0][0], plots[0][1], "g-+")
plot.plot(plots[1][0], plots[1][1], "b-x")
plot.plot(plots[2][0], plots[2][1], "r-")
plot.legend(["N=3", "N=7", "N=19"], loc="upper right")
plot.title("Current Distribution (L = λ)")
plot.xlabel("Z/λ")
plot.ylabel("|I| (mA)")
plot.savefig("graphs/wavelength-currents.pdf", bbox_inches="tight")

print("")
print("Wavelength dipole antenna input impedance:")
for i in range(len(segments_count)):
    print(format(segments_count[i], '02d') + " " + str(impedances[i]))
print("")


'''
Plot the inpute impedance active and reactive components for a
wavelength dipole antenna using 3..49 triangular base functions.
'''

actives = []
reactives = []

segments_count = [x for x in range(3, 50) if (x%2) != 0]
antenna_length = wavelength
antenna_radius = wavelength*math.pow(10, -4)

for count in segments_count:

    powered_segment = (count//2)+1
    half_base_width = float(antenna_length)/(count+1)

    I = implementation.solve_dipole_antenna_mom(
        frequency=frequency,
        antenna_length=antenna_length,
        antenna_radius=antenna_radius,
        segments_count=count,
        powered_segment=powered_segment
    )

    I_in = I[powered_segment-1]
    Z_in = 1/I_in

    actives.append([count, Z_in.real])
    reactives.append([count, Z_in.imag])

plot.figure()
plot.plot([x[0] for x in actives], [x[1] for x in actives], "g-")
plot.plot([x[0] for x in reactives], [x[1] for x in reactives], "b-")
plot.legend(["R (Ω)", "X (Ω)"], loc="center right")
plot.title("Input Impedance (L = λ/2)")
plot.xlabel("Z/λ")
plot.savefig("graphs/wavelength-impedances.pdf", bbox_inches="tight")


'''
Delete bytecode files and the cache folder.
'''

os.system("find . -type f -name \"*.pyc\" -delete")
os.system("rm -rf __pycache__")
