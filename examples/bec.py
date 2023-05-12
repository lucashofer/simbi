# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 01:48:34 2022

@author: hofer
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import physical_constants
import simbi


'''Some sensible experimental parameters'''
a0 = physical_constants['Bohr radius'][0]
amu = physical_constants['atomic mass constant'][0]
a = 100 * a0
m = 85.46 * amu
wx = 2 * np.pi * 350
wy = 2 * np.pi * 300
wz = 2 * np.pi * 30

nbec = 10000
t = .001
trap_freqs = (wx, wy, wz)
#initialize objects

mu = simbi.chemical_potential(trap_freqs, m, a, nbec)
_, _, tf_radii = simbi.get_radii(trap_freqs, m, a, mu, t)

#setup coordinates for calculating densities in one, two and three dimensions
_, coord_maps, diff_elements, clength = simbi.get_coordinates(tf_radii, 2, multi=True)
coords_3d, coords_2d, coords_1d = coord_maps
dv, da, dx = diff_elements

#get density and integrated densities
density_3d = simbi.bec_density(coords_3d, trap_freqs, m, a, mu, t)
int_density_2d = simbi.bec_density(coords_2d, trap_freqs, m, a, mu, t)
int_density_1d = simbi.bec_density(coords_1d, trap_freqs, m, a, mu, t)

#make sure analytic atom number matches atom number calculated from densities
print('Analytic Atom Number: ', simbi.bec_atom_number(trap_freqs, m, a, mu))
print('3D Integrated Atom Number: ', simbi.integrate_array(density_3d, dv))
print('2D Integrated Atom Number: ', simbi.integrate_array(int_density_2d, da))
print('1D Integrated Atom Number: ', simbi.integrate_array(int_density_1d, dx))

#plot the results
plt.figure()
plt.plot(density_3d[:, int(clength/2), int(clength/2)])
plt.title('Cut Through')

plt.figure()
plt.pcolormesh(*coords_2d, int_density_2d)
plt.title('Column Density')

plt.figure()
plt.plot(int_density_1d)
plt.title('Line Density')