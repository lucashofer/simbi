# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 14:24:14 2022

@author: hofer
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import physical_constants
import simbi


def print_check_densities(all_densities, atom_numbers, diff_elements):
    tlabels = ['Combined Atom Number', 'Thermal Atom Number', 'BEC Atom Number']
    ilabels = ['Integrated 3D', 'Integrated 2D', 'Integrated 1D']
    
    #rearrange densities by combined, thermal, bec rather than dimension
    all_densities = [[densities[i] for densities in all_densities] for i in range(3)]

    for tlabel, densities, atom_number in zip(tlabels, all_densities, atom_numbers):
        print(tlabel)
        print('Analytic: ', atom_number)
        for ilabel, density, diff_el in zip(ilabels, densities, diff_elements):
            print(ilabel, simbi.integrate_array(density, diff_el))
        print('-'*15)


'''Some sensible experimental parameters'''
a0 = physical_constants['Bohr radius'][0]
amu = physical_constants['atomic mass constant'][0]
mu = 7e-31
T = 200e-09
a = 90 * a0
m = 85.46 * amu
wx = 2 * np.pi * 350
wy = 2 * np.pi * 300
wz = 2 * np.pi * 30
trap_freqs = (wx, wy, wz)
nbec = 10000
t = .001

atom_numbers = simbi.bimodal_atom_numbers(trap_freqs, m, T, a, mu)
total_atom_number, ex_atom_number, bec_atom_number, condensed_fraction = atom_numbers

all_radii = simbi.bimodal_radii(trap_freqs, m, T, a, mu, t)
init_sigmas, t_sigmas, init_tf_radii, expans_scalars, tf_radii = all_radii

(x,y,z), coord_maps, diff_elements, clength = simbi.get_coordinates(t_sigmas, 
                                                                    6, 
                                                                    multi=True)
coords_3d, coords_2d, coords_1d = coord_maps
dv, da, dx = diff_elements

densities_3d = simbi.bimodal_density(coords_3d, trap_freqs, m, T, a, mu, t)
densities_2d = simbi.bimodal_density(coords_2d, trap_freqs, m, T, a, mu, t)
densities_1d = simbi.bimodal_density(coords_1d, trap_freqs, m, T, a, mu, t)
all_densities = [densities_3d, densities_2d, densities_1d]

#example here of how the these are actually lists containing numpy arrays for each
cd_3d, ex_3d, becd_3d = densities_3d
#%%
print_check_densities(all_densities, atom_numbers, diff_elements)

plabels = ['Combined', 'Excited', 'BEC']
slicev = int(clength / 2)

plt.figure()
plt.title('3D Cut-through')
for label, density in zip(plabels, densities_3d):
    plt.plot(x, density[slicev, slicev], label=label)
plt.legend()

fig, axes = plt.subplots(1, 3)
for ax, label, density in zip(axes, plabels, densities_2d):
    ax.imshow(density)
    ax.set_title(label)
    
plt.figure()
plt.title('1D Integrated')
for label, density in zip(plabels, densities_1d):
    plt.plot(x, density, label=label)
plt.legend()
