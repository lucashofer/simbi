# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 14:24:14 2022

@author: hofer
"""

import numpy as np
import matplotlib.pyplot as plt
import simbi

a0 = 5.29177210903 * 10 ** (-11)
m = 1.4192261 * 10 ** (-25)
wx = 2 * np.pi * 300
wy = 2 * np.pi * 30
wz = 2 * np.pi * 100
# mu = 8.60674437635897e-30
mu_0 = 1.1e-30
T_0 = 1.1e-06
t = .0001
a = 90 * a0
trap_freqs = (wx, wy, wz)

ntot_gt = 100000
cf_gt = 0.005
nbec = cf_gt * ntot_gt


mu, T = simbi.mu_temperature(trap_freqs, m, a, ntot_gt, cf_gt)
atom_numbers = simbi.bimodal_atom_numbers(trap_freqs, m, T, a, mu)


all_radii = simbi.get_radii(trap_freqs, m, T, a, mu, t)
init_sigmas, t_sigmas, init_tf_radii, expans_scalars, tf_radii = all_radii

(x,y), coords_2d, da, clength = simbi.get_coordinates(t_sigmas[:2], 6, multi=False)

densities_2d = simbi.bimodal_density(coords_2d, trap_freqs, m, T, a, mu, t)
#%%
plabels = ['Combined', 'Excited', 'BEC']
slicev = int(clength / 2)

fig, axes = plt.subplots(1, 3)
for ax, label, density in zip(axes, plabels, densities_2d):
    ax.imshow(density)
    ax.set_title(label)
plt.show()


plt.figure()
for label, density in zip(plabels, densities_2d):
    plt.plot(x, density[slicev], label=label)
plt.legend()
plt.show()


# Nbec = 5000
# cf = .1
# Ntot = Nbec / cf
# Nex = (1 - cf) * Ntot


# mu_calc = bec.get_chemical_potential(*trap_freqs, m, a, Nbec)
# T = ts.get_temperature(*trap_freqs, -mu_calc, Nex)

# from bimodal_cloud.thermal_cloud.enhanced_bose import ThermalCloud, TemperatureSolver
# from bimodal_cloud.thomas_fermi.thomas_fermi_bec import ThomasFermiBEC
