"""
System plot with derived master equation in linblad form
"""

import sys
sys.path.insert(0,'..')

import oqupy
import numpy as np
import matplotlib.pyplot as plt
import timeit


start = timeit.default_timer()

sigma_x = oqupy.operators.sigma("x")
sigma_y = oqupy.operators.sigma("y")
sigma_z = oqupy.operators.sigma("z")
sigma_plus = oqupy.operators.sigma("+")
sigma_minus = oqupy.operators.sigma("-")

up_density_matrix = oqupy.operators.spin_dm("z+")
down_density_matrix = oqupy.operators.spin_dm("z-")

gamma_1 = 1
gamma_2 = 1

        

Gamma = [gamma_1, gamma_1, gamma_2, gamma_2]


lindblad_operators  = [2*sigma_minus, 2*sigma_plus, 2*sigma_plus, 2*sigma_minus]
    
# # Same as in notes
# lindblad_operators = [sigma_minus, sigma_plus]


systemD = oqupy.System(
        sigma_z,
        gammas=Gamma,
        lindblad_operators=lindblad_operators)

# Correct?
initial_state = 0.5*(up_density_matrix + down_density_matrix)

dynamics = oqupy.compute_dynamics(system = systemD, initial_state=sigma_z,\
                                  start_time = 0,dt=0.01,num_steps=1400)


t, s_z = dynamics.expectations(sigma_z, real=True)

plt.plot(t, s_z, label=r'$\alpha=0.3$')
plt.xlabel(r'$t\,\Omega$')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()



# stop = timeit.default_timer()
# print('Time: ', stop - start)