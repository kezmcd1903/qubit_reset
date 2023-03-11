"""
OQuPy start
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
up_density_matrix = oqupy.operators.spin_dm("z+")
down_density_matrix = oqupy.operators.spin_dm("z-")



Omega = 1.0
omega_cutoff = 5.0
alpha = 0.3

system = oqupy.System(0.5 * Omega * sigma_x)
correlations = oqupy.PowerLawSD(alpha=alpha,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='exponential')
bath = oqupy.Bath(0.5 * sigma_z, correlations)
tempo_parameters = oqupy.TempoParameters(dt=0.1, dkmax=30, epsrel=10**(-4))

dynamics = oqupy.tempo_compute(system=system,
                               bath=bath,
                               initial_state=up_density_matrix,
                               start_time=0.0,
                               end_time=15.0,
                               parameters=tempo_parameters)
#and plot the result
t, s_z = dynamics.expectations(0.5*sigma_z, real=True)

plt.plot(t, s_z, label=r'$\alpha=0.3$')
plt.xlabel(r'$t\,\Omega$')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
##################################
#print(parameters)
parameters = oqupy.guess_tempo_parameters(system=system,
                                          bath=bath,
                                          start_time=0.0,
                                          end_time=5.0,
                                          tolerance=0.01)
print(parameters)

fig, ax = plt.subplots(1,1)
oqupy.helpers.plot_correlations_with_parameters(bath.correlations, parameters, ax=ax)
############

tempo_parameters = oqupy.TempoParameters(dt=0.1, dkmax=30, epsrel=10**(-4), name="my rough parameters")
print(tempo_parameters)

fig, ax = plt.subplots(1,1)
oqupy.helpers.plot_correlations_with_parameters(bath.correlations, tempo_parameters, ax=ax)
##############

tempo = oqupy.Tempo(system=system,
                      bath=bath,
                      parameters=tempo_parameters,
                      initial_state=up_density_matrix,
                      start_time=0.0)

tempo.compute(end_time=5.0)

dynamics_2 = tempo.get_dynamics()
plt.plot(*dynamics_2.expectations(0.5*sigma_z, real=True), label=r'$\alpha=0.3$')
plt.xlabel(r'$t\,\Omega$')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()

tempo.compute(end_time=15.0)

dynamics_2 = tempo.get_dynamics()
plt.plot(*dynamics_2.expectations(0.5*sigma_z, real=True), label=r'$\alpha=0.3$')
plt.xlabel(r'$t\,\Omega$')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()

stop = timeit.default_timer()
print('Time: ', stop - start)