"""

Helper function which suggests convergent TEMPO parameters.

Choose your own then plot to visualise convergence.

(Time indeoendent)

"""

import os
import sys
sys.path.insert(0,'..')

import oqupy
import numpy as np
import matplotlib.pyplot as plt
import timeit
import math


start = timeit.default_timer()

# Define operators
sigma_x = oqupy.operators.sigma("x")
sigma_y = oqupy.operators.sigma("y")
sigma_z = oqupy.operators.sigma("z")
sigma_plus = oqupy.operators.sigma("+")
sigma_minus = oqupy.operators.sigma("-")

up_density_matrix = oqupy.operators.spin_dm("z+")
down_density_matrix = oqupy.operators.spin_dm("z-")
mixed_density_matrix = oqupy.operators.spin_dm("mixed")


# Pick System Parameters
Omega = 1.0
omega_cutoff = 5.0
alpha = 0.0001
# Temperature in kelvin
T = 0 # K
# Temperature in units of 1/(ps kB)
# 1 K = 0.1309 1/(ps kB)
temperature = T*0.1309 # 1/(ps kB)
dur = 2000

# initial state = up_density_matrix or mixed_density_matrix
init_st = up_density_matrix


# Specify system
system = oqupy.System(Omega * sigma_z)

# Specify correlations
correlations = oqupy.PowerLawSD(alpha=alpha,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='exponential',
                                temperature=(temperature))

bath = oqupy.Bath(sigma_x, correlations)



# Find parameters!! ###########################
parameters = oqupy.guess_tempo_parameters(system=system,
                                          bath=bath,
                                          start_time=0.0,
                                          end_time=dur,
                                          tolerance=0.01)
print(parameters)

fig, ax = plt.subplots(1,1)
oqupy.helpers.plot_correlations_with_parameters(bath.correlations, parameters, ax=ax)


# Pick my params
dt = 0.05
dkmax = 300
epsrel = 10**(-6)

tempo_parameters = oqupy.TempoParameters(dt=dt, dkmax=dkmax, epsrel=epsrel, name="my rough parameters")
print(tempo_parameters)

fig, ax = plt.subplots(1,1)
oqupy.helpers.plot_correlations_with_parameters(bath.correlations, tempo_parameters, ax=ax)

desired_folder = 'C:/Users/sony/Desktop/Project/Plots/converg_plots'
file_name = 'alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Spect&Om.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()
# tempo = oqupy.Tempo(system=system,
#                       bath=bath,
#                       parameters=tempo_parameters,
#                       initial_state=up_density_matrix,
#                       start_time=0.0)

# tempo.compute(end_time=5.0)

# dynamics_2 = tempo.get_dynamics()
# plt.plot(*dynamics_2.expectations(sigma_z, real=True), label=r'$\alpha=0.3$')
# plt.xlabel(r'$t\,\Omega$')
# plt.ylabel(r'$\langle\sigma_z\rangle$')
# plt.legend()

# tempo.compute(end_time=15.0)

# dynamics_2 = tempo.get_dynamics()
# plt.plot(*dynamics_2.expectations(sigma_z, real=True), label=r'$\alpha=0.3$')
# plt.xlabel(r'$t\,\Omega$')
# plt.ylabel(r'$\langle\sigma_z\rangle$')
# plt.legend()
#####################


# tempo_parameters = oqupy.TempoParameters(dt=0.1, dkmax=100, epsrel=10**(-4))

# dynamics = oqupy.tempo_compute(system=system,
#                                 bath=bath,
#                                 initial_state=init_st,
#                                 start_time=0.0,
#                                 end_time=10.0,
#                                 parameters=tempo_parameters)
# #and plot the result
# t, s_z = dynamics.expectations(sigma_z, real=True)

# plt.plot(t, s_z, label=r'TEMPO $\alpha={0}$ '.format(alpha))
# plt.title('Decay from Up Density Matrix State at T = {0} K'.format(T)) 
# plt.xlabel(r'$t\,\Omega$')
# plt.ylabel(r'$\langle\sigma_z\rangle$')
# plt.legend()

stop = timeit.default_timer()
print('Time: ', stop - start)