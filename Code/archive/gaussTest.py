"""
Time dependent system hamiltonian with process tensor
"""

import sys
import os
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

# Pick parameters
Omega = 1.0
omega_cutoff = 1
alpha = 0.001
# Temperature in kelvin
T = 0 # K
# Temperature in units of 1/(ps kB). {1 K = 0.1309 1/(ps kB)}
temperature = T*0.1309 # 1/(ps kB)
# initial state = up_density_matrix or mixed_density_matrix
init_st = up_density_matrix

# TEMPO parameters
dt=0.05
dkmax=100
epsrel=10**(-6)
# duration (ps)
dur = 4.0 # can go up to current PT length

# For Master Equation (need to change if changing dtME)
num_steps = dur*10

def gaussian_shape(t, area = 1.0, tau = 1.0, t_0 = 0.0):
    return area/(tau*np.sqrt(np.pi)) * np.exp(-(t-t_0)**2/(tau**2))
# # Spectral density ohmic
def J(t):
    return 4*alpha*gaussian_shape(t)*np.exp(-2*gaussian_shape(t)/omega_cutoff)


# Define time dependent system
# (Time dependent system gaussian_shape(t, area = np.pi/2.0, tau = 0.245))
def hamiltonian_t(t):
    return gaussian_shape(t, area = np.pi/2.0, tau = 0.245)/2.0 * sigma_z

system = oqupy.TimeDependentSystem(hamiltonian_t)
correlations = oqupy.PowerLawSD(alpha=alpha,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='gaussian',
                                temperature=temperature)
bath = oqupy.Bath(sigma_x, correlations)

tempo_parameters = oqupy.TempoParameters(dt=dt, dkmax=dkmax, epsrel=epsrel)

tempo_sys = oqupy.Tempo(system=system,
                        bath=bath,
                        initial_state=init_st,
                        start_time=0,
                        parameters=tempo_parameters)
dynamics = tempo_sys.compute(end_time=dur)
    
# Plot the result
t1, s_z = dynamics.expectations(sigma_z, real=True)



# B-E Occupancy (for temp dependence)
def N(T):
    if T == 0:
        N = 0
    else:
        N = (1/(math.exp(((1.055e-34)*(10**12))/((1.381e-23)*T))-1))
    return N

# Decay factors
#gamma1 = gamma*(N+1)
gamma_1 = 2*np.pi*(N(T)+1)#*(N+1)
# gamma2 = gamma*(N)
gamma_2 = 2*np.pi*N(T) #*N

Gamma = [ lambda t: gamma_1*J(t), lambda t: gamma_2*J(t)]

lindblad_operators = [ lambda t: sigma_minus, lambda t: sigma_plus]

    
systemD = oqupy.TimeDependentSystem(hamiltonian_t,
                                   gammas=Gamma,
                                   lindblad_operators=lindblad_operators)

# systemD = oqupy.System(
#          gaussian_shape(t, area = np.pi/2.0, tau = 0.245) * sigma_z,
#         gammas=Gamma,
#         lindblad_operators=lindblad_operators)

# up_density_matrix

# Compute dynamics
dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                  start_time = 0,dt=0.1,num_steps=num_steps)

t2, s_z2 = dynamics.expectations(sigma_z, real=True)


# Plot
plt.plot(t1, s_z, label=r'TEMPO $\alpha={0}$ '.format(alpha))
plt.plot(t2, s_z2, label=r'QME $\gamma(t)$ ')
plt.title('Decay from Up Density Matrix State at T = {0} K'.format(T)) 
plt.xlabel(r'$t\,\Omega$')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)
