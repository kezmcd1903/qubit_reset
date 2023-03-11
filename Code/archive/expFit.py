"""
System plot giving TEMPO Hamiltonians and plot with derived master equation
 in linblad form
 
here trying to fit an exponent to find decay rate values
"""

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


# Pick parameters
Omega = 1.0
omega_cutoff = 5.0
alpha = 0.1
# Temperature in kelvin
T = 0 # K
# Temperature in units of 1/(ps kB)
# 1 K = 0.1309 1/(ps kB)
temperature = T*0.1309 # 1/(ps kB)

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
# Roos: dkmax=200,epsrel=2.5*10**-7 # dt = 0.06 dkmax=200,\
                                         # epsrel=10**(-7)
tempo_parameters = oqupy.TempoParameters(dt=0.06 , dkmax=200,\
                                          epsrel=10**(-7))

dynamics = oqupy.tempo_compute(system=system,
                                bath=bath,
                                initial_state=init_st,
                                start_time=0.0,
                                end_time=10.0,
                                parameters=tempo_parameters)
#and plot the result
t, s_z = dynamics.expectations(sigma_z, real=True)


##############################################################################
# Derived master equation

# B-E Occupancy
def N(T):
    if T == 0:
        N = 0
    else:
        N = (1/(math.exp(((1.055e-34)*(10**12))/((1.381e-23)*T))-1))
    return N

# Checks
print('T is', T)        
print('N is', N(T))

# Spectral density ohmic
J = 2*alpha*Omega*np.exp(-Omega/omega_cutoff) #replaced the 2 with 4
# # Spectral density super-ohmic
# J = 2*alpha*(Omega**(3)/omega_cutoff**(2))*np.exp(-Omega/omega_cutoff)
#gamma1 = gamma*(N+1)
gamma_1 = J*np.pi*(N(T)+1)#*(N+1)
# gamma2 = gamma*(N)
gamma_2 = J*np.pi*N(T) #*N

Gamma = [gamma_1, gamma_2]

# This is in the same form as in notes
lindblad_operators  = [2*sigma_minus, 2*sigma_plus] # Factors of 2 here!?
    

systemD = oqupy.System(
        Omega * sigma_z,
        gammas=Gamma,
        lindblad_operators=lindblad_operators)

# up_density_matrix
dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                  start_time = 0,dt=0.1,num_steps=100)


t2, s_z2 = dynamics.expectations(sigma_z, real=True)
#######################################################

# # Plot exponent
# t3 = np.linspace(0, 800, 800)
# a = 2
# # Change for rate
# b = 0.001
# c = 1
# y = a*np.exp(-b*t3)-c



plt.plot(t, s_z, label=r'TEMPO $\alpha={0}$ '.format(alpha))
plt.plot(t2, s_z2, label=r'QME $\gamma = {0} \pi$ '.format(J.round(3)))
# Plot exp
# plt.plot(t3, y, label='2exp(-{0}t)-1 '.format(b))
#
plt.title('Decay from Up Density Matrix State at T = {0} K'.format(T)) 
plt.xlabel(r'$t\,\Omega$')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()

stop = timeit.default_timer()
print('Time: ', stop - start)