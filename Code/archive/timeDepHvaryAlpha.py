"""
Time dependent system hamiltonian
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
omega_cutoff = 5
alpha = 1
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
dur = 15.0 # can go up to current PT length

# For Master Equation (need to change if changing dtME)
num_steps = dur*10

# Define time dependent eigenstate splitting
# # For centering gauss
# t_0 = 2.0
    
# def gaussian_shape(t, area = 1.0, tau = 1.0, t_0 = t_0): #area was 1
#     return area/(tau*np.sqrt(np.pi)) * np.exp(-(t-t_0)**2/(tau**2))
coeff = 1
def omega(t):
    return coeff*t

# # Shape of pulse
# detuning = lambda t: 0.0 * t
# Omega_t = gaussian_shape(t, area = np.pi/2.0, tau = 0.245)
# Delta_t = detuning(t)

# # Spectral density ohmic
def J(t):
    return 4*alpha*omega(t)*np.exp(-2*omega(t)/omega_cutoff)
# # Super-ohmic
# def J(t):
#     return 4*alpha*(omega(t)**3/omega_cutoff**(2))*\
#         np.exp(-2*omega(t)/omega_cutoff)

# Make plot of spectral density and time dependent splitting
t_plot = np.linspace(0,int(dur),int(dur))
plt.title('Omega(t) and Spectral Density')
plt.plot(t_plot, J(t_plot),label=r"$J(\Omega(t))$") 
plt.plot(t_plot, omega(t_plot),label=r"$\Omega(t)$={0}t".format(coeff))
plt.xlabel("Time (ps)")
plt.ylabel(r"$\omega(t) \mathrm{ps}^{-1}$")
plt.legend()
plt.show()

# Define time dependent system
# (Time dependent system gaussian_shape(t, area = np.pi/2.0, tau = 0.245))
def hamiltonian_t(t):
    return omega(t)*sigma_z

system = oqupy.TimeDependentSystem(hamiltonian_t)
correlations = oqupy.PowerLawSD(alpha=alpha,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='gaussian',
                                temperature=temperature)
bath = oqupy.Bath(sigma_x, correlations)

# Compute dynamics with TEMPO
tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
                                          epsrel=epsrel)

tempo_sys = oqupy.Tempo(system=system,
                        bath=bath,
                        initial_state=init_st,
                        start_time=0.0,
                        parameters=tempo_parameters)
dynamics = tempo_sys.compute(end_time=dur)
    
# Plot the result
t, s_z = dynamics.expectations(sigma_z, real=True)

# Save
desired_folder_data = 'C:/Users/sony/Desktop/Project/PlotData'

file_name_t = 'Om_'+str(coeff)+'t_alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'
    
file_name_s_z = 'Om_'+str(coeff)+'t_alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'
    
full_path_t = os.path.join(desired_folder_data, file_name_t)
full_path_s_z = os.path.join(desired_folder_data, file_name_s_z)

np.savetxt(full_path_t, t)
np.savetxt(full_path_s_z, s_z)

# Define time dependent system
# (Time dependent system gaussian_shape(t, area = np.pi/2.0, tau = 0.245))
alpha2 = 0.1

def hamiltonian_t(t):
    return omega(t)*sigma_z

system = oqupy.TimeDependentSystem(hamiltonian_t)
correlations = oqupy.PowerLawSD(alpha=alpha2,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='gaussian',
                                temperature=temperature)
bath = oqupy.Bath(sigma_x, correlations)

# Compute dynamics with TEMPO
tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
                                          epsrel=epsrel)

tempo_sys = oqupy.Tempo(system=system,
                        bath=bath,
                        initial_state=init_st,
                        start_time=0.0,
                        parameters=tempo_parameters)
dynamics = tempo_sys.compute(end_time=dur)
    
# Plot the result
t2, s_z2 = dynamics.expectations(sigma_z, real=True)

# Define time dependent system
# (Time dependent system gaussian_shape(t, area = np.pi/2.0, tau = 0.245))
alpha3 = 0.01

def hamiltonian_t(t):
    return omega(t)*sigma_z

system = oqupy.TimeDependentSystem(hamiltonian_t)
correlations = oqupy.PowerLawSD(alpha=alpha2,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='gaussian',
                                temperature=temperature)
bath = oqupy.Bath(sigma_x, correlations)

# Compute dynamics with TEMPO
tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
                                          epsrel=epsrel)

tempo_sys = oqupy.Tempo(system=system,
                        bath=bath,
                        initial_state=init_st,
                        start_time=0.0,
                        parameters=tempo_parameters)
dynamics = tempo_sys.compute(end_time=dur)
    
# Plot the result
t3, s_z3 = dynamics.expectations(sigma_z, real=True)

# Derived master equation

# # B-E Occupancy (for temp dependence)
# def N(T):
#     if T == 0:
#         N = 0
#     else:
#         N = (1/(math.exp(((1.055e-34)*(10**12))/((1.381e-23)*T))-1))
#     return N

# # Decay factors
# #gamma1 = gamma*(N+1)
# gamma_1 = 2*np.pi*(N(T)+1)#*(N+1)
# # gamma2 = gamma*(N)
# gamma_2 = 2*np.pi*N(T) #*N

# Gamma = [ lambda t: gamma_1*J(t), lambda t: gamma_2*J(t)]

# lindblad_operators = [ lambda t: sigma_minus, lambda t: sigma_plus]

    
# systemD = oqupy.TimeDependentSystem(hamiltonian_t,
#                                    gammas=Gamma,
#                                    lindblad_operators=lindblad_operators)

# # systemD = oqupy.System(
# #          gaussian_shape(t, area = np.pi/2.0, tau = 0.245) * sigma_z,
# #         gammas=Gamma,
# #         lindblad_operators=lindblad_operators)

# # up_density_matrix

# # Compute dynamics
# dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
#                                   start_time = 0,dt=0.1,num_steps=num_steps)

# t2, s_z2 = dynamics.expectations(sigma_z, real=True)


# Plot
plt.plot(t, s_z, label=r'TEMPO $\alpha={0}$ '.format(alpha))
plt.plot(t2, s_z2, label=r'TEMPO $\alpha={0}$ '.format(alpha2))
plt.plot(t3, s_z3, label=r'TEMPO $\alpha={0}$ '.format(alpha3))
#plt.plot(t2, s_z2, label=r'QME $\gamma(t)$ ')
plt.title('Decay from Up Density Matrix State at T = {0} K'.format(T)) 
plt.xlabel(r'$t\,\Omega$')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
# Save plots
desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'Om_'+str(coeff)+'t_alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)
