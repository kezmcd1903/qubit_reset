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
omega_cutoff = 5
alpha3 = 0.001
# Temperature in kelvin
T = 0 # K
# Temperature in units of 1/(ps kB). {1 K = 0.1309 1/(ps kB)}
temperature = T*0.1309   # 1/(ps kB)
# initial state = up_density_matrix or mixed_density_matrix
init_st = up_density_matrix

# TEMPO parameters
dt=0.05
dkmax=200
epsrel=10**(-6)
# duration (ps)
#dur = 3000 # can go up to current PT length
dur = [2,20,200]
# Test values
#dur = [1,1,1]
# # For Master Equation (need to change if changing dtME)
# num_steps = dur
coeff1 = 1
coeff2 = 0.1
coeff3 = 0.01

# Change coefficient?

def omega(t):
    return coeff1*t

# # Spectral density ohmic
def J(t):
    return 4*alpha3*omega(t)*np.exp(-2*omega(t)/omega_cutoff)
# # Super-ohmic
# def J(t):
#     return 4*alpha*(omega(t)**3/omega_cutoff**(2))*\
#         np.exp(-2*omega(t)/omega_cutoff)

# # Make plot of spectral density and time dependent splitting
# t = np.linspace(0,int(dur),int(dur)*100) #
# plt.title('Energy Splitting and Spectral Density')
# plt.plot(t, J(t),label=r"Sup $J(\Omega(t))$") 
# plt.plot(t, omega(t),label=r"$\Omega(t)$={0}t".format(coeff1))
# plt.xlabel("Time (ps)")
# plt.ylabel(r"$\omega(t) \mathrm{ps}^{-1}$")
# plt.legend()
# # Save plots
# desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
# file_name = 'Om_'+str(coeff)+'t_alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Spect&Om.pdf'
# full_path = os.path.join(desired_folder, file_name)
# plt.savefig(full_path)
plt.show()

# Define time dependent system
# (Time dependent system gaussian_shape(t, area = np.pi/2.0, tau = 0.245))
def hamiltonian_t(t):
    return omega(t)*sigma_z

system = oqupy.TimeDependentSystem(hamiltonian_t)
correlations = oqupy.PowerLawSD(alpha=alpha3,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='gaussian',
                                temperature=temperature)
bath = oqupy.Bath(sigma_x, correlations)

# #Generate PT
# # Compute dynamics with TEMPO
# tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
#                                           epsrel=epsrel)
    
# process_tensor = oqupy.pt_tempo_compute(bath=bath,
#                                         start_time=0,
#                                         end_time=dur[0],
#                                         parameters=tempo_parameters)
# # ###### Save PT ######## Doesn't work if already saved
# desired_folder_PT = 'C:/Users/sony/Desktop/Project/process_tensor/'

# file_name_PT = 'Om_'+str(coeff1)+'t_alpha='+str(alpha3)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_PT.txt'

# full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

# oqupy.SimpleProcessTensor.export(self=process_tensor,
#                                   filename=full_path_PT)

# Read in PT
full_path_PT = 'C:/Users/sony/Desktop/Project/process_tensor/Om_0.05t_alpha=0.001_T=0_dt=0.05_dkmax=200_epsrel=1e-06_dur=200_PT.txt'

process_tensor = oqupy.process_tensor.FileProcessTensor(filename=full_path_PT,
                                                        mode = 'read')



dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
                                system=system,
                                initial_state=init_st,
                                start_time=0.0 # needed?, progress_type
                                ,num_steps=dur[0]/dt) #NUm steps correct? YES (PT has certain size when made)
    
# Plot the result
t1, s_z1 = dynamics.expectations(sigma_z, real=True)

# # Save
# desired_folder_data = 'C:/Users/sony/Desktop/Project/PlotData'

# file_name_t = 'Om_'+str(coeff)+'t_alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'
    
# file_name_s_z = 'Om_'+str(coeff)+'t_alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'
    
# full_path_t = os.path.join(desired_folder_data, file_name_t)
# full_path_s_z = os.path.join(desired_folder_data, file_name_s_z)

# np.savetxt(full_path_t, t)
# np.savetxt(full_path_s_z, s_z)

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

# alpha
# up_density_matrix
num_steps = dur[0]
# Compute dynamics
dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                  start_time = 0,dt=1,num_steps=num_steps)

tME1, s_zME1 = dynamics.expectations(sigma_z, real=True)
##################################################################
#################################################################
# # #######################################################################
# # ######################################################################


# Change coefficient?

def omega(t):
    return coeff2*t

# # Spectral density ohmic
def J(t):
    return 4*alpha3*omega(t)*np.exp(-2*omega(t)/omega_cutoff)


# Define time dependent system
# (Time dependent system gaussian_shape(t, area = np.pi/2.0, tau = 0.245))
def hamiltonian_t(t):
    return omega(t)*sigma_z

system = oqupy.TimeDependentSystem(hamiltonian_t)
correlations = oqupy.PowerLawSD(alpha=alpha3,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='gaussian',
                                temperature=temperature)
bath = oqupy.Bath(sigma_x, correlations)

# #Generate PT
# # Compute dynamics with TEMPO
# tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
#                                           epsrel=epsrel)
    
# process_tensor = oqupy.pt_tempo_compute(bath=bath,
#                                         start_time=0,
#                                         end_time=dur,
#                                         parameters=tempo_parameters)
# # ###### Save PT ######## Doesn't work if already saved
# desired_folder_PT = 'C:/Users/sony/Desktop/Project/process_tensor/'

# file_name_PT = 'Om_'+str(coeff2)+'t_alpha='+str(alpha3)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_PT.txt'

# full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

# oqupy.SimpleProcessTensor.export(self=process_tensor,
#                                   filename=full_path_PT)


# Read in PT
full_path_PT = 'C:/Users/sony/Desktop/Project/process_tensor/Om_0.05t_alpha=0.001_T=0_dt=0.05_dkmax=200_epsrel=1e-06_dur=200_PT.txt'

process_tensor = oqupy.process_tensor.FileProcessTensor(filename=full_path_PT,
                                                        mode = 'read')

dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
                                system=system,
                                initial_state=init_st,
                                start_time=0.0 # needed?, progress_type
                                ,num_steps=dur[1]/dt) #NUm steps correct? YES (PT has certain size when made)
    
# Plot the result
t2, s_z2 = dynamics.expectations(sigma_z, real=True)


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

num_steps = dur[1]
# Compute dynamics
dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                  start_time = 0,dt=1,num_steps=num_steps)

tME2, s_zME2 = dynamics.expectations(sigma_z, real=True)
##################################################################
#################################################################
# ##############################################
# # #################################################################





# Change coefficient?
def omega(t):
    return coeff3*t

# # Spectral density ohmic
def J(t):
    return 4*alpha3*omega(t)*np.exp(-2*omega(t)/omega_cutoff)


# Define time dependent system
# (Time dependent system gaussian_shape(t, area = np.pi/2.0, tau = 0.245))
def hamiltonian_t(t):
    return omega(t)*sigma_z

system = oqupy.TimeDependentSystem(hamiltonian_t)
correlations = oqupy.PowerLawSD(alpha=alpha3,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='gaussian',
                                temperature=temperature)
bath = oqupy.Bath(sigma_x, correlations)

# #Generate PT
# # Compute dynamics with TEMPO
# tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
#                                           epsrel=epsrel)
    
# process_tensor = oqupy.pt_tempo_compute(bath=bath,
#                                         start_time=0,
#                                         end_time=dur,
#                                         parameters=tempo_parameters)
# # ###### Save PT ######## Doesn't work if already saved
# desired_folder_PT = 'C:/Users/sony/Desktop/Project/process_tensor/'

# file_name_PT = 'Om_'+str(coeff3)+'t_alpha='+str(alpha3)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_PT.txt'

# full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

# oqupy.SimpleProcessTensor.export(self=process_tensor,
#                                   filename=full_path_PT)

# Read in PT
full_path_PT = 'C:/Users/sony/Desktop/Project/process_tensor/Om_0.05t_alpha=0.001_T=0_dt=0.05_dkmax=200_epsrel=1e-06_dur=200_PT.txt'

process_tensor = oqupy.process_tensor.FileProcessTensor(filename=full_path_PT,
                                                        mode = 'read')

dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
                                system=system,
                                initial_state=init_st,
                                start_time=0.0 # needed?, progress_type
                                ,num_steps=dur[2]/dt) #NUm steps correct? YES (PT has certain size when made)
    
# Plot the result
t3, s_z3 = dynamics.expectations(sigma_z, real=True)

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

# alpha
# up_density_matrix
num_steps = dur[2]
# Compute dynamics
dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                  start_time = 0,dt=1,num_steps=num_steps)

tME3, s_zME3 = dynamics.expectations(sigma_z, real=True)
##################################################################
#################################################################
# ########################################

# # Change coefficient?
# coeff4 = 0.001
# def omega(t):
#     return coeff4*t

# # # Spectral density ohmic
# def J(t):
#     return 4*alpha3*omega(t)*np.exp(-2*omega(t)/omega_cutoff)


# # Define time dependent system
# # (Time dependent system gaussian_shape(t, area = np.pi/2.0, tau = 0.245))
# def hamiltonian_t(t):
#     return omega(t)*sigma_z

# system = oqupy.TimeDependentSystem(hamiltonian_t)
# correlations = oqupy.PowerLawSD(alpha=alpha3,
#                                 zeta=1,
#                                 cutoff=omega_cutoff,
#                                 cutoff_type='gaussian',
#                                 temperature=temperature)
# bath = oqupy.Bath(sigma_x, correlations)

# # #Generate PT
# # # Compute dynamics with TEMPO
# # tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
# #                                           epsrel=epsrel)
    
# # process_tensor = oqupy.pt_tempo_compute(bath=bath,
# #                                         start_time=0,
# #                                         end_time=dur,
# #                                         parameters=tempo_parameters)
# # # ###### Save PT ######## Doesn't work if already saved
# # desired_folder_PT = 'C:/Users/sony/Desktop/Project/process_tensor/'

# # file_name_PT = 'Om_'+str(coeff4)+'t_alpha='+str(alpha3)+'_T='+str(T)+'_dt='+str(dt)+\
# #     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_PT.txt'

# # full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

# # oqupy.SimpleProcessTensor.export(self=process_tensor,
# #                                   filename=full_path_PT)

# # Read in PT
# full_path_PT = 'C:/Users/sony/Desktop/Project/process_tensor/Om_0.001t_alpha=0.001_T=0_dt=0.05_dkmax=200_epsrel=1e-06_dur=200_PT.txt'

# process_tensor = oqupy.process_tensor.FileProcessTensor(filename=full_path_PT,
#                                                         mode = 'read')                                                     

# dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
#                                 system=system,
#                                 initial_state=init_st,
#                                 start_time=0.0 # needed?, progress_type
#                                 ,num_steps=dur/dt) #NUm steps correct? YES (PT has certain size when made)
    
# # Plot the result
# t4, s_z4 = dynamics.expectations(sigma_z, real=True)

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
#                                     gammas=Gamma,
#                                     lindblad_operators=lindblad_operators)

# # alpha
# # up_density_matrix

# # Compute dynamics
# dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
#                                   start_time = 0,dt=1,num_steps=num_steps)

# tME4, s_zME4 = dynamics.expectations(sigma_z, real=True)


################################################################
###################################################

# Plot
plt.plot(t1*coeff1, s_z1, label=r'TEMPO $\Omega(t)$={1}t'.format(alpha3,coeff1))
plt.plot(tME1*coeff1, s_zME1, label=r'QME $\gamma(t),\Omega(t)$={0}t'.format(coeff1))
plt.plot(t2*coeff2, s_z2, label=r'TEMPO $\Omega(t)$={1}t'.format(alpha3,coeff2))
plt.plot(tME2*coeff2, s_zME2, label=r'QME $\gamma(t),\Omega(t)$={0}t'.format(coeff2))
plt.plot(t3*coeff3, s_z3, label=r'TEMPO $\Omega(t)$={1}t'.format(alpha3,coeff3))
plt.plot(tME3*coeff3, s_zME3, label=r'QME $\gamma(t),\Omega(t)$={0}t'.format(coeff3))
# plt.plot(t4, s_z4, label=r'TEMPO $\Omega(t)$={1}t'.format(alpha3,coeff4))
# plt.plot(tME4, s_zME4, label=r'QME $\gamma(t),\Omega(t)$={0}t'.format(coeff4))

plt.title(r'Decay at T = {0} K for $\alpha = {1} $'.format(T,alpha3)) 
plt.xlabel(r'$\Omega$')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend(loc='lower left') #loc='upper right'
# Save plots
desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'Scal_varyingCoeff_a=0.001'+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)
