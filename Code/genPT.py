"""

Generate Process Tensors

- For differing alpha and duration
- Can be used for different system hamiltonians - time dependent of independent

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
# can go up to current PT length
#dur = [1000,1000,1000,1000]
dur = 7
alpha = 0.1
# Pick Ramping Coefficients
coeff = 1


# Change coefficient?
def omega(t):
    return coeff*t


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

#Generate PT
# Compute dynamics with TEMPO
tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
                                          epsrel=epsrel)
    
process_tensor = oqupy.pt_tempo_compute(bath=bath,
                                        start_time=0,
                                        end_time=dur,
                                        parameters=tempo_parameters)
# ###### Save PT ######## Doesn't work if already saved
desired_folder_PT = 'C:/Users/sony/Desktop/Project/process_tensor/'

file_name_PT = 'Om_'+str(coeff)+'t_alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_PT.txt'

full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

oqupy.SimpleProcessTensor.export(self=process_tensor,
                                  filename=full_path_PT)

###################################################################
# ###################################################################



# # Change coefficient?
# def omega(t):
#     return coeff[1]*t


# # Define time dependent system
# # (Time dependent system gaussian_shape(t, area = np.pi/2.0, tau = 0.245))
# def hamiltonian_t(t):
#     return omega(t)*sigma_z

# system = oqupy.TimeDependentSystem(hamiltonian_t)
# correlations = oqupy.PowerLawSD(alpha=alpha[1],
#                                 zeta=1,
#                                 cutoff=omega_cutoff,
#                                 cutoff_type='gaussian',
#                                 temperature=temperature)
# bath = oqupy.Bath(sigma_x, correlations)

# #Generate PT
# # Compute dynamics with TEMPO
# tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
#                                           epsrel=epsrel)
    
# process_tensor = oqupy.pt_tempo_compute(bath=bath,
#                                         start_time=0,
#                                         end_time=dur[1],
#                                         parameters=tempo_parameters)
# # ###### Save PT ######## Doesn't work if already saved
# desired_folder_PT = 'C:/Users/sony/Desktop/Project/process_tensor/'

# file_name_PT = 'Om_'+str(coeff[1])+'t_alpha='+str(alpha[1])+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur[1])+'_PT.txt'

# full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

# oqupy.SimpleProcessTensor.export(self=process_tensor,
#                                   filename=full_path_PT)

# #########################################################################
# #########################################################################


# # Change coefficient?
# def omega(t):
#     return coeff[2]*t


# # Define time dependent system
# # (Time dependent system gaussian_shape(t, area = np.pi/2.0, tau = 0.245))
# def hamiltonian_t(t):
#     return omega(t)*sigma_z

# system = oqupy.TimeDependentSystem(hamiltonian_t)
# correlations = oqupy.PowerLawSD(alpha=alpha[2],
#                                 zeta=1,
#                                 cutoff=omega_cutoff,
#                                 cutoff_type='gaussian',
#                                 temperature=temperature)
# bath = oqupy.Bath(sigma_x, correlations)

# #Generate PT
# # Compute dynamics with TEMPO
# tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
#                                           epsrel=epsrel)
    
# process_tensor = oqupy.pt_tempo_compute(bath=bath,
#                                         start_time=0,
#                                         end_time=dur[2],
#                                         parameters=tempo_parameters)
# # ###### Save PT ######## Doesn't work if already saved
# desired_folder_PT = 'C:/Users/sony/Desktop/Project/process_tensor/'

# file_name_PT = 'Om_'+str(coeff[2])+'t_alpha='+str(alpha[2])+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur[2])+'_PT.txt'

# full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

# oqupy.SimpleProcessTensor.export(self=process_tensor,
#                                   filename=full_path_PT)

# ############################################################################
# ############################################################################

# # Change coefficient?
# def omega(t):
#     return coeff[3]*t

# # Define time dependent system
# # (Time dependent system gaussian_shape(t, area = np.pi/2.0, tau = 0.245))
# def hamiltonian_t(t):
#     return omega(t)*sigma_z

# system = oqupy.TimeDependentSystem(hamiltonian_t)
# correlations = oqupy.PowerLawSD(alpha=alpha[3],
#                                 zeta=1,
#                                 cutoff=omega_cutoff,
#                                 cutoff_type='gaussian',
#                                 temperature=temperature)
# bath = oqupy.Bath(sigma_x, correlations)

# #Generate PT
# # Compute dynamics with TEMPO
# tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
#                                           epsrel=epsrel)
    
# process_tensor = oqupy.pt_tempo_compute(bath=bath,
#                                         start_time=0,
#                                         end_time=dur[3],
#                                         parameters=tempo_parameters)
# # ###### Save PT ######## Doesn't work if already saved
# desired_folder_PT = 'C:/Users/sony/Desktop/Project/process_tensor/'

# file_name_PT = 'Om_'+str(coeff[3])+'t_alpha='+str(alpha[3])+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur[3])+'_PT.txt'

# full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

# oqupy.SimpleProcessTensor.export(self=process_tensor,
#                                   filename=full_path_PT)