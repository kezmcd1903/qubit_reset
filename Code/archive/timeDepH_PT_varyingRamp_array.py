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
dur = 2 # can go up to current PT length

# For Master Equation (need to change if changing dtME)
num_steps = dur/dt

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

# file_name_PT = 'Om_'+str(coeff1)+'t_alpha='+str(alpha3)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_PT.txt'

# full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

# oqupy.SimpleProcessTensor.export(self=process_tensor,
#                                   filename=full_path_PT)

# Read in PT
full_path_PT = 'C:/Users/sony/Desktop/Project/process_tensor/Om_1t_alpha=0.001_T=0_dt=0.05_dkmax=200_epsrel=1e-06_dur=200_PT.txt'

process_tensor = oqupy.process_tensor.FileProcessTensor(filename=full_path_PT,
                                                        mode = 'read')

coeffs = [1,0.1,0.01,0.001,0.0001]
systems = []

for coeff in coeffs:
    # NOTE: omitting "delta=delta" in the parameter definition below
    #       would lead to all systems having the same detuning.
    #       This is a common python pitfall. Check out
    #       https://docs.python-guide.org/writing/gotchas/#late-binding-closures
    #       for more information on this.
    def hamiltonian_t(t, coeff=coeff):
        return coeff*t*sigma_z
    
    system = oqupy.TimeDependentSystem(hamiltonian_t)
    systems.append(system)

s_z_list = []
t_list = []
for system in systems:
    dynamics = oqupy.compute_dynamics(
        process_tensor=process_tensor,
        system=system,
        initial_state=init_st,
        start_time=0,
        num_steps=num_steps,
        progress_type="silent")
    t, s_z = dynamics.expectations(sigma_z, real=True)
    s_z_list.append(s_z)
    t_list.append(t)
    print(".", end="", flush=True)
print(" done.", flush=True)


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


##################################################################
#################################################################


# Plot
for t, s_z, coeff in zip(t_list, s_z_list, coeffs):
    plt.plot(t, s_z, label="coeff ="+f"{coeff:0.1f}")
    plt.xlabel(r'$t/$ps')
    plt.ylabel(r'$<\sigma_{z}>$')
plt.ylim((-1.0,1.0))
plt.legend()
plt.show()


# ##############
# plt.plot(t1, s_z1, label=r'TEMPO $\Omega(t)$={1}t'.format(alpha3,coeff1))


# plt.title(r'Decay at T = {0} K for $\alpha = {1} $'.format(T,alpha3)) 
# plt.xlabel(r'$t\,\Omega$')
# plt.ylabel(r'$\langle\sigma_z\rangle$')
# plt.legend()
# # Save plots
# desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
# file_name = 'varyingCoeff_a=0.001'+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
# full_path = os.path.join(desired_folder, file_name)
# plt.savefig(full_path)
# plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)
