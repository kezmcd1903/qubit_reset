"""

Solves Time Dependent System Hamiltonians

- Input process tensor
- Solve for 4 different system hamiltonian with different time
dependent energy splittings
- Solves using TEMPO and master equation
- Plots both.

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



# Pick system parameters
Omega = 1
omega_cutoff = 5
alpha = 10
# Temperature in kelvin
T = 0 # K
# Temperature in units of 1/(ps kB). {1 K = 0.1309 1/(ps kB)}
temperature = T*0.1309   # 1/(ps kB)
# initial state = up_density_matrix or mixed_density_matrix
init_st = up_density_matrix

# TEMPO parameters
dt=0.01
dkmax=200
epsrel=10**(-6)
# duration (ps)
dur = 0.5 # can go up to current PT length

# Time step for Master Equation
dtME = 0.01
num_stepsME = dur/dtME


# Change coefficient for different time dependent systems
coeff = [1,1.5,2,2.5]
b = 1

# Different time dependent energy splittings
def omega1(t):
    return coeff[0]*t + b

def omega2(t):
    return coeff[1]*t + b

def omega3(t):
    return coeff[2]*t + b

def omega4(t):
    return coeff[3]*t + b



# Make plot of spectral density
w = np.linspace(0,25,250)
plt.title('Spectral Density')
J = 4*alpha*w*np.exp(-2*w/omega_cutoff)
plt.plot(w, J)
plt.axvline(x = 2.5,color='r',linestyle='dashed',
            label=r"$\Omega_{0}={1}$".format(0,2.5))
plt.xlabel(r'$\Omega (ps^{-1})$')
plt.ylabel(r"$J(\Omega) \mathrm{ps}^{-1}$")
# Save plots
desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'J(Om)_'+str(coeff)+'_b='+str(b)+'t_alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Spect&Om.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()

# Make plot of time dependent splitting
t = np.linspace(0,4,400) #int(dur),int(dur)*100)
plt.title('Varying Ramp - Time Dependent Energy Splitting')
plt.axhline(y = 2.5,color='c',linestyle='dashed',
            label=r"$\Omega_{0}={1}$".format('p',2.5))
plt.plot(t, omega1(t),color='r',label=r"$\Omega(t)$={0}t+{1}".format(coeff[0],b))
plt.plot(t, omega2(t),color='b',label=r"$\Omega(t)$={0}t+{1}".format(coeff[1],b))
plt.plot(t, omega3(t),color='g',label=r"$\Omega(t)$={0}t+{1}".format(coeff[2],b))
plt.plot(t, omega4(t),color='m',label=r"$\Omega(t)$={0}t+{1}".format(coeff[3],b))
plt.xlabel("Time (ps)")
plt.ylabel(r'$\Omega \; (ps^{-1})$')
plt.legend()
# Save plots
desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'Om_'+str(coeff)+'_b='+str(b)+'t_alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Spect&Om.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()

# # Spectral density ohmic
def J(t):
    return 4*alpha*omega1(t)*np.exp(-2*omega1(t)/omega_cutoff)

# B-E Occupancy (for temp dependence in master equation)
def N(T):
    if T == 0:
        N = 0
    else:
        N = (1/(math.exp(((1.055e-34)*(10**12))/((1.381e-23)*T))-1))
    return N

# Time dependent system
def hamiltonian_t(t):
    return omega1(t)*sigma_z

system = oqupy.TimeDependentSystem(hamiltonian_t)
correlations = oqupy.PowerLawSD(alpha=alpha,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='gaussian',
                                temperature=temperature)
bath = oqupy.Bath(sigma_x, correlations)

#Generate process tensor (PT) if one is not already made 

# # Compute dynamics with TEMPO
# tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
#                                           epsrel=epsrel)
    
# process_tensor = oqupy.pt_tempo_compute(bath=bath,
#                                         start_time=0,
#                                         end_time=dur,
#                                         parameters=tempo_parameters)
# # ###### Save PT ######## Doesn't work if already saved
# desired_folder_PT = 'C:/Users/sony/Desktop/Project/process_tensor/'

# file_name_PT = 'Om_'+str(coeff)+'t_alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_PT.txt'

# full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

# oqupy.SimpleProcessTensor.export(self=process_tensor,
#                                   filename=full_path_PT)

# Read in PT
full_path_PT = 'C:/Users/sony/Desktop/Project/process_tensor/Om_[5, 20, 80, 100]t_alpha=10_T=0_dt=0.01_dkmax=200_epsrel=1e-06_dur=0.5_PT.txt'

process_tensor = oqupy.process_tensor.FileProcessTensor(filename=full_path_PT,
                                                        mode = 'read')

dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
                                system=system,
                                initial_state=init_st,
                                start_time=0.0 # needed?, progress_type
                                ,num_steps=dur/dt) #NUm steps correct? YES (PT has certain size when made)
    
# Plot the result
t1, s_z1 = dynamics.expectations(sigma_z, real=True)

# # Can save plot data
# desired_folder_data = 'C:/Users/sony/Desktop/Project/PlotData'

# file_name_t = 'Om_'+str(coeff)+'t_alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'
    
# file_name_s_z = 'Om_'+str(coeff)+'t_alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'
    
# full_path_t = os.path.join(desired_folder_data, file_name_t)
# full_path_s_z = os.path.join(desired_folder_data, file_name_s_z)

# np.savetxt(full_path_t, t)
# np.savetxt(full_path_s_z, s_z)


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


# Compute dynamics
dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                  start_time = 0,dt=dtME,num_steps=num_stepsME)

tME1, s_zME1 = dynamics.expectations(sigma_z, real=True)
##################################################################
# #################################################################
# # # #######################################################################
# # # ######################################################################


# Next time dependent system

# # Spectral density ohmic
def J(t):
    return 4*alpha*omega2(t)*np.exp(-2*omega2(t)/omega_cutoff)


# Define time dependent system
# (Time dependent system gaussian_shape(t, area = np.pi/2.0, tau = 0.245))
def hamiltonian_t(t):
    return omega2(t)*sigma_z

system = oqupy.TimeDependentSystem(hamiltonian_t)
correlations = oqupy.PowerLawSD(alpha=alpha,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='gaussian',
                                temperature=temperature)
bath = oqupy.Bath(sigma_x, correlations)


dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
                                system=system,
                                initial_state=init_st,
                                start_time=0.0 # needed?, progress_type
                                ,num_steps=dur/dt) #NUm steps correct? YES (PT has certain size when made)
    
# Plot the result
t2, s_z2 = dynamics.expectations(sigma_z, real=True)


# Master Equation

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


# Compute dynamics
dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                  start_time = 0,dt=dtME,num_steps=num_stepsME)

tME2, s_zME2 = dynamics.expectations(sigma_z, real=True)
##################################################################
#################################################################
# ##############################################
# # #################################################################

# Repeat for different system

# # Spectral density ohmic
def J(t):
    return 4*alpha*omega3(t)*np.exp(-2*omega3(t)/omega_cutoff)


# Define time dependent system
# (Time dependent system gaussian_shape(t, area = np.pi/2.0, tau = 0.245))
def hamiltonian_t(t):
    return omega3(t)*sigma_z

system = oqupy.TimeDependentSystem(hamiltonian_t)
correlations = oqupy.PowerLawSD(alpha=alpha,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='gaussian',
                                temperature=temperature)
bath = oqupy.Bath(sigma_x, correlations)

dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
                                system=system,
                                initial_state=init_st,
                                start_time=0.0 # needed?, progress_type
                                ,num_steps=dur/dt)


# Plot the result
t3, s_z3 = dynamics.expectations(sigma_z, real=True)

# Master equation

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


# Compute dynamics
dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                  start_time = 0,dt=dtME,num_steps=num_stepsME)

tME3, s_zME3 = dynamics.expectations(sigma_z, real=True)
##################################################################
################################################################
########################################

# Repeat

# # Spectral density ohmic
def J(t):
    return 4*alpha*omega4(t)*np.exp(-2*omega4(t)/omega_cutoff)


# Define time dependent system
# (Time dependent system gaussian_shape(t, area = np.pi/2.0, tau = 0.245))
def hamiltonian_t(t):
    return omega4(t)*sigma_z

system = oqupy.TimeDependentSystem(hamiltonian_t)
correlations = oqupy.PowerLawSD(alpha=alpha,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='gaussian',
                                temperature=temperature)
bath = oqupy.Bath(sigma_x, correlations)
                                                   

dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
                                system=system,
                                initial_state=init_st,
                                start_time=0.0 # needed?, progress_type
                                ,num_steps=dur/dt) #NUm steps correct? YES (PT has certain size when made)
    
# Plot the result
t4, s_z4 = dynamics.expectations(sigma_z, real=True)

# Master equation

# Decay factors
#gamma1 = gamma*(N+1)
gamma_1 = 2*np.pi*(N(T)+1)
# gamma2 = gamma*(N)
gamma_2 = 2*np.pi*N(T) 

Gamma = [ lambda t: gamma_1*J(t), lambda t: gamma_2*J(t)]

lindblad_operators = [ lambda t: sigma_minus, lambda t: sigma_plus]


systemD = oqupy.TimeDependentSystem(hamiltonian_t,
                                    gammas=Gamma,
                                    lindblad_operators=lindblad_operators)


# Compute dynamics
dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                  start_time = 0,dt=dtME,num_steps=num_stepsME)

tME4, s_zME4 = dynamics.expectations(sigma_z, real=True)


#################################################################
###################################################

#  Plot TEMPO Solution
plt.plot(t1, s_z1,color='r',label=r'TEMPO $\Omega_{0}={1}t+{2}$ '.format(0,coeff[0],b))
plt.plot(t2, s_z2,color='b',label=r'TEMPO $\Omega_{0}={1}t+{2}$ '.format(1,coeff[1],b))
plt.plot(t3, s_z3,color='g',label=r'TEMPO $\Omega_{0}={1}t+{2}$ '.format(2,coeff[2],b))
plt.plot(t4, s_z4,color='m',label=r'TEMPO $\Omega_{0}={1}t+{2}$ '.format(3,coeff[3],b))

plt.title(r'Decay from Up Density Matrix State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
plt.grid()
# Save plot
desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'varyingCoeff_a='+str(alpha)+'c='+str(coeff)+'_b='+str(b)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()

#  Plot TEMPO Solution Zoomed in
plt.plot(t1, s_z1,color='r',label=r'TEMPO $\Omega_{0}={1}t+{2}$ '.format(0,coeff[0],b))
plt.plot(t2, s_z2,color='b',label=r'TEMPO $\Omega_{0}={1}t+{2}$ '.format(1,coeff[1],b))
plt.plot(t3, s_z3,color='g',label=r'TEMPO $\Omega_{0}={1}t+{2}$ '.format(2,coeff[2],b))
plt.plot(t4, s_z4,color='m',label=r'TEMPO $\Omega_{0}={1}t+{2}$ '.format(3,coeff[3],b))

plt.title(r'Decay from Up Density Matrix State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
plt.grid()
plt.ylim([-0.2,0.1])
plt.xlim([dur/2, dur])

desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'varyingCoeffZm_a='+str(alpha)+'c='+str(coeff)+'_b='+str(b)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()


# Plot QME solution
plt.plot(tME1, s_zME1,color='r', label=r'QME $\Omega_{0}={1}t+{2}$ '.format(0,coeff[0],b))
plt.plot(tME2, s_zME2,color='b', label=r'QME $\Omega_{0}={1}t+{2}$ '.format(1,coeff[1],b))
plt.plot(tME3, s_zME3,color='g', label=r'QME $\Omega_{0}={1}t+{2}$ '.format(2,coeff[2],b))
plt.plot(tME4, s_zME4,color='m', label=r'QME $\Omega_{0}={1}t+{2}$ '.format(3,coeff[3],b))

plt.title(r'Decay from Up Density Matrix State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
plt.grid()

desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'varyingCoeffQME_a='+str(alpha)+'c='+str(coeff)+'_b='+str(b)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)
