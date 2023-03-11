"""
Time dependent system hamiltonian with process tensor

Looking at fastest decay time with PT-TEMPO (only), scaling x-axis as omega.
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
alpha = 0.1
# Temperature in kelvin
T = 0 # K
# Temperature in units of 1/(ps kB). {1 K = 0.1309 1/(ps kB)}
temperature = T*0.1309   # 1/(ps kB)
# initial state = up_density_matrix or mixed_density_matrix
init_st = up_density_matrix

dur = 4

# TEMPO parameters
dt=0.05
dkmax=70
epsrel=10**(-6)
# duration (ps) [can go up to current PT length]
# For Master Equation (need to change if changing dtME)
dtME = 0.1
num_stepsME = dur/dtME


#Define time dependent eigenstate splitting
# For centering gauss
t_0 = 2
area = 6.0
tau = 1.3
b = 1
 
def gaussian_shape(t, area = area, tau = tau, t_0 = t_0): #area was 1
    return area/(tau*np.sqrt(np.pi)) * np.exp(-(t-t_0)**2/(tau**2)) + b

# # For tanh shape time dependence
# D = 10
# a = 2
# b = 6
# c = 14

# def tanh(t,a=a,b=b,c=c,D=D):
#     return D*np.tanh(a*t-b)+c

# # Shape of pulse
# detuning = lambda t: 0.0 * t


# Make plot of spectral density and time dependent splitting
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
file_name = 'J(Om)_alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Spect&Om.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()

#t_plot0 = np.linspace(0,int(dur[0]),int(dur[0])*10) #
t_plot1 = np.linspace(0,int(dur),int(dur)*10)
plt.title('Gaussian Profile - Time Dependent Energy Splitting')
#plt.plot(t_plot0, omega(t_plot0),label=r"$\Omega(t)$={0}t".format(coeff[0])) 
plt.plot(t_plot1, gaussian_shape(t_plot1),label=r"$\Omega(t)= \frac{A}{\tau \sqrt{\pi}} \: e^{-\frac{(t-t_{o})^{2}}{\tau^{2}}+b} $")#.format(coeff) tau = 0.245
# #plt.axhline(y = Omega, xmin = 0, xmax = int(dur[2]), color = 'g',\
#             label=r"$\Omega={0}$".format(Omega))
plt.axhline(y = 2.5,color='r',linestyle='dashed',
            label=r"$\Omega_{0}={1}$".format('p',2.5))

plt.xlabel('Time (ps)')
plt.ylabel(r"$\Omega \; (\mathrm{ps}^{-1})$")
plt.legend()
# Save plots
desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
# pulse
# file_name = 'pulse_a='+str(alpha)+'_t_0='+str(t_0)+'_A='+str(area)+'_tau='+str(tau)+'_b='+str(b)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Spect&Om.pdf'
# tanh
file_name = 'pulse_al='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Spect&Om.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()


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

# file_name_PT = 'Om_'+str(coeff[0])+'t_alpha='+str(alpha3)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur[3])+'_PT.txt'

# full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

# oqupy.SimpleProcessTensor.export(self=process_tensor,
#                                   filename=full_path_PT)

# Read in PT
full_path_PT = 'C:/Users/sony/Desktop/Project/process_tensor/alpha=0.1_T=0_dt=0.05_dkmax=70_epsrel=1e-06_dur=5.0_PT.txt'

process_tensor = oqupy.process_tensor.FileProcessTensor(filename=full_path_PT,
                                                        mode = 'read')


# # Spectral density ohmic
def J(t):
    return 4*alpha*gaussian_shape(t)*np.exp(-2*gaussian_shape(t)/omega_cutoff)


# Define time dependent system
# (Time dependent system gaussian_shape(t, area = np.pi/2.0, tau = 0.245))
def hamiltonian_t(t):
    return gaussian_shape(t)*sigma_z

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
t0, s_z0 = dynamics.expectations(sigma_z, real=True)
#################################################################

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



# Compute dynamics
dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                  start_time = 0,dt=dtME,num_steps=num_stepsME)

tME0, s_zME0 = dynamics.expectations(sigma_z, real=True)

# Plot
# plt.plot(t1, s_z1, label=r'TEMPO $\Omega(t)$={0}t'.format(coeff[0]))
# How to scale gauss?? coeff = 1?
# plt.plot(t2, s_z2, label=r'TEMPO $\Omega(t)$ pulse'.format())
#plt.plot(t3, s_z3, label=r'TEMPO $\Omega$ = {0}'.format(Omega))

plt.plot(t0, s_z0, label=r'TEMPO $\Omega(t)$ pulse'.format())

plt.title(r'Decay from Up Density Matrix State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
plt.grid()

desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'pulse_al='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()

# Plot Zoom
plt.plot(t0, s_z0, label=r'TEMPO $\Omega(t)$ pulse'.format())

plt.title(r'Decay from Up Density Matrix State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
plt.grid()
plt.ylim([-1,-0.65])
plt.xlim([dur/2, dur])

desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'pulseZm_al='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()


# Plot QME
plt.plot(tME0, s_zME0,color='r', label='QME Pulse '.format())


plt.title(r'Decay from Up Density Matrix State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
plt.grid()

desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'pulseQME_al='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)
