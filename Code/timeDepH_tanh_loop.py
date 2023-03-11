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

# Time step for Master Equation
dtME = 0.01
num_stepsME = dur/dtME


def tanh(t,A,b,c,D):
    return D*np.tanh(A*t-b)+c

As = [1, 2, 3, 3 ]
bs = [3, 6, 8, 10]
cs = [4,12.5,12.5,16.5]
Ds = [0,10,10,14 ]


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


for A,b,c,D in zip(As,bs,cs,Ds):
    plt.title('Time Dependent Energy Splitting')
    plt.axhline(y = 2.5,color='m',linestyle='dashed')
    t_plot1 = np.linspace(0,int(dur),int(dur)*10)
    plt.plot(t_plot1, tanh(t_plot1,A=A,b=b,c=c,D=D),label=r'$\Omega(t)={3} tanh({0}t-{1})+{2}$'.format(A,b,c,D))
    plt.xlabel("Time (ps)")
    plt.ylabel(r'$\Omega \;(\mathrm{ps}^{-1})$')
    
plt.legend()
# Save plots
desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'Split_tanh_alpha='+str(alpha)+'A='+str(As)+'b='+str(bs)+'c='+str(cs)+'D='+str(Ds)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()


##########################

systems = []
for A,b,c,D in zip(As,bs,cs,Ds):
    def hamiltonian_t(t,A=A,b=b,c=c,D=D):
        return tanh(t,A=A,b=b,c=c,D=D)*sigma_z
    system = oqupy.TimeDependentSystem(hamiltonian_t)
    systems.append(system)
########



s_z_list = []
t_list = []
for system in systems:
    dynamics = oqupy.compute_dynamics(
        process_tensor=process_tensor,
        system=system,
        initial_state=init_st,
        num_steps=dur/dt,
        start_time=0)
    t, s_z = dynamics.expectations(sigma_z, real=True)
    s_z_list.append(s_z)
    t_list.append(t)
    print(".", end="", flush=True)
print(" done.", flush=True)


############# Master equation
def tanh(t,A,b,c,D):
    return D*np.tanh(A*t-b)+c

# # Spectral density ohmic
def J(t):
    return 4*alpha*tanh(t,A=A,b=b,c=c,D=D)*np.exp(-2*tanh(t,A=A,b=b,c=c,D=D)/omega_cutoff)


# B-E Occupancy (for temp dependence in master equation)
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


systemMEs = []
for A,b,c,D in zip(As,bs,cs,Ds):
    def hamiltonian_t(t,A=A,b=b,c=c,D=D):
        return tanh(t,A=A,b=b,c=c,D=D)*sigma_z
    systemME = oqupy.TimeDependentSystem(hamiltonian_t,
                                    gammas=Gamma,
                                    lindblad_operators=lindblad_operators)
    systemMEs.append(systemME)

# Compute dynamics
s_zME_list = []
tME_list = []
for systemME in systemMEs:
    dynamics = oqupy.compute_dynamics(system = systemME, initial_state = init_st,\
                                  start_time = 0,dt=dtME,num_steps=num_stepsME)
    tME, s_zME = dynamics.expectations(sigma_z, real=True)
    s_zME_list.append(s_zME)
    tME_list.append(tME)
    print(".", end="", flush=True)
print(" done.", flush=True)
##################


for t, s_z, A, b, c, D in zip(t_list, s_z_list, As,bs,cs,Ds):
    plt.title(r'TEMPO Decay ($\alpha$={0})'.format(alpha)) 
    plt.plot(t, s_z, label=r'$\Omega(t)={3} tanh({0}t-{1})+{2}$'.format(A,b,c,D))
    plt.xlabel('Time (ps)')
    plt.ylabel(r'$<\sigma_{z}>$')

plt.legend()
plt.grid()
# Save plots
desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'TEMPOdecay_tanh_alpha='+str(alpha)+'A='+str(As)+'b='+str(bs)+'c='+str(cs)+'D='+str(Ds)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()


for t, s_z, A, b, c, D in zip(t_list, s_z_list, As,bs,cs,Ds):
    plt.title(r'TEMPO Decay ($\alpha$={0})'.format(alpha)) 
    plt.plot(t, s_z, label=r'$\Omega(t)={3} tanh({0}t-{1})+{2}$'.format(A,b,c,D))
    plt.xlabel('Time (ps)')
    plt.ylabel(r'$<\sigma_{z}>$')

plt.legend()
plt.grid()
plt.ylim([-1,-0.8])
plt.xlim([dur/2, dur])
# Save plots
desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'TEMPOdecayZM_tanh_alpha='+str(alpha)+'A='+str(As)+'b='+str(bs)+'c='+str(cs)+'D='+str(Ds)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()

# Plot FOMs
for t, s_z, A, b, c, D in zip(t_list, s_z_list, As,bs,cs,Ds):
    plt.title(r'TEMPO Decay ($\alpha$={0})'.format(alpha)) 
    plt.plot(t, (((s_z-1)/(-2))*100), label=r'$\Omega(t)={3} tanh({0}t-{1})+{2}$'.format(A,b,c,D))
    plt.xlabel('Time (ps)')
    plt.ylabel(r'$<\sigma_{z}>$')


plt.axhline(y = 99,color='c',linestyle='dashed',
            label="99%")
plt.title(r'Time Taken to Reach Fidelity at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel('Fidelity (%)')
plt.legend()
plt.minorticks_on()
plt.grid(which='both')
# Save plots
desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'FOM_tanh_alpha='+str(alpha)+'A='+str(As)+'b='+str(bs)+'c='+str(cs)+'D='+str(Ds)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()


# Plot FOMs zoom
for t, s_z, A, b, c, D in zip(t_list, s_z_list, As,bs,cs,Ds):
    plt.title(r'TEMPO Decay ($\alpha$={0})'.format(alpha)) 
    plt.plot(t, (((s_z-1)/(-2))*100), label=r'$\Omega(t)={3} tanh({0}t-{1})+{2}$'.format(A,b,c,D))
    plt.xlabel('Time (ps)')
    plt.ylabel(r'$<\sigma_{z}>$')


plt.axhline(y = 99,color='c',linestyle='dashed',
            label="99%")
plt.title(r'Time Taken to Reach Fidelity at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel('Fidelity (%)')
plt.legend()
plt.minorticks_on()
plt.grid(which='both')
plt.ylim([90,100])
plt.xlim([dur/2, dur])
# Save plots
desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'FOMzm_tanh_alpha='+str(alpha)+'A='+str(As)+'b='+str(bs)+'c='+str(cs)+'D='+str(Ds)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()


# Master equation plots
for tME, s_zME, A, b, c, D in zip(tME_list, s_zME_list, As,bs,cs,Ds):
    plt.title(r'Master Equation Decay ($\alpha$={0})'.format(alpha)) 
    plt.plot(tME, s_zME, label=r'$\Omega(t)={3} tanh({0}t-{1})+{2}$'.format(A,b,c,D))
    plt.xlabel('Time (ps)')
    plt.ylabel(r'$<\sigma_{z}>$')

plt.legend()
plt.grid()
# Save plots
desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'MEdecay_tanh_alpha='+str(alpha)+'A='+str(As)+'b='+str(bs)+'c='+str(cs)+'D='+str(Ds)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()




stop = timeit.default_timer()
print('Time: ', stop - start)
