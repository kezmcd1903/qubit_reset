"""

Solving Time-Dependent System Hamiltonians

- Generate/load process tensor
- Solves using PT-TEMPO and master equations
- Plots
- Loops through for many different time dependent energy splitting shapes

"""

import sys
import os
sys.path.insert(0,'..')

import oqupy
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})
from matplotlib.pyplot import cm
import timeit
import math

# Start timer to print the time it takes to run the full script
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

# Set initial state - excited state (represented by up density matrix) for two level system
init_st = up_density_matrix

# Pick system parameters
omega_cutoff = 5
alpha = 0.1
# Temperature in kelvin
T = 0 # K
# Temperature in units of 1/(ps kB). {1 K = 0.1309 1/(ps kB)}
temperature = T*0.1309   # 1/(ps kB)

dur = 4

# TEMPO parameters
dt=0.05
dkmax=200
epsrel=10**(-6)

# Time step for Master Equation
dtME = 0.1
num_stepsME = dur/dtME


##### Define time-dependent energy splitting 

# Linear dependence

# Define parameters to change lin ramp
# coeff = [1,2,3,4]
# shift = [1,2,3,4]

# def lin_ramp(t,coeff,shift):
#     return coeff*t + shift


# Gaussian profile

# Define parameters to change gaussian
# areas = [1,2,3,4]
# taus = [1,2,3,4]
# t_0s = [1,2,3,4]

# def gauss(t, area, tau, t_0): 
#     return area/(tau*np.sqrt(np.pi)) * np.exp(-(t-t_0)**2/(tau**2))



# Hyperbolic tanh profile

# Define parameters to change tanh shape
# Refered to as the grid
# grid can be extended to search as many combinations as necessary
# Set D = 0 if wanting to look at a time-independent spacing
Ds = [0,10,10,10]
As = [1,4,3.5,3]
bs = [1,10,9,8]
cs = [8.6,12.5,12.5,12.5]

# Depending on the results of these combinations the grids can then be cut so as to compare particular values using the following 
f = 0
l = 28

# Ds1 = Ds[11:18]
# As1 = As[11:18]
# bs1 = bs[11:18]
# cs1 = cs[11:18]

# Ds2 = Ds[25:27]
# As2 = As[25:27]
# bs2 = bs[25:27]
# cs2 = cs[25:27]

# Ds = Ds1 + Ds2
# As = As1 + As2
# bs = bs1 + bs2
# cs = cs1 + cs2

Ds = Ds[f:l]
As = As[f:l]
bs = bs[f:l]
cs = cs[f:l]


def tanh(t,A,b,c,D):
    return D*np.tanh(A*t-b)+c




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

# file_name_PT = 'alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_PT.txt'

# full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

# oqupy.SimpleProcessTensor.export(self=process_tensor,
#                                   filename=full_path_PT)

# Read in PT
full_path_PT = 'process_tensor/alpha=0.1_T=0_dt=0.05_dkmax=200_epsrel=1e-06_dur=30_PT.txt'

process_tensor = oqupy.process_tensor.FileProcessTensor(filename=full_path_PT,
                                                        mode = 'read')


# Define colour wheel to use for plots

# split viridis colour scheme evenly
n = int(len(As))#2
c = cm.viridis(np.linspace(0, 1, n))
#colors = np.concatenate((c, c))

# or use custom list
colors = ['limegreen','darkcyan','navy','purple']

# Changing line thickess can be desireable if plotting lots of different tanh shape 
# It can be difficult to distinguish which colour relates to which tanh.

#lines = [1,1,1,1]


# Plot all the different time dependent shapes
for A,b,c,D,color in zip(As,bs,cs,Ds,colors):#,line,lines
    #plt.title('Time Dependent Energy Splitting')
    t_plot1 = np.linspace(0,int(dur),int(dur)*10)
    # Plot tanh shapes
    plt.plot(t_plot1, tanh(t_plot1,A=A,b=b,c=c,D=D),c=color,label=r'${3} tanh({0}t-{1})+{2}$'.format(A,b,c,D))#,color=colour,linewidth=line,
    plt.xlabel("Time (ps)")
    plt.ylabel(r'Spacing Frequency $\Omega \;(\mathrm{ps}^{-1})$')

# Plot horizontal line signifying the peak in spectral density
plt.axhline(y = 2.5,color='black',linestyle='dashed',label=r'$\Omega_{p}=2.5$')
plt.legend()# Setting to remove legend from plot frameif there's a lot of values to plot - bbox_to_anchor=(1.1, 1.05)
# Save plots
desired_folder = 'Plots'
file_name = 'Split_tanh_alpha='+str(alpha)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path,bbox_inches="tight")
plt.show()
##########################

# Define each system relating to a different time-dependent shape
systems = []
for A,b,c,D in zip(As,bs,cs,Ds):
    def hamiltonian_t(t,A=A,b=b,c=c,D=D):
        return tanh(t,A=A,b=b,c=c,D=D)*sigma_z
    system = oqupy.TimeDependentSystem(hamiltonian_t)
    systems.append(system)


# Calculate dynamics for each system

# Defining empty lists to append results to
s_z_list = []
t_list = []
for system in systems:
    # Calculate
    dynamics = oqupy.compute_dynamics(
        process_tensor=process_tensor,
        system=system,
        initial_state=init_st,
        num_steps=dur/dt,
        start_time=0)
    # Saving results
    t, s_z = dynamics.expectations(sigma_z, real=True)
    # Appending to lists
    s_z_list.append(s_z)
    t_list.append(t)
    print(".", end="", flush=True)
print(" done.", flush=True)   

# Save Data
desired_folder_data = 'PlotData'

file_name_t = 'tanh_PltDat'+'a='+str(alpha)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'

file_name_s_z = 'tanh_PltDat'+'a='+str(alpha)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'

full_path_t = os.path.join(desired_folder_data, file_name_t)
full_path_s_z = os.path.join(desired_folder_data, file_name_s_z)

np.savetxt(full_path_t, t_list)
np.savetxt(full_path_s_z, s_z_list)



# # Spectral density in ohmic form
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

# Time dependent decay rates
Gamma = [ lambda t: gamma_1*J(t), lambda t: gamma_2*J(t)]

# Specifying Lindblad operators
lindblad_operators = [ lambda t: sigma_minus, lambda t: sigma_plus]


# Define each system relating to a different time-dependent shape
systemMEs = []
for A,b,c,D in zip(As,bs,cs,Ds):
    def hamiltonian_t(t,A=A,b=b,c=c,D=D):
        return tanh(t,A=A,b=b,c=c,D=D)*sigma_z
    systemME = oqupy.TimeDependentSystem(hamiltonian_t,
                                    gammas=Gamma,
                                    lindblad_operators=lindblad_operators)
    systemMEs.append(systemME)

    
# Solve dynamics with master equation   
s_zME_list = []
tME_list = []
for systemME in systemMEs:
    # Calculate
    dynamics = oqupy.compute_dynamics(system = systemME, initial_state = init_st,\
                                  start_time = 0,dt=dtME,num_steps=num_stepsME)
    # Save results
    tME, s_zME = dynamics.expectations(sigma_z, real=True)
    s_zME_list.append(s_zME)
    tME_list.append(tME)
    print(".", end="", flush=True)
print(" done.", flush=True)

#### Save data
file_name_tME = 'tanh_PltDat'+'a='+str(alpha)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_TimeME.txt'

file_name_s_zME = 'tanh_PltDat'+'a='+str(alpha)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_zME.txt'

full_path_tME = os.path.join(desired_folder_data, file_name_tME)
full_path_s_zME = os.path.join(desired_folder_data, file_name_s_zME)

np.savetxt(full_path_tME, tME_list)
np.savetxt(full_path_s_zME, s_zME_list)




######## Plot Solutions

# Plot decays given by PT-TEMPO
#,line,lines
for t, s_z, A, b, c, D,color in zip(t_list, s_z_list, As,bs,cs,Ds,colors):
    plt.title(r'PT-TEMPO Decay ($\alpha$={0})'.format(alpha)) 
    plt.plot(t, s_z,c=color, label=r'$\Omega(t)={3} tanh({0}t-{1})+{2}$'.format(A,b,c,D))#color=colour,
    plt.xlabel('Time (ps)')
    plt.ylabel(r'$<\sigma_{z}>$')

plt.legend()
# plt.grid()
# Save plots
desired_folder = 'Plots'
file_name = 'TEMPOdecay_tanh_alpha='+str(alpha)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path,bbox_inches="tight")
plt.show()



# Plot Fidelity Evolution
#,line,lines
for t, s_z, A, b, c, D,color in zip(t_list, s_z_list, As,bs,cs,Ds,colors): 
    #plt.title(r'PT-TEMPO Decay ($\alpha$={0})'.format(alpha)) 
    plt.plot(t, (((s_z-1)/(-2))*100),c=color,  label=r'$\Omega(t)={3} tanh({0}t-{1})+{2}$'.format(A,b,c,D))#linewidth=line,
    plt.xlabel('Time (ps)')
    plt.ylabel(r'$<\sigma_{z}>$')

# Plot horizontal line corresponding to the threshold fidelity
plt.axhline(y = 99,color='black',linestyle='dashed',
            label="99%")
#plt.title(r'Time Taken to Reach Fidelity at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel('Fidelity (%)')
plt.legend()
# plt.minorticks_on()
# plt.grid(which='both')
# Save plots
desired_folder = 'Plots'
file_name = 'FOM_tanh_alpha='+str(alpha)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path,bbox_inches="tight")
plt.show()


# Plot Fidelity Evolution Zoomed In
#,line,lines
for t, s_z, A, b, c, D,color in zip(t_list, s_z_list, As,bs,cs,Ds,colors): 
    #plt.title(r'TEMPO Decay ($\alpha$={0})'.format(alpha)) 
    plt.plot(t, (((s_z-1)/(-2))*100),c=color,  label=r'${3} tanh({0}t-{1})+{2}$'.format(A,b,c,D))#linewidth=line,\Omega(t)=
    plt.xlabel('Time (ps)')
    plt.ylabel(r'$<\sigma_{z}>$')


plt.axhline(y = 99,color='black',linestyle='dashed',
            label="99%")
#plt.title(r'Fidelity Evolution at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel('Fidelity (%)')
plt.legend()
# plt.minorticks_on()
# plt.grid(which='both')

# Select zoom range
plt.ylim([88,100])
plt.xlim([2, int(dur)])

# Save plots
desired_folder = 'Plots'
file_name = 'FOMzm_tanh_alpha='+str(alpha)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path,bbox_inches="tight")
plt.show()


# Master Equation Plots
#,line,lines
for t, s_z, A, b, c, D,color in zip(t_list, s_z_list, As,bs,cs,Ds,colors): 
    plt.title(r'Master Equation Decay ($\alpha$={0})'.format(alpha)) 
    plt.plot(tME, s_zME,c=color,  label=r'$\Omega(t)={3} tanh({0}t-{1})+{2}$'.format(A,b,c,D))#linewidth=line,
    plt.xlabel('Time (ps)')
    plt.ylabel(r'$<\sigma_{z}>$')

plt.legend()
# plt.grid()
# Save plots
desired_folder = 'Plots'
file_name = 'MEdecay_tanh_alpha='+str(alpha)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path,bbox_inches="tight")
plt.show()




stop = timeit.default_timer()
print('Time: ', stop - start)
