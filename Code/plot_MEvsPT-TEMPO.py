"""

Solves Time Independent Systems via PT-TEMPO and Master Equation Methods for Comparison

- Calculations carried out a specific coupling and energy spacing
- Steps follow the same procedure as 'difOmeags.py' until plotting
- Generate/load process tensor
- Solves using PT-TEMPO and master equation solver (which uses the TEMPO method)
- Plots both for comparison

"""

import sys
import os
sys.path.insert(0,'..')

import oqupy
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})
from matplotlib.pyplot import cm
import timeit
import math
from scipy.fft import fft, fftshift, fftfreq
import multiprocessing as mp
import ipyparallel as ipp

# Start timer
start = timeit.default_timer()

# Set up parallel processing
cluster = ipp.Cluster(n=12)
cluster


# Define operators
sigma_x = oqupy.operators.sigma("x")
sigma_y = oqupy.operators.sigma("y")
sigma_z = oqupy.operators.sigma("z")
sigma_plus = oqupy.operators.sigma("+")
sigma_minus = oqupy.operators.sigma("-")

up_density_matrix = oqupy.operators.spin_dm("z+")
down_density_matrix = oqupy.operators.spin_dm("z-")
mixed_density_matrix = oqupy.operators.spin_dm("mixed")

# Set excited state as intial state
init_st = up_density_matrix

# Pick System parameters
fid_thres = 99
Om_spacing = .1
start_om = 4
end_om = 9
# Omega = np.arange(start_om,end_om+Om_spacing,Om_spacing)
Omega = np.array([1])
# print('Omega =',Omega)
omega_cutoff = 5.0
alpha = 0.1
# Temperature in kelvin
T = 0 # K
# Temperature in units of 1/(ps kB)
# 1 K = 0.1309 1/(ps kB)
temperature = T*0.1309 # 1/(ps kB)

# TEMPO parameters
dt=0.05
dkmax=200 #200
epsrel=1e-06
# duration (ps)
dur = 5

# Time steps for Master Equation
dtME = 0.01#0.05
num_stepsME = dur/dtME

# Save TEMPO parameters
tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
                                          epsrel=epsrel)


# B-E Occupancy for master equation
def N(T):
    if T == 0:
        N = 0
    else:
        N = (1/(math.exp(((1.055e-34)*(10**12))/((1.381e-23)*T))-1))
    return N

#n = int(len(Omega))
#c = cm.viridis(np.linspace(0, 1, n))

# Define variable names that are getting iterated over
t=[]
s_z=[]
tME=[]
s_zME=[]
J=[]
fid=[]
t_fidThres=[]
for i in range(len(Omega)):
    t.append('t'+str(i))
    s_z.append('s_z'+str(i))
    tME.append('tME'+str(i))
    s_zME.append('s_zME'+str(i))
    J.append('J'+str(i))
    fid.append('FOM'+str(i))
    t_fidThres.append('t_fidThres'+str(i))
    

# Specify correlations
correlations = oqupy.PowerLawSD(alpha=alpha,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='exponential',
                                temperature=(temperature))


bath = oqupy.Bath(sigma_x, correlations)


# # #############
#Compute process tensor
process_tensor = oqupy.pt_tempo_compute(bath=bath,
                                        start_time=0.0,
                                        end_time=dur,
                                        parameters=tempo_parameters)

# ###### Save PT ######## Doesn't work if already saved
# desired_folder_PT = 'process_tensor'

# file_name_PT = 'alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_PT.txt'

# full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

# oqupy.SimpleProcessTensor.export(self=process_tensor,
#                                     filename=full_path_PT)
####################################

# # # #Read in PT
# full_path_PT = 'process_tensor/alpha=0.1_T=0_dt=0.05_dkmax=200_epsrel=1e-06_dur=5_PT.txt'

# process_tensor = oqupy.process_tensor.FileProcessTensor(filename=full_path_PT,
#                                                         mode = 'read')

#############################################
##############################################    
    
    
system = oqupy.System(Omega[i]*sigma_z)
#     # Compute
dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
                                system=system,
                                initial_state=init_st,
                                start_time=0.0 # needed?, progress_type
                                ,num_steps=dur/dt)

t0, s_z0 = dynamics.expectations(sigma_z, real=True)

# Calculating system dynamics. Here, PT-TEMPO caculations are commented out so as to only solve master equation dynamics

for i in range(len(Omega)):
    # Specify system
    system = oqupy.System(Omega[i]*sigma_z)
    
# #     # Compute
#     dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
#                                     system=system,
#                                     initial_state=init_st,
#                                     start_time=0.0 # needed?, progress_type
#                                     ,num_steps=dur/dt)
    
#     t[i], s_z[i] = dynamics.expectations(sigma_z, real=True)
    
#     # Save
    desired_folder_data = 'PlotData'
    
#     file_name_t = 'alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
#         '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'
        
#     file_name_s_z = 'alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
#         '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'
        
#     full_path_t = os.path.join(desired_folder_data, file_name_t)
#     full_path_s_z = os.path.join(desired_folder_data, file_name_s_z)
    
#     np.savetxt(full_path_t, t[i])
#     np.savetxt(full_path_s_z, s_z[i])
    
#     ###########
    
    J[i] = 4*alpha*Omega[i]*np.exp(-2*Omega[i]/omega_cutoff)
    
    # Derived master equation
    #gamma1 = gamma*(N+1)
    gamma_1 = 2*J[i]*np.pi*(N(T)+1)#*(N+1)
    # gamma2 = gamma*(N)
    gamma_2 = 2*J[i]*np.pi*N(T) #*N
    
    Gamma = [gamma_1, gamma_2]
    
    # This is in the same form as in notes
    lindblad_operators  = [sigma_minus, sigma_plus]
        
    
    systemD = oqupy.System(
              Omega[i] * sigma_z,
            gammas=Gamma,
            lindblad_operators=lindblad_operators)
    
    # Compute
    dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                      start_time = 0,dt=dtME,num_steps=num_stepsME)
    
    tME[i], s_zME[i] = dynamics.expectations(sigma_z, real=True)
    
    #save
    file_name_tME = 'alpha='+str(alpha)+'_T=0_dt=0.05_dkmax=70_epsrel=1e-06_dur='+str(dur)+'_TimeME.txt'
        
    file_name_s_zME = 'alpha='+str(alpha)+'_T=0_dt=0.05_dkmax=70_epsrel=1e-06_dur='+str(dur)+'_S_zME.txt'
        
    full_path_t = os.path.join(desired_folder_data, file_name_tME)
    full_path_s_z = os.path.join(desired_folder_data, file_name_s_zME)
    
    np.savetxt(full_path_t, tME[i])
    np.savetxt(full_path_s_z, s_zME[i])



#####################################

# Read in Data not already calculated

#####################################

# Load data
folder_data = 'PlotData'


file_name_t = 'alpha='+str(alpha)+'_T=0_dt=0.05_dkmax=70_epsrel=1e-06_dur=2000.0_Time.txt'

file_name_s_z = 'alpha='+str(alpha)+'_T=0_dt=0.05_dkmax=70_epsrel=1e-06_dur=2000.0_S_z.txt'


# file_name_t = 'alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
#         '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'

# file_name_s_z = 'alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
#         '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'


# full_path_t = os.path.join(folder_data, file_name_t)
# full_path_s_z = os.path.join(folder_data, file_name_s_z)

# t0 = np.loadtxt(full_path_t)
# s_z0 = np.loadtxt(full_path_s_z)


file_name_tME = 'alpha='+str(alpha)+'_T=0_dt=0.05_dkmax=70_epsrel=1e-06_dur='+str(dur)+'_TimeME.txt'

file_name_s_zME = 'alpha='+str(alpha)+'_T=0_dt=0.05_dkmax=70_epsrel=1e-06_dur='+str(dur)+'_S_zME.txt'


full_path_tME = os.path.join(folder_data, file_name_tME)
full_path_s_zME = os.path.join(folder_data, file_name_s_zME)

tME0 = np.loadtxt(full_path_tME)
s_zME0 = np.loadtxt(full_path_s_zME)


# ######################################################


## PLOT for comparison

c = ['limegreen','midnightblue']#,'purple']

# Plot TEMPO and Master Equation Solutions

plt.plot(t0, s_z0,color=c[0],label='PT-TEMPO')
plt.plot(tME0, s_zME0,color=c[1], label='ME')

#plt.title(r'Decay from Excited State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
# plt.grid()

desired_folder = 'Plots'
file_name = 'constOmsTEMPO_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)


plt.subplots_adjust(bottom=0.15)


plt.savefig(full_path,bbox_inches="tight")

plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)
