"""

Solves Time Independent Systems with Varying System Hamiltonians

- Generate/load process tensor
- Solves using PT-TEMPO and master equation solver (which uses the TEMPO method)
- Plots both
- Can repeat for many differeny system with varying system parameter
- Here the system parameter changed is the energy spacing frequency

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
from scipy.fft import fft, fftshift, fftfreq
import multiprocessing as mp
import ipyparallel as ipp
import mpmath as mp

# Start a timer which times the total time taken to run script and prints this time after completion
start = timeit.default_timer()

# To use parallel computation set to the number of core avialable to your machine
# For JuypterHub 12 cores, n = 12. For Cirrus cluster computer n = 18.
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

# Initial state is the excited state represented by the up density matrix
init_st = up_density_matrix

# Fidelity threshold value at which the time taken to reach this value is recorded and used in OptimalEnegySpacing_final_plot.py
fid_thres = 99

# Set energy spacing frequencies to sample
Om_spacing = 2
start_om = 1
end_om = 12
#Omega = np.arange(start_om,end_om+Om_spacing,Om_spacing)
Omega = np.array([2.5])
# print('Omega =',Omega)

# Set environment and interaction parameters

# Set cutoff frequency which sets the shape of the spectral density
omega_cutoff = 5.0
# Set coupling strength alpha
alpha = 0.05
# Set temperature in kelvin
T = 0 # K
# Convert temperature to units of 1/(ps kB)
# 1 K = 0.1309 1/(ps kB)
temperature = T*0.1309 # 1/(ps kB)


# Set TEMPO parameters
dt=0.05
dkmax=200 
epsrel=1e-06
# Set duration of system dynamics. In (picoseconds).
dur = 4

# Set time step for Master Equation
dtME = .05
# Calculate number of steps required for master equation dynamics
num_stepsME = dur/dtME

# Save TEMPO parameters
tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
                                          epsrel=epsrel)


# Define B-E Occupancy for master equation
def N(T):
    if T == 0:
        N = 0
    else:
        N = (1/(mp.exp(((1.055e-34)*(10**12))/((1.381e-23)*T))-1))
    return N


# # Make plot of spectral density function. Energy spacing frequencies are shown as vertical lines

# plt.title(r'Sampling the Spectral Density Function')

# Set frequency range over which to plot 
w = np.linspace(0,20,200)
# # Set spectral density function - here of ohmic form
J = 4*alpha*w*np.exp(-2*w/omega_cutoff)
# Plot
plt.plot(w,J,color='black')

# Define colour wheel for different colours of energy spacing frequency
n = int(len(Omega))
c = cm.viridis(np.linspace(0, 1, n))

# Plot energy spacing frequencies as vertical lines
for i in range(len(Omega)):
    plt.axvline(x = Omega[i],color='gold',linestyle='dashed',label=r'$\Omega_{p}=2.5$')

plt.xlabel(r'Energy Spacing Frequency $\Omega \; (ps^{-1})$')
plt.ylabel(r"Spectral Density $J(\Omega) \; (\mathrm{ps}^{-1})$")
plt.legend()
# Save plots
desired_folder = 'Plots'
file_name = 'Spect&Om_constOm_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path,bbox_inches="tight") # Specifying setting which ensures y-label isn't chopped off in saved pdf
plt.show()

# Save plot data for replotting
desired_folder_data = 'PlotData'

file_name_J = 'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.txt'

file_name_w = 'w='+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.txt'    

full_path_J = os.path.join(desired_folder_data, file_name_J)
full_path_w = os.path.join(desired_folder_data, file_name_w)

np.savetxt(full_path_J, J)
np.savetxt(full_path_w, w)
####


# Save system-environment correlations
correlations = oqupy.PowerLawSD(alpha=alpha,
                                zeta=1,
                                cutoff=omega_cutoff,
                                # Set cutoff type and therefore the spectral density. The profile the spectral density function follows
                                # after the cutoff. 
                                # Exponential is used for ohmic spectral density.
                                cutoff_type='exponential',
                                temperature=(temperature))

# Specify the environment bath
# takes the coefficient of the interaction Hamiltonian and the system-environment correlations
bath = oqupy.Bath(sigma_x, correlations)

## Generate process tensor (PT)

# ####################################
# Compute PT
# process_tensor = oqupy.pt_tempo_compute(bath=bath,
#                                         start_time=0.0,
#                                         end_time=dur,
#                                         parameters=tempo_parameters)

## Save PT. Retruns error if the exact PT (with the same name) has already been saved

# desired_folder_PT = 'C:/Users/sony/Desktop/Project/process_tensor/'

# file_name_PT = 'alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_PT.txt'

# full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

# oqupy.SimpleProcessTensor.export(self=process_tensor,
#                                     filename=full_path_PT)
# ####################################

#############################################
# #Read in PT
full_path_PT = 'process_tensor/alpha=0.05_T=5_dt=0.05_dkmax=200_epsrel=1e-06_dur=15_PT.txt'

process_tensor = oqupy.process_tensor.FileProcessTensor(filename=full_path_PT,
                                                        mode = 'read')
##############################################

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
    

# Compute PT-TEMPO and ME dynamics for each system
for i in range(len(Omega)):
    
    # Specify system
    system = oqupy.System(Omega[i]*sigma_z)
    
    # Compute
    dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
                                    system=system,
                                    initial_state=init_st,
                                    start_time=0.0 
                                    ,num_steps=dur/dt)
    
    # Save output of computation as variables
    t[i], s_z[i] = dynamics.expectations(sigma_z, real=True)
    
    # Save the outputs
    desired_folder_data = 'PlotData'
    
    file_name_t = 'PltDat_constOm='+str(Omega[i])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
        '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'
        
    file_name_s_z = 'PltDat_constOm='+str(Omega[i])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
        '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'
        
    full_path_t = os.path.join(desired_folder_data, file_name_t)
    full_path_s_z = os.path.join(desired_folder_data, file_name_s_z)
    
    np.savetxt(full_path_t, t[i])
    np.savetxt(full_path_s_z, s_z[i])
    
    ###########
    
    # Set spectral density (ohmic) for master equation dynamics
    J[i] = 4*alpha*Omega[i]*np.exp(-2*Omega[i]/omega_cutoff)
    
    # Define decay rates in the derived master equation
    #gamma1 = gamma*(N+1)
    gamma_1 = 2*J[i]*np.pi*(N(T)+1)#*(N+1)
    # gamma2 = gamma*(N)
    gamma_2 = 2*J[i]*np.pi*N(T) #*N
    
    Gamma = [gamma_1, gamma_2]
    
    # Specify the first two operators in each line of the master equation in lindblad form 
    # Details in (Section 3.2 of report and OQuPy documentation)
    lindblad_operators  = [sigma_minus, sigma_plus]
        
    # Save each system
    systemD = oqupy.System(
              Omega[i] * sigma_z,
            gammas=Gamma,
            lindblad_operators=lindblad_operators)
    
    # Compute
    dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                      start_time = 0,dt=dtME,num_steps=num_stepsME)
    
    tME[i], s_zME[i] = dynamics.expectations(sigma_z, real=True)
    
    # Save
    file_name_tME = 'PltDat_constOm='+str(Omega[i])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
        '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_TimeME.txt'
        
    file_name_s_zME = 'PltDat_constOm='+str(Omega[i])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
        '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_zME.txt'
        
    full_path_t = os.path.join(desired_folder_data, file_name_tME)
    full_path_s_z = os.path.join(desired_folder_data, file_name_s_zME)
    
    np.savetxt(full_path_t, tME[i])
    np.savetxt(full_path_s_z, s_zME[i])

#####################################################################
####################################################################


####### PLOTS ########

# Plot PT-TEMPO Solution
for i in range(len(Omega)):
    plt.plot(t[i], s_z[i],color=c[i],label=r'$\Omega_{0}={1}$ '.format(i,Omega[i]))

plt.title(r'PT-TEMPO - Decay from Excited State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
# plt.grid()

# Save plot
desired_folder = 'Plots'
file_name = 'constOmsTEMPO_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path,bbox_inches="tight")

plt.show()



# Plot PT-TEMPO solution zoomed in
for i in range(len(Omega)):
    plt.plot(t[i], s_z[i],color=c[i],label=r'$\Omega_{0}={1}$ '.format(i,Omega[i]))
    

plt.title(r'PT-TEMPO - Decay from Excited State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
# plt.minorticks_on()
# plt.grid(which='both')

# Specify which part you want to zoom in on
plt.ylim([-1,-0.7])
plt.xlim([dur/2, dur])

desired_folder = 'Plots'
file_name = 'constOmsTEMPOzoom_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path,bbox_inches="tight")

plt.show()


# Now calculate the fourier transform of the decay
# This is interesting at strong coupling alpha = 1 when oscillations are observed in the decay.
# Not so interesting for weaker coupling

# Save variable names
s_zft=[]
f =[]
for i in range(len(Omega)):
    s_zft.append('s_zft'+str(i))
    f.append('f'+str(i))
    
# Fourier transform
for i in range(len(Omega)):
    s_zft[i] = fftshift(fft(fftshift(s_z[i])))

n = t[0].size
timestep = dt

# Create a masking filter that removes the central peak in the FT due to the DC shift
filt = []
#print('length=',len(s_z0ft))

# Adjust the range of the loop to change the width of this filter
for i in range(int(len(s_zft[0])*(4/10))):
    filt.append(1)
for i in range(int(len(s_zft[0])*(2/10)+1)):
    filt.append(0)
for i in range(int(len(s_zft[0])*(4/10))):
    filt.append(1)

# Now apply mask
for i in range(len(Omega)):
    s_zft[i] = filt*s_zft[i]
    # Shift each FT to center
    f[i] = fftshift(fftfreq(n, d=timestep))
    
# Plot FT
for i in range(len(Omega)):
    plt.plot(f[i], np.absolute(s_zft[i]),color=c[i],label=r'$\Omega_{0}={1}$ '.format(i,Omega[i]))

# Plot black dashed vertical lines that relate to the frequency of observed oscillations as predicted by the system eigenvalues
for i in range(len(Omega)):
    plt.axvline(x = (np.pi/Omega[i])**-1,color='black',linestyle='dashed')    
    
plt.title(r'Fourier Transform at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel(r'Frequency $ \;(ps^{-1})$')
plt.ylabel(r'FT of $\langle\sigma_z\rangle$')
plt.legend()
# plt.grid()

desired_folder = 'Plots'
file_name = 'constOmsTEMPOft_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path,bbox_inches="tight")

plt.show()



# Plot FT zoomed in
for i in range(len(Omega)):
    plt.plot(f[i], np.absolute(s_zft[i]),color=c[i],label=r'$\Omega_{0}={1}$ '.format(i,Omega[i]))

for i in range(len(Omega)):
    plt.axvline(x = (np.pi/Omega[i])**-1,color='black',linestyle='dashed')    
    
plt.title(r'Fourier Transform of Decay Oscillations at $\alpha$={0}'.format(alpha)) 
plt.xlabel(r'Frequency $ \;(ps^{-1})$')
plt.ylabel(r'FT of $\langle\sigma_z\rangle$')
plt.legend()

# Specify where you want to zoom
plt.ylim([0.5,2])
plt.xlim([4, 8])
# plt.grid()

desired_folder = 'Plots'
file_name = 'constOmsTEMPOft_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path,bbox_inches="tight")

plt.show()



# Plot Fidelity Evolution
for i in range(len(Omega)):
    # Calculate fidelity from the sigma_z expectation value
    fid[i] = (((s_z[i]-1)/(-2))*100)
    # plot
    plt.plot(t[i], fid[i],color=c[i],label=r'$\Omega_{0}={1}$ '.format(i,Omega[i]))

# Plot threshold fidelity as black dashed horizontal line
plt.axhline(y = 99,color='black',linestyle='dashed',
            label="99%")

plt.title(r'PT-TEMPO - Fidelity Evolution at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel('Fidelity (%)')
plt.legend()
# plt.minorticks_on()
# plt.grid(which='both')

desired_folder = 'Plots'
file_name = 'FOMs_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path,bbox_inches="tight")

plt.show()



# Plot Fidelity Evolution zoomed in
for i in range(len(Omega)):
    plt.plot(t[i], fid[i],color=c[i],label=r'TEMPO $\Omega_{0}={1}$ '.format(i,Omega[i]))

plt.axhline(y = 99,color='black',linestyle='dashed',
            label="99%")


plt.title(r'Time Taken to Reach Fidelity at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel('Fidelity (%)')
plt.legend() #loc='lower left'
# plt.minorticks_on()
# plt.grid(which='both')

# Zoom limits
plt.ylim([92,100])
plt.xlim([dur/2, dur])

desired_folder = 'Plots'
file_name = 'FOMsZm_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path,bbox_inches="tight")
plt.show()



# Plot master equation dynamics
for i in range(len(Omega)):
    plt.plot(tME[i], s_zME[i],color=c[i], label=r'$\Omega_{0}={1}$ '.format(i,Omega[i]))

plt.title(r'Master Equation - Decay from Excited State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
# plt.grid()

desired_folder = 'Plots'
file_name = 'constOmsQME_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path,bbox_inches="tight")

plt.show()


#########################################
# Save data for final energy spacing plot
# where each system evolution reaches the threshold fidelity is saved for further plotting at different coupling strengths


for i in range(len(Omega)):
    # If there's a value over 99 save the time at which this first occurs
    if max(fid[i]) > fid_thres:
        # Find the first value where the threshold fidelity is achieved
        idx_fidThres = next(idx for idx, value in enumerate(fid[i]) if value > fid_thres)
        # The time this freshold fidelity is first reached is
        tCur = t[i]
        # Save the time when fidelity first over 99
        t_fidThres[i] = tCur[idx_fidThres]
    
    # If fidelity isn't reached
    else:
        # set mask value to remove at later stage
        t_fidThres[i] = -999



Omega = Omega.tolist()
# For debugging
# print(Omega)
# print(t_fidThres)#.remove(-999)

# Removing the energy splitting frequencies that don't achieve the threshold fidelity

# # Entering the indices at which items are to be deleted
ind2rem =[]
for i in range(len(t_fidThres)):
    if t_fidThres[i] == -999:
        ind2rem.append(i)
    else:
        continue
    
# # Reversing Indices List
indicesList = sorted(ind2rem, reverse=True)

# Traversing in the indices list to delete each item
for indx in indicesList:

    # checking whether the corresponding iterator index is less than the list length
    if indx < len(t_fidThres):

      # removing element by index using pop() function
      t_fidThres.pop(indx)
      Omega.pop(indx)

# # printing the list after deleting items at the given indices
print("Omega after deleting items at the given indices:\n", Omega)
print("t_fidThres after deleting items at the given indices:\n", t_fidThres)


#  Save data for final FOM plot

desired_folder_data = 'FOMdata'

file_name_Omega = 'fid_thres'+str(fid_thres)+'a='+str(alpha)+'_Omega.txt'
    
file_name_t_fidThres = 'fid_thres'+str(fid_thres)+'a='+str(alpha)+'_t_fidThres.txt'
    
full_path_Omega = os.path.join(desired_folder_data, file_name_Omega)
full_path_t_fidThres = os.path.join(desired_folder_data, file_name_t_fidThres)

np.savetxt(full_path_Omega, Omega)
np.savetxt(full_path_t_fidThres, t_fidThres)



# Stop timer
stop = timeit.default_timer()
print('Time: ', stop - start)
