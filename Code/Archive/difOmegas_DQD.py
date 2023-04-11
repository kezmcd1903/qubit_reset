"""

Solves Time Independent Systems with Varying System Hamiltonians

- Input process tensor
- Solves using TEMPO and master equation 
- Plots both
- Can repeat for 4 different systems

"""

import sys
import os
sys.path.insert(0,'..')

import oqupy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import timeit
import math
from scipy.fft import fft, fftshift, fftfreq
import multiprocessing as mp
import ipyparallel as ipp
import mpmath as mp

start = timeit.default_timer()

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

# initial state = up_density_matrix or mixed_density_matrix
init_st = up_density_matrix

# Pick System parameters
fid_thres = 99
Om_spacing = 0.05
start_om = 0.05
end_om = 1
#Omega = np.arange(start_om,end_om+Om_spacing,Om_spacing)
Omega = np.array([5.2])
# print('Omega =',Omega)
omega_cutoff = 5.0
alpha = 0.05
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
dur = 40


# Time steps for Master Equation
dtME = .05#0.05
num_stepsME = dur/dtME


tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
                                          epsrel=epsrel)


# B-E Occupancy for master equation
def N(T):
    if T == 0:
        N = 0
    else:
        N = (1/(mp.exp(((1.055e-34)*(10**12))/((1.381e-23)*T))-1))
    return N

# # Make plot of spectral density and time dependent splitting
plt.title(r'Spectral Density as a Function of Energy Splitting')
w = np.linspace(0,20,200)#20,200

# # Spectral density ohmic
J = 4*alpha*w*np.exp(-2*w/omega_cutoff)

plt.plot(w,J,color='black')#.format(coeff) tau = 0.245


n = int(len(Omega))
c = cm.viridis(np.linspace(0, 1, n))

for i in range(len(Omega)):
    plt.axvline(x = Omega[i],color=c[i],linestyle='dashed')


plt.xlabel(r'$\Omega \; (ps^{-1})$')
plt.ylabel(r"$J(\Omega) \; (\mathrm{ps}^{-1})$")

#plt.ylim([0, Omega[0]+0.1])
# Save plots
desired_folder = 'Plots'
file_name = 'DQD_Spect&Om_constOm_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path,bbox_inches="tight")
plt.show()

# Save
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


# Specify correlations
correlations = oqupy.PowerLawSD(alpha=alpha,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='exponential',
                                temperature=(temperature))


bath = oqupy.Bath(sigma_x, correlations)

# # #############
# process_tensor = oqupy.pt_tempo_compute(bath=bath,
#                                         start_time=0.0,
#                                         end_time=dur,
#                                         parameters=tempo_parameters)

# ###### Save PT ######## Doesn't work if already saved
# desired_folder_PT = 'C:/Users/sony/Desktop/Project/process_tensor/'

# file_name_PT = 'alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_PT.txt'

# full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

# oqupy.SimpleProcessTensor.export(self=process_tensor,
#                                     filename=full_path_PT)
""
# #Read in PT
full_path_PT = 'process_tensor/alpha=0.05_T=0_dt=0.05_dkmax=200_epsrel=1e-06_dur=40_PT.txt'

process_tensor = oqupy.process_tensor.FileProcessTensor(filename=full_path_PT,
                                                        mode = 'read')

#############################################
# #############################################
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
    



for i in range(len(Omega)):
    # Specify system
    system = oqupy.System(Omega[i]*sigma_z+2*sigma_x)
    
    # Compute
    dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
                                    system=system,
                                    initial_state=init_st,
                                    start_time=0.0 # needed?, progress_type
                                    ,num_steps=dur/dt)
    
    t[i], s_z[i] = dynamics.expectations(sigma_z, real=True)
    
    # Save
    desired_folder_data = 'PlotData'
    
    file_name_t = 'DQD_PltDat_constOm='+str(Omega[i])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
        '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'
        
    file_name_s_z = 'DQD_PltDat_constOm='+str(Omega[i])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
        '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'
        
    full_path_t = os.path.join(desired_folder_data, file_name_t)
    full_path_s_z = os.path.join(desired_folder_data, file_name_s_z)
    
    np.savetxt(full_path_t, t[i])
    np.savetxt(full_path_s_z, s_z[i])
    
    ###########
    
    # Spectral density ohmic
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
              Omega[i] * sigma_z+2*sigma_x,
            gammas=Gamma,
            lindblad_operators=lindblad_operators)
    
    # Compute
    dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                      start_time = 0,dt=dtME,num_steps=num_stepsME)
    
    tME[i], s_zME[i] = dynamics.expectations(sigma_z, real=True)
    
    # Save
    file_name_tME = 'DQD_PltDat_constOm='+str(Omega[i])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
        '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_TimeME.txt'
        
    file_name_s_zME = 'DQD_PltDat_constOm='+str(Omega[i])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
        '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_zME.txt'
        
    full_path_t = os.path.join(desired_folder_data, file_name_tME)
    full_path_s_z = os.path.join(desired_folder_data, file_name_s_zME)
    
    np.savetxt(full_path_t, tME[i])
    np.savetxt(full_path_s_z, s_zME[i])

#####################################################################
# ###################################################################
# ######################################################

# Plot TEMPO Solution
for i in range(len(Omega)):
    plt.plot(t[i], s_z[i],color=c[i],label=r'TEMPO $\Omega_{0}={1}$ '.format(i,Omega[i]))

plt.title(r'Decay from Up Density Matrix State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
# plt.grid()

desired_folder = 'Plots'
file_name = 'DQD_constOmsTEMPO_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()



# Plot TEMPO solution zoomed in
for i in range(len(Omega)):
    plt.plot(t[i], s_z[i],color=c[i],label=r'TEMPO $\Omega_{0}={1}$ '.format(i,Omega[i]))
    

plt.title(r'Decay from Up Density Matrix State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
# plt.minorticks_on()
# plt.grid(which='both')
plt.ylim([-1,-0.7])
plt.xlim([dur/2, dur])

desired_folder = 'Plots'
file_name = 'DQD_constOmsTEMPOzoom_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()



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

filt = []
#print('length=',len(s_z0ft))
for i in range(int(len(s_zft[0])*(4/10))):
    filt.append(1)
for i in range(int(len(s_zft[0])*(2/10)+1)):
    filt.append(0)
for i in range(int(len(s_zft[0])*(4/10))):
    filt.append(1)

 
for i in range(len(Omega)):
    s_zft[i] = filt*s_zft[i]
    # shift
    f[i] = fftshift(fftfreq(n, d=timestep))
    
# Plot
for i in range(len(Omega)):
    plt.plot(f[i], np.absolute(s_zft[i]),color=c[i],label=r'TEMPO $\Omega_{0}={1}$ '.format(i,Omega[i]))

for i in range(len(Omega)):
    plt.axvline(x = (np.pi/Omega[i])**-1,color='black',linestyle='dashed')    
    
plt.title(r'Fourier Transform at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel(r'Frequency $ \;(ps^{-1})$')
plt.ylabel(r'FT of $\langle\sigma_z\rangle$')
plt.legend()
# plt.grid()

desired_folder = 'Plots'
file_name = 'DQD_constOmsTEMPOft_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()



# Plot FT zoomed in
for i in range(len(Omega)):
    plt.plot(f[i], np.absolute(s_zft[i]),color=c[i],label=r'TEMPO $\Omega_{0}={1}$ '.format(i,Omega[i]))

for i in range(len(Omega)):
    plt.axvline(x = (np.pi/Omega[i])**-1,color='black',linestyle='dashed')    
    
plt.title(r'Fourier Transform at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel(r'Frequency $ \;(ps^{-1})$')
plt.ylabel(r'FT of $\langle\sigma_z\rangle$')
plt.legend()
plt.ylim([0.5,2])
plt.xlim([2.5, 10])
# plt.grid()

desired_folder = 'Plots'
file_name = 'DQD_constOmsTEMPOft_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()







# Plot FOMs
for i in range(len(Omega)):
    fid[i] = (((s_z[i]-1)/(-2))*100)
    plt.plot(t[i], fid[i],color=c[i],label=r'TEMPO $\Omega_{0}={1}$ '.format(i,Omega[i]))

plt.axhline(y = 99,color='black',linestyle='dashed',
            label="99%")


plt.title(r'Time Taken to Reach Fidelity at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel('Fidelity (%)')
plt.legend()
# plt.minorticks_on()
# plt.grid(which='both')

desired_folder = 'Plots'
file_name = 'DQD_FOMs_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()


# Make example FOM  plot with current omega values


# Plot FOMs zoomed in
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
plt.ylim([92,100])
plt.xlim([dur/2, dur])

desired_folder = 'Plots'
file_name = 'DQD_FOMsZm_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()



# Plot master equation solution
for i in range(len(Omega)):
    plt.plot(tME[i], s_zME[i],color=c[i], label=r'QME $\Omega_{0}={1}$ '.format(i,Omega[i]))

plt.title(r'Decay from Up Density Matrix State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
# plt.grid()

desired_folder = 'Plots'
file_name = 'DQD_constOmsQME_a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()

#########################################
# Save data for final FOM plot

# When = 99 (Fidelity Threshold)
for i in range(len(Omega)):
    # if there's a value over 99 do
    if max(fid[i]) > fid_thres:
        idx_fidThres = next(idx for idx, value in enumerate(fid[i]) if value > fid_thres)
        # Current t
        tCur = t[i]
        # Define time when fidelity first over 99
        t_fidThres[i] = tCur[idx_fidThres]
    
    else:
        # set mask value to remove at later stage
        t_fidThres[i] = -999



Omega = Omega.tolist()
# print(Omega)
# print(t_fidThres)#.remove(-999)

# # Entering the indices list at which the items are to be deleted
ind2rem =[]
for i in range(len(t_fidThres)):
    if t_fidThres[i] == -999:
        ind2rem.append(i)
    else:
        continue
    
# print(ind2rem)



# # Reversing Indices List
indicesList = sorted(ind2rem, reverse=True)

# Traversing in the indices list
for indx in indicesList:

    # checking whether the corresponding iterator index is less than the list length
    if indx < len(t_fidThres):

      # removing element by index using pop() function
      t_fidThres.pop(indx)
      Omega.pop(indx)

# # printing the list after deleting items at the given indices
print("Om after deleting items at the given indices:\n", Omega)
print("t_fidThres after deleting items at the given indices:\n", t_fidThres)


#  Save data for final FOM plot
 # Save
desired_folder_data = 'FOMdata'

file_name_Omega = 'fid_thres'+str(fid_thres)+'a='+str(alpha)+'_Omega.txt'
    
file_name_t_fidThres = 'fid_thres'+str(fid_thres)+'a='+str(alpha)+'_t_fidThres.txt'
    
full_path_Omega = os.path.join(desired_folder_data, file_name_Omega)
full_path_t_fidThres = os.path.join(desired_folder_data, file_name_t_fidThres)

np.savetxt(full_path_Omega, Omega)
np.savetxt(full_path_t_fidThres, t_fidThres)




stop = timeit.default_timer()
print('Time: ', stop - start)
