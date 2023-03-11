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


# Pick System parameters
Omega = [1,7,8,9]
omega_cutoff = 2
alpha = 1
# Temperature in kelvin
T = 0 # K
# Temperature in units of 1/(ps kB)
# 1 K = 0.1309 1/(ps kB)
temperature = T*0.1309 # 1/(ps kB)

# initial state = up_density_matrix or mixed_density_matrix
init_st = up_density_matrix

# TEMPO parameters
dt=0.01
dkmax=200
epsrel=10**(-6)
# duration (ps)
dur = 1.4

# Time steps for Master Equation
dtME = 0.01
num_stepsME = dur/dtME


tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
                                          epsrel=epsrel)

# Spectral density ohmic
J0 = 4*alpha*Omega[0]*np.exp(-2*Omega[0]/omega_cutoff)
J1 = 4*alpha*Omega[1]*np.exp(-2*Omega[1]/omega_cutoff)
J2 = 4*alpha*Omega[2]*np.exp(-2*Omega[2]/omega_cutoff)
J3 = 4*alpha*Omega[3]*np.exp(-2*Omega[3]/omega_cutoff)

# B-E Occupancy for master equation
def N(T):
    if T == 0:
        N = 0
    else:
        N = (1/(math.exp(((1.055e-34)*(10**12))/((1.381e-23)*T))-1))
    return N

# # Make plot of spectral density and time dependent splitting
plt.title(r'Spectral Density as a Function of Energy Splitting, $\alpha={0}$'.format(alpha))
w = np.linspace(0,30,300)#20,200
J = 4*alpha*w*np.exp(-2*w**2/omega_cutoff**2)
plt.plot(w, J)#.format(coeff) tau = 0.245

plt.axvline(x = Omega[0],color='r',linestyle='dashed',
            label=r"$\Omega_{0}={1}$".format(0,Omega[0]))
plt.axvline(x = Omega[1],color='b',linestyle='dashed',
            label=r"$\Omega_{0}={1}$".format(1,Omega[1]))
plt.axvline(x = Omega[2],color='g',linestyle='dashed',
            label=r"$\Omega_{0}={1}$".format(2,Omega[2]))
plt.axvline(x = Omega[3],color='m',linestyle='dashed',
            label=r"$\Omega_{0}={1}$".format(3,Omega[3]))

plt.xlabel(r'$\Omega (ps^{-1})$')
plt.ylabel(r"$J(\Omega) \mathrm{ps}^{-1}$")
plt.legend()
#plt.ylim([0, Omega[0]+0.1])
# Save plots
desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'Spect&Om_constOm='+str(Omega)+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()






# Specify correlations
correlations = oqupy.PowerLawSD(alpha=alpha,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='gaussian',
                                temperature=(temperature))


bath = oqupy.Bath(sigma_x, correlations)








# #############
# process_tensor = oqupy.pt_tempo_compute(bath=bath,
#                                         start_time=0.0,
#                                         end_time=dur,
#                                         parameters=tempo_parameters)

# ###### Save PT ######## Doesn't work if already saved
# desired_folder_PT = 'C:/Users/sony/Desktop/Project/process_tensor/'

# file_name_PT = 'alpha='+str(alpha)+'w_c='+str(omega_cutoff)+'_T='+str(T)+'_dt='+str(dt)+\
#     '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_PT.txt'

# full_path_PT = os.path.join(desired_folder_PT, file_name_PT)

# oqupy.SimpleProcessTensor.export(self=process_tensor,
#                                   filename=full_path_PT)
# # ######################################

# Read in PT
full_path_PT = 'C:/Users/sony/Desktop/Project/process_tensor/alpha=1w_c=0.5_T=0_dt=0.01_dkmax=200_epsrel=1e-06_dur=1.4_PT.txt'

process_tensor = oqupy.process_tensor.FileProcessTensor(filename=full_path_PT,
                                                        mode = 'read')

# ################################################
#################################################

# Specify system
system = oqupy.System(Omega[0]*sigma_z)

# Compute
dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
                                system=system,
                                initial_state=init_st,
                                start_time=0.0 # needed?, progress_type
                                ,num_steps=dur/dt)

t0, s_z0 = dynamics.expectations(sigma_z, real=True)

# Save
desired_folder_data = 'C:/Users/sony/Desktop/Project/PlotData'

file_name_t = 'PltDat_constOm='+str(Omega)+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'
    
file_name_s_z = 'PltDat_constOm='+str(Omega)+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'
    
full_path_t = os.path.join(desired_folder_data, file_name_t)
full_path_s_z = os.path.join(desired_folder_data, file_name_s_z)

np.savetxt(full_path_t, t0)
np.savetxt(full_path_s_z, s_z0)
#np.loadtxt
###########



# Derived master equation
#gamma1 = gamma*(N+1)
gamma_1 = 2*J0*np.pi*(N(T)+1)#*(N+1)
# gamma2 = gamma*(N)
gamma_2 = 2*J0*np.pi*N(T) #*N

Gamma = [gamma_1, gamma_2]

# This is in the same form as in notes
lindblad_operators  = [sigma_minus, sigma_plus]
    

systemD = oqupy.System(
          Omega[0] * sigma_z,
        gammas=Gamma,
        lindblad_operators=lindblad_operators)

# Compute
dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                  start_time = 0,dt=dtME,num_steps=num_stepsME)

tME0, s_zME0 = dynamics.expectations(sigma_z, real=True)

#####################################################################
#####################################################################
################################################
# Repeat for different system

# Specify system
system = oqupy.System(Omega[1]*sigma_z)

dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
                                system=system,
                                initial_state=init_st,
                                start_time=0.0 # needed?, progress_type
                                ,num_steps=dur/dt)
#and plot the result
t1, s_z1 = dynamics.expectations(sigma_z, real=True)

# Save
desired_folder_data = 'C:/Users/sony/Desktop/Project/PlotData'

file_name_t = 'PltDat_constOm='+str(Omega[1])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'
    
file_name_s_z = 'PltDat_constOm='+str(Omega[1])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'
    
full_path_t = os.path.join(desired_folder_data, file_name_t)
full_path_s_z = os.path.join(desired_folder_data, file_name_s_z)

np.savetxt(full_path_t, t1)
np.savetxt(full_path_s_z, s_z1)
#np.loadtxt
###########

# Derived master equation

#gamma1 = gamma*(N+1)
gamma_1 = 2*J1*np.pi*(N(T)+1)#*(N+1)
# gamma2 = gamma*(N)
gamma_2 = 2*J1*np.pi*N(T) #*N

Gamma = [gamma_1, gamma_2]

# This is in the same form as in notes
lindblad_operators  = [sigma_minus, sigma_plus]
    

systemD = oqupy.System(
          Omega[1] * sigma_z,
        gammas=Gamma,
        lindblad_operators=lindblad_operators)

# up_density_matrix
dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                  start_time = 0,dt=dtME,num_steps=num_stepsME)

tME1, s_zME1 = dynamics.expectations(sigma_z, real=True)

######################################################################
######################################################################
###################################

# Repeat for different system

# Specify system
system = oqupy.System(Omega[2]*sigma_z)

dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
                                system=system,
                                initial_state=init_st,
                                start_time=0.0 # needed?, progress_type
                                ,num_steps=dur/dt)
#and plot the result
t2, s_z2 = dynamics.expectations(sigma_z, real=True)

# Save
desired_folder_data = 'C:/Users/sony/Desktop/Project/PlotData'

file_name_t = 'PltDat_constOm='+str(Omega[2])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'
    
file_name_s_z = 'PltDat_constOm='+str(Omega[2])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'
    
full_path_t = os.path.join(desired_folder_data, file_name_t)
full_path_s_z = os.path.join(desired_folder_data, file_name_s_z)

np.savetxt(full_path_t, t2)
np.savetxt(full_path_s_z, s_z2)
#np.loadtxt
###########

# Derived master equation

#gamma1 = gamma*(N+1)
gamma_1 = 2*J2*np.pi*(N(T)+1)#*(N+1)
# gamma2 = gamma*(N)
gamma_2 = 2*J2*np.pi*N(T) #*N

Gamma = [gamma_1, gamma_2]

# This is in the same form as in notes
lindblad_operators  = [sigma_minus, sigma_plus]
    

systemD = oqupy.System(
          Omega[2] * sigma_z,
        gammas=Gamma,
        lindblad_operators=lindblad_operators)

# up_density_matrix
dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                  start_time = 0,dt=dtME,num_steps=num_stepsME)

tME2, s_zME2 = dynamics.expectations(sigma_z, real=True)

######################################################################
######################################################################
# ###################################

# Repeat for different system

# Specify system
system = oqupy.System(Omega[3]*sigma_z)

dynamics = oqupy.compute_dynamics(process_tensor=process_tensor,
                                system=system,
                                initial_state=init_st,
                                start_time=0.0 # needed?, progress_type
                                ,num_steps=dur/dt)
#and plot the result
t3, s_z3 = dynamics.expectations(sigma_z, real=True)

# Save
desired_folder_data = 'C:/Users/sony/Desktop/Project/PlotData'

file_name_t = 'PltDat_constOm='+str(Omega[3])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'
    
file_name_s_z = 'PltDat_constOm='+str(Omega[3])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'
    
full_path_t = os.path.join(desired_folder_data, file_name_t)
full_path_s_z = os.path.join(desired_folder_data, file_name_s_z)

np.savetxt(full_path_t, t2)
np.savetxt(full_path_s_z, s_z2)
#np.loadtxt
# ###########

# Derived master equation

#gamma1 = gamma*(N+1)
gamma_1 = 2*J3*np.pi*(N(T)+1)#*(N+1)
# gamma2 = gamma*(N)
gamma_2 = 2*J3*np.pi*N(T) #*N

Gamma = [gamma_1, gamma_2]

# This is in the same form as in notes
lindblad_operators  = [sigma_minus, sigma_plus]
    

systemD = oqupy.System(
          Omega[3] * sigma_z,
        gammas=Gamma,
        lindblad_operators=lindblad_operators)


dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                  start_time = 0,dt=dtME,num_steps=num_stepsME)

tME3, s_zME3 = dynamics.expectations(sigma_z, real=True)



#######################################################

# Plot TEMPO Solution
plt.plot(t0, s_z0,color='r',label=r'TEMPO $\Omega_{0}={1}$ '.format(0,Omega[0]))
plt.plot(t1, s_z1,color='b',label=r'TEMPO $\Omega_{0}={1}$ '.format(1,Omega[1]))
plt.plot(t2, s_z2,color='g',label=r'TEMPO $\Omega_{0}={1}$ '.format(2,Omega[2]))
plt.plot(t3, s_z3,color='m',label=r'TEMPO $\Omega_{0}={1}$ '.format(3,Omega[3]))

plt.title(r'Decay from Up Density Matrix State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
plt.grid()

desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'constOmsTEMPO='+str(Omega)+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()



# Plot TEMPO solution zoomed in
plt.plot(t0, s_z0,color='r',label=r'TEMPO $\Omega_{0}={1}$ '.format(0,Omega[0]))
plt.plot(t1, s_z1,color='b',label=r'TEMPO $\Omega_{0}={1}$ '.format(1,Omega[1]))
plt.plot(t2, s_z2,color='g',label=r'TEMPO $\Omega_{0}={1}$ '.format(2,Omega[2]))
plt.plot(t3, s_z3,color='m',label=r'TEMPO $\Omega_{0}={1}$ '.format(3,Omega[3]))

plt.title(r'Decay from Up Density Matrix State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
plt.grid()
plt.ylim([0.9,1.1])
plt.xlim([dur/2, dur])

desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'constOmsTEMPOzoom='+str(Omega)+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()



# Plot master equation solution
plt.plot(tME0, s_zME0,color='r', label=r'QME $\Omega_{0}={1}$ '.format(0,Omega[0]))
plt.plot(tME1, s_zME1,color='b', label=r'QME $\Omega_{0}={1}$ '.format(1,Omega[1]))
plt.plot(tME2, s_zME2,color='g', label=r'QME $\Omega_{0}={1}$ '.format(2,Omega[2]))
plt.plot(tME3, s_zME3,color='m', label=r'QME $\Omega_{0}={1}$ '.format(3,Omega[3]))


plt.title(r'Decay from Up Density Matrix State at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()
plt.grid()

desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'constOmsQME='+str(Omega)+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()


stop = timeit.default_timer()
print('Time: ', stop - start)