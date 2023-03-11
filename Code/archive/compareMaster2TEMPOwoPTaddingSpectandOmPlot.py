"""
System plot giving TEMPO Hamiltonians and plot with derived master equation
 in linblad form
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
Omega = [1,1.5,0.5]
omega_cutoff = 5.0
alpha = 0.11
# Temperature in kelvin
T = 0 # K
# Temperature in units of 1/(ps kB)
# 1 K = 0.1309 1/(ps kB)
temperature = T*0.1309 # 1/(ps kB)

# initial state = up_density_matrix or mixed_density_matrix
init_st = up_density_matrix


# Specify system
system = oqupy.System(Omega*sigma_z)

# Specify correlations
correlations = oqupy.PowerLawSD(alpha=alpha,
                                zeta=1,
                                cutoff=omega_cutoff,
                                cutoff_type='exponential',
                                temperature=(temperature))

bath = oqupy.Bath(sigma_x, correlations)

# TEMPO parameters
dt=0.05
dkmax=70
epsrel=10**(-6)
# duration
dur = 5.0

# For Master Equation (need to change if changing dtME)
num_steps = dur*100

tempo_parameters = oqupy.TempoParameters(dt=dt , dkmax=dkmax,\
                                          epsrel=epsrel)


dynamics = oqupy.tempo_compute(system=system,
                                bath=bath,
                                initial_state=init_st,
                                start_time=0.0,
                                end_time=dur,
                                parameters=tempo_parameters)

#and plot the result
t, s_z = dynamics.expectations(sigma_z, real=True)

# Save
desired_folder_data = 'C:/Users/sony/Desktop/Project/PlotData'

file_name_t = 'alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'
    
file_name_s_z = 'alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'
    
full_path_t = os.path.join(desired_folder_data, file_name_t)
full_path_s_z = os.path.join(desired_folder_data, file_name_s_z)

np.savetxt(full_path_t, t)
np.savetxt(full_path_s_z, s_z)

#np.loadtxt

##############################################################################
# Derived master equation

# B-E Occupancy
def N(T):
    if T == 0:
        N = 0
    else:
        N = (1/(math.exp(((1.055e-34)*(10**12))/((1.381e-23)*T))-1))
    return N

# Checks
print('T is', T)        
print('N is', N(T))

# Spectral density ohmic
J = 4*alpha*Omega*np.exp(-2*Omega/omega_cutoff)
# # Spectral density super-ohmic
# J = 2*alpha*(Omega**(3)/omega_cutoff**(2))*np.exp(-Omega/omega_cutoff)
#gamma1 = gamma*(N+1)
gamma_1 = 2*J*np.pi*(N(T)+1)#*(N+1)
# gamma2 = gamma*(N)
gamma_2 = 2*J*np.pi*N(T) #*N

Gamma = [gamma_1, gamma_2]

# This is in the same form as in notes
lindblad_operators  = [sigma_minus, sigma_plus]



# Make plot of spectral density and time dependent splitting
print("J =",J)
t_plot = np.linspace(0,20,int(dur)*100) #
plt.title('Energy Splitting and Spectral Density')
plt.axhline(y=J,label=r"Sup $J$")
plt.axhline(y=Omega,label=r"$\Omegat") 
plt.xlabel("Time (ps)")
plt.ylabel(r"$\omega(t) \mathrm{ps}^{-1}$")
plt.legend()
plt.show()
    

systemD = oqupy.System(
         Omega * sigma_z,
        gammas=Gamma,
        lindblad_operators=lindblad_operators)

# up_density_matrix
dynamics = oqupy.compute_dynamics(system = systemD, initial_state = init_st,\
                                  start_time = 0,dt=0.01,num_steps=num_steps)

t2, s_z2 = dynamics.expectations(sigma_z, real=True)


# Plot
plt.plot(t, s_z, label=r'TEMPO $\alpha={0}$ '.format(alpha))
plt.plot(t2, s_z2, label=r'QME $\gamma = {0} \pi$ '.format(J.round(3)))

plt.title('Decay from Up Density Matrix State at T = {0} K'.format(T)) 
plt.xlabel(r'$t\,\Omega$')
plt.ylabel(r'$\langle\sigma_z\rangle$')
plt.legend()

desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'alpha='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()


#fig1, fig2 = plt.subplots() 



stop = timeit.default_timer()
print('Time: ', stop - start)