# -*- coding: utf-8 -*-
"""

Generate Plots from saved plot data

"""

import sys
import os
sys.path.insert(0,'..')

import oqupy
import numpy as np
import matplotlib.pyplot as plt
import timeit
import math
from scipy.fft import fft, fftshift

start = timeit.default_timer()


# Pick parameters you want to plot
# Pick System parameters
Omega = [2.5,8,8.5,9]
omega_cutoff = 5.0
alpha = 0.1
# Temperature in kelvin
T = 0 # K
# Temperature in units of 1/(ps kB)
# 1 K = 0.1309 1/(ps kB)
temperature = T*0.1309 # 1/(ps kB)



# TEMPO parameters
dt=0.05
dkmax=200
epsrel=10**(-6)
# duration (ps)
dur = 18


# Load data
folder_data = 'C:/Users/sony/Desktop/Project/PlotData'


#J and w

# Load data
desired_folder_data = 'C:/Users/sony/Desktop/Project/PlotData'

file_name_J = 'J='+str(Omega[0])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.txt'

file_name_w = 'w='+str(Omega[0])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.txt'    

    
full_path_J = os.path.join(desired_folder_data, file_name_J)
full_path_w = os.path.join(desired_folder_data, file_name_w)

J = np.loadtxt(full_path_J)
w = np.loadtxt(full_path_w)

#Plot

# # Make plot of spectral density and time dependent splitting
plt.title(r'Spectral Density as a Function of Energy Splitting')
# w = np.linspace(0,20,200)#20,200
# J = 4*alpha*w*np.exp(-2*w/omega_cutoff)
plt.plot(w, J)#.format(coeff) tau = 0.245

plt.axvline(x = Omega[0],color='r',linestyle='dashed',
            label=r"$\Omega_{0}={1}$".format(0,Omega[0]))
plt.axvline(x = Omega[1],color='b',linestyle='dashed',
            label=r"$\Omega_{0}={1}$".format(1,Omega[1]))
plt.axvline(x = Omega[2],color='g',linestyle='dashed',
            label=r"$\Omega_{0}={1}$".format(2,Omega[2]))
plt.axvline(x = Omega[3],color='m',linestyle='dashed',
            label=r"$\Omega_{0}={1}$".format(3,Omega[3]))

plt.xlabel(r'$\Omega \; (ps^{-1})$')
plt.ylabel(r"$J(\Omega) \; (\mathrm{ps}^{-1})$")
plt.legend()
#plt.ylim([0, Omega[0]+0.1])
# Save plots
desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'Spect&Om_constOm='+str(Omega)+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)
plt.show()



# 1st
# !!! change to Omega[0]
file_name_t = 'PltDat_constOm='+str(Omega[0])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'
    
file_name_s_z = 'PltDat_constOm='+str(Omega[0])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'
    

full_path_t = os.path.join(folder_data, file_name_t)
full_path_s_z = os.path.join(folder_data, file_name_s_z)

t0 = np.loadtxt(full_path_t)
s_z0 = np.loadtxt(full_path_s_z)

# MASTER EQUATION
file_name_tME = 'PltDat_constOm='+str(Omega[0])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_TimeME.txt'
    
file_name_s_zME = 'PltDat_constOm='+str(Omega[0])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_zME.txt'
    

full_path_tME = os.path.join(folder_data, file_name_tME)
full_path_s_zME = os.path.join(folder_data, file_name_s_zME)

tME0 = np.loadtxt(full_path_t)
s_zME0 = np.loadtxt(full_path_s_z)



# 2nd
file_name_t = 'PltDat_constOm='+str(Omega[1])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'
    
file_name_s_z = 'PltDat_constOm='+str(Omega[1])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'
    

full_path_t = os.path.join(folder_data, file_name_t)
full_path_s_z = os.path.join(folder_data, file_name_s_z)

t1 = np.loadtxt(full_path_t)
s_z1 = np.loadtxt(full_path_s_z)

# MASTER EQUATION
file_name_tME = 'PltDat_constOm='+str(Omega[1])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_TimeME.txt'
    
file_name_s_zME = 'PltDat_constOm='+str(Omega[1])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_zME.txt'
    

full_path_tME = os.path.join(folder_data, file_name_tME)
full_path_s_zME = os.path.join(folder_data, file_name_s_zME)

tME1 = np.loadtxt(full_path_t)
s_zME1 = np.loadtxt(full_path_s_z)




# 3rd
file_name_t = 'PltDat_constOm='+str(Omega[2])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'
    
file_name_s_z = 'PltDat_constOm='+str(Omega[2])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'
    

full_path_t = os.path.join(folder_data, file_name_t)
full_path_s_z = os.path.join(folder_data, file_name_s_z)

t2 = np.loadtxt(full_path_t)
s_z2 = np.loadtxt(full_path_s_z)

# MASTER EQUATION
file_name_tME = 'PltDat_constOm='+str(Omega[2])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_TimeME.txt'
    
file_name_s_zME = 'PltDat_constOm='+str(Omega[2])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_zME.txt'
    

full_path_tME = os.path.join(folder_data, file_name_tME)
full_path_s_zME = os.path.join(folder_data, file_name_s_zME)

tME2 = np.loadtxt(full_path_t)
s_zME2 = np.loadtxt(full_path_s_z)




# 4th
file_name_t = 'PltDat_constOm='+str(Omega[3])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_Time.txt'
    
file_name_s_z = 'PltDat_constOm='+str(Omega[3])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_z.txt'
    

full_path_t = os.path.join(folder_data, file_name_t)
full_path_s_z = os.path.join(folder_data, file_name_s_z)

t3 = np.loadtxt(full_path_t)
s_z3 = np.loadtxt(full_path_s_z)


# MASTER EQUATION
file_name_tME = 'PltDat_constOm='+str(Omega[3])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_TimeME.txt'
    
file_name_s_zME = 'PltDat_constOm='+str(Omega[3])+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'_S_zME.txt'
    

full_path_tME = os.path.join(folder_data, file_name_tME)
full_path_s_zME = os.path.join(folder_data, file_name_s_zME)

tME3 = np.loadtxt(full_path_t)
s_zME3 = np.loadtxt(full_path_s_z)




#######################################################

# PLOT DATA

#######################################################

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
plt.minorticks_on()
plt.grid(which='both')
plt.ylim([-1,-0.8])
plt.xlim([dur/2, dur])

desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'constOmsTEMPOzoom='+str(Omega)+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()




# Fourier transform


s_z0ft = fftshift(fft(s_z0))
s_z1ft = fftshift(fft(s_z1))
s_z2ft = fftshift(fft(s_z2))
s_z3ft = fftshift(fft(s_z3))
# s_z0ft = np.fft.fft(s_z0)
# s_z1ft = np.fft.fft(s_z1)
# s_z2ft = np.fft.fft(s_z2)
# s_z3ft = np.fft.fft(s_z3)

# Plot
plt.plot(t0, s_z0ft,color='r',label=r'TEMPO $\Omega_{0}={1}$ '.format(0,Omega[0]))
plt.plot(t1, s_z1ft,color='b',label=r'TEMPO $\Omega_{0}={1}$ '.format(1,Omega[1]))
plt.plot(t2, s_z2ft,color='g',label=r'TEMPO $\Omega_{0}={1}$ '.format(2,Omega[2]))
plt.plot(t3, s_z3ft,color='m',label=r'TEMPO $\Omega_{0}={1}$ '.format(3,Omega[3]))

plt.title(r'Fourier Transform at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel(r'FT of $\langle\sigma_z\rangle$')
plt.legend()
plt.grid()


desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'constOmsTEMPOft='+str(Omega)+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()


# Plot FOMs
plt.plot(t0, (((s_z0-1)/(-2))*100),color='r',label=r'TEMPO $\Omega_{0}={1}$ '.format(0,Omega[0]))
plt.plot(t1, (((s_z1-1)/(-2))*100),color='b',label=r'TEMPO $\Omega_{0}={1}$ '.format(1,Omega[1]))
plt.plot(t2, (((s_z2-1)/(-2))*100),color='g',label=r'TEMPO $\Omega_{0}={1}$ '.format(2,Omega[2]))
plt.plot(t3, (((s_z3-1)/(-2))*100),color='m',label=r'TEMPO $\Omega_{0}={1}$ '.format(3,Omega[3]))
plt.axhline(y = 99,color='c',linestyle='dashed',
            label="99%")


plt.title(r'Time Taken to Reach Fidelity at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel('Fidelity (%)')
plt.legend()
plt.minorticks_on()
plt.grid(which='both')

desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'FOMs='+str(Omega)+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
    '_dkmax='+str(dkmax)+'_epsrel='+str(epsrel)+'_dur='+str(dur)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)

plt.show()


# Plot FOMs zoomed in
plt.plot(t0, (((s_z0-1)/(-2))*100),color='r',label=r'TEMPO $\Omega_{0}={1}$ '.format(0,Omega[0]))
plt.plot(t1, (((s_z1-1)/(-2))*100),color='b',label=r'TEMPO $\Omega_{0}={1}$ '.format(1,Omega[1]))
plt.plot(t2, (((s_z2-1)/(-2))*100),color='g',label=r'TEMPO $\Omega_{0}={1}$ '.format(2,Omega[2]))
plt.plot(t3, (((s_z3-1)/(-2))*100),color='m',label=r'TEMPO $\Omega_{0}={1}$ '.format(3,Omega[3]))
plt.axhline(y = 99,color='c',linestyle='dashed',
            label="99%")


plt.title(r'Time Taken to Reach Fidelity at $\alpha$={0}, T = {1} K'.format(alpha,T)) 
plt.xlabel('Time (ps)')
plt.ylabel('Fidelity (%)')
plt.legend()
plt.minorticks_on()
plt.grid(which='both')
plt.ylim([95.5,99.2])
plt.xlim([dur/2, dur])

desired_folder = 'C:/Users/sony/Desktop/Project/Plots'
file_name = 'FOMsZm='+str(Omega)+'a='+str(alpha)+'_T='+str(T)+'_dt='+str(dt)+\
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