"""

Generates Final Optimal Energy Spacing Plots

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
from brokenaxes import brokenaxes
import numpy as np
import matplotlib as mpl
from matplotlib.pyplot import cm

start = timeit.default_timer()


# Choose threshold fidelity
fid_thres = 99

# Specify which coupling strengths to plot
alpha = [0.001,0.01,0.05,0.1]



Omega=[]
t_fidThres=[]
for i in range(len(alpha)):
    Omega.append('Omega'+str(i))
    t_fidThres.append('t_fidThres'+str(i))

# Load data
folder_data = '/Project/PlotData'


# Load data for each coupling strength
for i in range(len(alpha)):
    
    # Load data
    desired_folder_data = '/FOMdata/select'
    
    file_name_Omega = 'fid_thres'+str(fid_thres)+'a='+str(alpha[i])+'_Omega.txt'
        
    file_name_t_fidThres = 'fid_thres'+str(fid_thres)+'a='+str(alpha[i])+'_t_fidThres.txt'
        
    full_path_Omega = os.path.join(desired_folder_data, file_name_Omega)
    full_path_t_fidThres = os.path.join(desired_folder_data, file_name_t_fidThres)
    
    Omega[i] = np.loadtxt(full_path_Omega)
    t_fidThres[i] = np.loadtxt(full_path_t_fidThres)



#######################################################

# PLOT DATA

#######################################################

# Plot 
fig = plt.figure()
# Use broken axis
bax = brokenaxes( ylims=((8, 30), (220, 300)), yscale='log')#,  xlims=((0, 5), (8, 10)),

# Set colour wheel for parabolas
# n = int(len(Omega))
# c = cm.viridis(np.linspace(0, 1, n))

# Manually set colours for parabolas
c = ['limegreen','darkcyan','navy','purple']

# Plot data for each alpha
for i in range(len(alpha)):
    bax.plot(Omega[i], t_fidThres[i],color=c[i],label=r'$\alpha={0}$'.format(alpha[i]))

# Plot horizontal lines to locate minima
# bax.axvline(x = 5.2,color='red',linestyle='dashed')
# bax.axvline(x = 8.6,color='red',linestyle='dashed')


bax.first_col[0].set_yticks([225,250,300])
bax.first_col[1].set_yticks([8,10,15,20,25,30])

# Trying to add an axis on the right as well
# bax.secondary_yaxis()
# bax.secondary_yaxis().set_yticks([],[])
# bax.secondary_xaxis().set_xticks([],[])
#bax.tick_params(top=False, labeltop=False, right=False, labelright=False)


bax.legend(loc=1)
bax.set_title('Optimal Energy Spacing')
bax.set_xlabel(r'Energy Splitting Frequency $\Omega, \; (ps^{-1})$')
bax.set_ylabel('Time to reach Threshold Fidelity (ps)')

desired_folder = '/Plots/FOMthres'
file_name = 'finalFOM_alpha='+str(alpha)+'.pdf'
full_path = os.path.join(desired_folder, file_name)
plt.savefig(full_path)


fig.show()


stop = timeit.default_timer()
print('Time: ', stop - start)