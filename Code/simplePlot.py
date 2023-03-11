# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 10:20:53 2023

@author: sony
"""

import matplotlib.pyplot as plt

alpha = [0.0001,0.001,0.01,0.1,1,10]
omega = [2,2,2,4,10,43]


plt.title('Optimal Energy Splitting for Different Coupling Strengths')
plt.xlabel('Alpha')
plt.ylabel(r"Omega $(ps^{-1})$")
plt.plot(alpha,omega)
plt.show()