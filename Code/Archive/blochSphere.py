from qutip import *
import matplotlib.pyplot as plt

import numpy as np

b = Bloch()

# Define the initial qubit state as |1⟩
psi = Qobj([[1/2], [1]])

psi = psi / np.sqrt((psi.dag() * psi)[0, 0])  # Normalize the state vector

# Add the state vector to the Bloch sphere
b.add_states(psi)

# Add labels for the standard qubit states
labels = ['|0⟩', '|1⟩']
b.vector_labels = [labels, ['-'.center(len(l), '-') for l in labels]]

# Set the view and plot the Bloch sphere
b.show()
#b.view = [30,30]
#b.show(figsize=[8,8])
# Save the plot as a PDF file
plt.savefig("bloch_sphere.pdf")