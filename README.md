# Optimising Qubit Reset in Open Quantum Systems


This README file describes the general file structure and also serves as a bit of documentation.


# Contents
## File structure
- [**Code**]
	- **Archive**
	- `difOmegas.py` 
 	- `OptimalEnegySpacing_final_plot.py`  
	- `plot_difOmegas.py`
 	- `plot_MEvsPT-TEMPO.py`  
	- `timeDepH_tanh_loop.py`
- [**FOMdata**]
	
- [**PlotData**]

- [**Plots**]

- [**process_tensor**]

- `README.txt` 



## Description of contents


### Code
Contains Python files for solving system dynamics with OQuPy and producing plots.
- `difOmegas.py` : Main file for solving dynamics with the PT-TEMPO method and master equation solver.
- `OptimalEnegySpacing_final_plot.py` : Produces the final plot for use in report (Figure 11b), Optimal Energy Spacing for Varying Coupling Strength. Uses data produced from difOmegas.py. 
- `plot_difOmegas.py` : Produces the plot for use in report, Figure 7. Uses data produced from difOmegas.py to plot spectral density, decay from excited state and fidelity evolution.
- `plot_MEvsPT-TEMPO.py` : Produces the plot for use in report, Figure 8. Uses data produced from difOmegas.py to plot a comparison of dynamics for the PT-TEMPO method and master equation.
- `timeDepH_tanh_loop.py` : For solving dynamics with time-dependent Hamiltonians. 



### FOMdata
Where select data produced in difOmegas.py can be stored and used for the plots in OptimalEnegySpacing_final_plot.py.


### PlotData


### Plots


### process_tensor




## Tutorial
### Brief tutorial of .
