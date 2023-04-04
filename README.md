# qubit_reset
Optimising Qubit Reset in Open Quantum Systems

This README file describes the general file structure and also serves as a bit of documentation.


# Contents
## File structure
- [**data**]
	empty
- [**Keras_NN**]
	- **_pycache**
	- `architecture.py` 
 	- `hypOp.py`  
	- `TrainModel.py` 
- [**models**]
	- `FINALmodelRealMatched_30.pth`
- [**Plots**]
	- `CSV4histOfTrackCount.py` 
 	- `CSV4histOfVariables.py`  
	- `PlotTrackCount.py`
	- `PlotVariables.py` 
- [**Preprocessing**]

- [**PyTorch_NN**]
	- `hyperparam_search.py` 
 	- `TrainNN.py` 

- `README.txt` 
- `XGboostClassf.py` 


## Description of contents

### data
Where data from preprocessing is output for use in training.

### Keras_NN
Contains Python files for training, testing and optimising a Keras NN.
- `architecture.py`: Defines architecture for use in hyperparameter search.
- `hypOp.py`: Random Search for Keras Hyperparameter Optimisation. 
- `TrainModel.py`: Train and Test Keras NN on Dataset.

### models
Where trained PyTorch models are output.
- `FINALmodelRealMatched_30.pth`: Final trained model with optimal hyperparameters.

### Plots
- `CSV4histOfTrackCount.py`: Produce a CSVs to use for plotting the number of tracks per PV.
- `CSV4histOfVariables.py` : Produce CSVs for plotting input variable histograms
- `PlotTrackCount.py`: Plot the number of tracks per PV.
	Takes the two CSVs created in CSV4histOfTrackCount.py and plots desired histogram.
	This python file was ran on my windows laptop using a spyder environment to view the plots.
- `PlotVariables.py`: Plot the  input variables per PV.
	Takes the two CSVs created in CSV4histOfVariables.py and plots desired histogram.
	This python file was also ran on my windows laptop using a spyder environment to view the plots.

### Preprocessing
- `preprocessing.py`: Preprocessing converting ROOT file to pandas df and outputting as CSV.

### PyTorch_NN
- `hyperparam_search.py`: Tune PyTorch hyperparameters with Ray Tune.
- `predictions.py`: Make predictions on an alternative dataset using a trained model.
- `TrainNN.py`: Train and test PyTorch Neural Network for Binary Classification.


- `XGboostClassf.py`: Simple Binary Classifier using XGboost.



## Tutorial
### Brief tutorial of final preprocessing and NN used.

Dataset CSVs are produced using /Preprocessing/prepocessing.py (need to specify where the root file is that you want to convert).
The NN is then built, trained and tested in /PyTorch_NN/TrainNN.py where a trained model is also output (this model can then be used
to make predictions on a new dataset).
