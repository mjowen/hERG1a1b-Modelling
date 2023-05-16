#hERG 1a/1b independence model - Fitting to data and test cases
#%%
# Imports
import matplotlib.pyplot as plt
import myokit
import numpy as np
import scipy

#Load model
model = myokit.load_model('.\Models\hERG1a1bIndependence.mmt')
