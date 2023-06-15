import myokit
import numpy as np
import matplotlib.pyplot as plt
from Code import myokitPlot

#Plots of M2 for validation
myokitPlot.multiProtocolPlot('perissinotti 2018 m2','Perissinotti 2018\steadyStateActivation', factorVars = ['states.hERG1b','environment.T'], factorVals = [0,296])
myokitPlot.multiProtocolPlot('perissinotti 2018 m2','Perissinotti 2018\steadyStateActivation', factorVars = ['states.hERG1b','environment.T'], factorVals = [1,296])

#%% Plot output of VGC-KiMo for M2 validation
data = vgckimossaM21bcsv #Import one of vgckimo-ssaM2.csv, vgckimo-ssaM2-1a.csv, or vgckimo-ssaM2-1b.csv
time = np.asarray(data[1:,0], dtype=float)
currents = np.asarray(data[1:,1:], dtype=float) #60:-1:-90 (0-15)
plt.figure()
for i in range(16):
    plt.plot(time,currents[:,i])
plt.show()

#%% State occupancy for validatating against supp S9 and S11
use1a = True
modelName = 'perissinotti 2018 m1'
if use1a:
    protocol = 'Perissinotti 2018\AP M1'
    factorVars = ['states.hERG1b']
    factorVals = [0]
    initialValues = [0.947,0.001,0.021,0.030,0.001]
else:
    protocol = 'Perissinotti 2018\AP M1 1b'
    factorVars = ['states.hERG1b']
    factorVals = [1]
    initialValues = [1,0,0,0,0]

myokitPlot.simplePlot(modelName = modelName, protocol=protocol, factorVars = factorVars, factorVals = factorVals, plotProtocol = True, initialValues = initialValues)
states = ['states.C3', 'states.C2', 'states.C1', 'states.O', 'states.I']
for i in range(len(states)):
    myokitPlot.simplePlot(modelName = modelName, protocol=protocol, factorVars = [], factorVals = [], plotProtocol = False, ylabel = states[i], output = states[i], initialValues = initialValues)
