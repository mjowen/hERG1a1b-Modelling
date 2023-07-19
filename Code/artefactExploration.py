import os
import numpy as np
import matplotlib.pyplot as plt
import csv
if os.getcwd()[-9:] == 'Modelling':
    os.chdir('Code')
import myokitPlot
if os.getcwd()[-4:] == 'Code':
    os.chdir('..')
modelName = 'fink2008-artefactFullModel'
protocol = 'Additional Protocols/staircase-ramp'
factorVars = [['voltageclamp.rseries']]*5
factorVals = [[0.1*5e-3], [5e-3], [2*5e-3], [5*5e-3], [10*5e-3]]

myokitPlot.multiSettingPlot(modelName, protocol, protocolName = 'voltage', factorVars = factorVars, factorVals = factorVals, plotProtocol = True, output = 'membrane.V', title = 'Rseries')
#Higher series resistances result in voltages resisting being moved away from EK. This can be enough to change flat voltage steps into slopes or exponential decays. Increasing compension reduces this effect.

factorVars = [['membrane.cm']]*3
factorVals = [[0.001*15], [15], [5000*15]]
myokitPlot.multiSettingPlot(modelName, protocol, protocolName = 'voltage', factorVars = factorVars, factorVals = factorVals, plotProtocol = False, output = 'membrane.V', title = 'Membrane capacitance')
#For very large cell membrane capacitance (~15000), then all voltages look like exponential decays to their intended values

factorVars = [['voltageclamp.alpha', 'voltageclamp.rseries', 'membrane.cm']]*3
factorVals = [[0.001, 5e-2, 3000*15], [0.5, 5e-2, 3000*15], [0.999, 5e-2, 3000*15]]
myokitPlot.multiSettingPlot(modelName, protocol, protocolName = 'voltage', factorVars = factorVars, factorVals = factorVals, plotProtocol = False, output = 'membrane.V', title = 'Supercharging')
#I cant find nice settings which show the supercharging making a difference
