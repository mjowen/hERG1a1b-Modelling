#Current trace plots
#%% Imports and function definitions - Only need to run once
import matplotlib.pyplot as plt
import myokit
import numpy as np

def currentTracePlots(modelNames, protocol):
    """
    Simulates the models described by 'modelNames' under its possible specifications (hand-coded per model). Designed for protocol files with a single potential function. It will plot model current for model in modelNames on a new graph which explores all their possible specifications and then plot the potential on a new graph.
    Parameters
    ----------
    modelNames : list of strings
        Name of myokit model files in the Models folder. For example, 'larsen2010' for the Models/larsen2010.mmt model. 
    protocol : string
        Filename of a protocol to load from the 'Protocols' folder. For example, 'Clerx 2019\Traditional\Pr2' loads protocol 2 from the 4 ways to fit an ion channel paper (full filename: .\Protocols\Clerx 2019\Traditional\Pr2.csv). Default is the staircase-ramp.

    Returns
    -------
    None.
    """
    for i in range(len(modelNames)):
        modelName = modelNames[i]
        #Load model
        model = myokit.load_model('.\\Models\\'+modelName+'.mmt')
        
        #Load Protocol
        log = myokit.DataLog.load_csv('.\\Protocols\\'+protocol+'.csv')
        times = log.time()
        voltages = log['voltage']
        
        #Set up simulation
        sim = myokit.Simulation(model)
        sim.set_fixed_form_protocol(times, voltages)
        
        # plot variables
        tmax = (times[-1] + 1)
        plotTimes = np.arange(0,tmax)
        
        if modelName=='larsen2010':
            extraPlots = True
            factorVar = 'iKr_Markov.hERG1b'
            factorVal = [0,20,40,60,80,100]
            factorLabels = ['0% hERG1b', '20% hERG1b', '40% hERG1b', '60% hERG1b', '80% hERG1b', '100% hERG1b']
            factorCount = len(factorVal)
        elif modelName=='herg1a1bIndependence':
            extraPlots = True
            factorVar = 'subunits.a'
            factorVal = [0,1,2,3,4]
            factorLabels = ['0:4','1:3','2:2','3:1','4:0']
            factorCount = len(factorVal)
        else:
            extraPlots = False
            factorCount = 1
        
        if extraPlots:
            currents = np.zeros((len(plotTimes),factorCount))
        else:
            currents = np.zeros((len(plotTimes),1))
        
        plt.figure()
        plt.xlabel('Time (ms)')
        plt.ylabel('Current (nA)')
        plt.title(modelName)
        
        for i in range(factorCount):
            sim.reset()
            if extraPlots:
                sim.set_constant(factorVar, factorVal[i])
            log2 = sim.run(tmax,log_times=plotTimes)
            currents[:,i] = log2['environment.IKr']
            plt.plot(plotTimes,currents[:,i])
        
        if extraPlots:
            plt.legend(factorLabels)
        
        plt.show()
        
    plt.figure()
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (mV)')
    plt.plot(times,voltages)

def currentTraceMultiProtocol(modelName, protocol = 'Additional Protocols\staircase-ramp', factorVar = [], factorVal = []):
    """
    Simulates the model described by 'modelName' under changes specified by factorVar and factorVal for a given protocol. Designed for protocol files with multiple potential functions in the same csv. It will plot model current for each potential on the same graph and then plot the potentials, also on the same graph.
    Parameters
    ----------
    modelName : string
        Name of myokit model file in the Models folder. For example, 'larsen2010' for the Models/larsen2010.mmt model.
    protocol : string
        Filename of a protocol to load from the 'Protocols' folder. For example, 'Clerx 2019\Traditional\Pr2' loads protocol 2 from the 4 ways to fit an ion channel paper (full filename: .\Protocols\Clerx 2019\Traditional\Pr2.csv). Default is the staircase-ramp.
    factorVar : list of strings, optional
        List of names of variables (constants) in the myokit model to change before running. Default is [], ie no changes from the .mmt file
    factorVal : list of values, optional
        List of values to set the variables described in factorVar to in the .mmt file before simulation. Should be the same length as factorVar. Default is [], ie no changes from the .mmt file

    Returns
    -------
    None.
    """
    #Load model
    model = myokit.load_model('.\\Models\\'+modelName+'.mmt')
    
    #Load Protocol
    log = myokit.DataLog.load_csv('.\\Protocols\\'+protocol+'.csv')
    times = log.time()
    voltageKeys = list(log.keys())[1:]
    
    for i in range(len(factorVal)):
        model.set_value(factorVar[i],factorVal[i])
    
    #Set up simulation
    sim = myokit.Simulation(model)
    
    # plot variables
    tmax = (times[-1] + 1)
    plotTimes = np.arange(0,tmax)
    
    currents = np.zeros((len(plotTimes),len(voltageKeys)))
    
    plt.figure()
    plt.xlabel('Time (ms)')
    plt.ylabel('Current (nA)')
    plt.title(modelName)
    for i in range(len(voltageKeys)):
        sim.reset()
        sim.set_fixed_form_protocol(times, log[voltageKeys[i]])
        log2 = sim.run(tmax,log_times=plotTimes)
        currents[:,i] = log2['environment.IKr']
        plt.plot(plotTimes,currents[:,i])
    
    plt.show()
    
    plt.figure()
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (mV)')
    for i in range(len(voltageKeys)):
        plt.plot(times,log[voltageKeys[i]])

#%%
modelNames = ['larsen2010', 'sale2008-1a1b']

currentTracePlots(modelNames, 'Clerx 2019\Traditional\Pr6')

modelName = 'larsen2010'
factorVar = ['iKr_Markov.hERG1b']
factorVal = [40]

currentTraceMultiProtocol(modelName, 'Clerx 2019\Traditional\Pr5', factorVar, factorVal)
currentTraceMultiProtocol('sale2008-1a1b', 'Clerx 2019\Traditional\Pr5')
