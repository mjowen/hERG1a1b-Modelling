import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

if os.getcwd()[-9:] == 'Modelling':
    os.chdir('Code')
from myokitGeneral import *
if os.getcwd()[-4:] == 'Code':
    os.chdir('..')

def simplePlot(modelName, protocol, protocolName = 'voltage', factorVars = [], factorVals = [], plotProtocol=False, xlabel = 'Time (ms)', ylabel = 'Current (nA)', title = '', output = 'environment.IKr', initialValues = []):
    """
    Plots a current trace predicted by a model for a given protocol. Model parameters can be changed by using factorVars and factorVals
    Parameters
    ----------
    modelName : string
        filname of Myokit model
    protocol : string
        A filename/path to a csv protocol.
    protocolName : string, optional
        Name of protocol in the csv to use. Only relevent if the csv defines multiple protocols. Default is 'voltage', the case for single potential protocols
    factorVars : list of strings, optional
        A list of names of varaibles in the model to change. Default is []
    factorVals : list of floats, optional
        A list of floats corresponding to the variables in factorVars. Default is []
    plotProtocol : bool, optional
        Whether or not to plot the protocol in an additional plot. Default is False
    xlabel : string
        Label for the x-axis of the plot. Default is Time (ms)
    ylabel : string
        Label for the y-axis of the plot. Default is Current (nA)
    title : string
        Label for the title of the plot. Default is modelName
    output : string
        Model variable to plot. Default is 'environment.IKr'
    Returns
    -------
    None.
    """
    model = loadModel(modelName)
    sim = compileModel(model, factorVars, factorVals, initialValues = initialValues)
    sim, tmax = addProtocol(sim, protocol, protocolName, plotProtocol)
    current = simulate(sim, tmax, output=[output])[0]
    plt.figure()
    plt.plot(np.arange(tmax),current)
    if len(title)==0:
        title = modelName
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()

def multiModelPlot(modelNames, protocol, protocolName = 'voltage', factorVars = [], factorVals = [], legend = [], plotProtocol = False, xlabel = 'Time (ms)', ylabel = 'Current (nA)', title = ''):
    """
    Plots multiple models on a single figure for a single protocol. The protocolName can be specified in the given protocol features multiple voltage potentials. The model settings can be changed with factorVars and factorVals
    Parameters
    ----------
    modelNames : list of strings
        list of model filenames
    protocol : string
        Filename/path to a csv protocol file
    protocolName : string, optional
        Name of protocol in the csv to use. Only relevent if the csv defines multiple protocols.
    factorVars : list of lists of strings
        A list of factorVars (a single factorVar is a list of names of varaibles in the model to change)
    factorVals : list of lists of floats
        A list of factorVals (a single factorVal is a list of values of varaibles specified in factorVars in the model to change)
    legend : list of strings, optional
        Legend names to use in the plot. Default is the labels of the potential in the csv protocol.
    plotProtocol : bool, optional
        Whether or not to plot the protocol in an additional plot. Default is False
    xlabel : string
        Label for the x-axis of the plot. Default is Time (ms)
    ylabel : string
        Label for the y-axis of the plot. Default is Current (nA)
    title : string
        Label for the title of the plot. Default is modelName
    Returns
    -------
    None.
    """
    models = loadModels(modelNames)
    plt.figure()
    if len(factorVars)==0:
        factorVars = [[]]*len(modelNames)
        factorVals = factorVars
    for i in range(len(modelNames)):
        sim = compileModel(models[i], factorVars[i], factorVals[i])
        sim, tmax = addProtocol(sim, protocol, protocolName)
        current = simulate(sim, tmax)[0]
        plt.plot(np.arange(tmax),current)
    if len(legend)==0:
        legend = modelNames
    plt.legend(legend)
    if len(title)==0:
        title = modelName
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()
    if plotProtocol:
        sim, tmax = addProtocol(sim, protocol, protocolName, plotProtocol)

def multiProtocolPlot(modelName, protocol, factorVars = [], factorVals = [], legend = [], xlabel = 'Time (ms)', ylabel = 'Current (nA)', title = ''):
    """
    Plots the output of a model on a single figure for protocols which specify multiple voltage potentials. The model settings can be changed with factorVars and factorVals.
    Parameters
    ----------
    modelName : string
        model filename
    protocol : string
        Filename/path to a csv protocol file. Should be a file containing multiple voltage potentials.
    protocolName : string, optional
        Name of protocol in the csv to use. Only relevent if the csv defines multiple protocols.
    factorVars : list of strings
        A list of names of varaibles in the model to change. Each setting in factorVars will be looped through to produce a seperate curve.
    factorVals : list of floats
        A list of factorVals (a single factorVal is a list of values of varaibles specified in factorVars in the model to change)
    legend : list of strings, optional
        Legend names to use in the plot. Default is modelNames.
    xlabel : string
        Label for the x-axis of the plot. Default is Time (ms)
    ylabel : string
        Label for the y-axis of the plot. Default is Current (nA)
    title : string
        Label for the title of the plot. Default is modelName
    Returns
    -------
    None.
    """
    model = loadModel(modelName)
    plt.figure()
    sim = compileModel(model, factorVars, factorVals)
    voltageKeys = getProtocolKeys(protocol)
    for i in range(len(voltageKeys)):
        sim.reset()
        sim, tmax = addProtocol(sim, protocol, protocolName=voltageKeys[i])
        current = simulate(sim, tmax)
        plt.plot(np.arange(tmax),current[0])
    if len(legend)==0:
        legend = voltageKeys
    plt.legend(legend)
    if len(title)==0:
        title = modelName
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()

def multiSettingPlot(modelName, protocol, protocolName = 'voltage', factorVars = [[]], factorVals = [[]], legend = [], plotProtocol = False, xlabel = 'Time (ms)', ylabel = 'Current (nA)', title = '', output = 'environment.IKr'):
    """
    Plots the output of a model on a single figure for a protocol while varying model settings. 
    Parameters
    ----------
    modelName : string
        model filename
    protocol : string
        Filename/path to a csv protocol file. 
    factorVars : list of list of strings
        A list of factorVars (a single factorVar is a list of names of variables in the model to change)
    factorVals : list of floats
        A list of factorVals (a single factorVal is a list of values of variables specified in factorVars in the model to change)
    legend : list of strings, optional
        Legend names to use in the plot. Default is factorVals.
    plotProtocol : bool, optional
        Whether or not to plot the protocol in an additional plot. Default is False
    xlabel : string
        Label for the x-axis of the plot. Default is Time (ms)
    ylabel : string
        Label for the y-axis of the plot. Default is Current (nA)
    title : string
        Label for the title of the plot. Default is modelName
    Returns
    -------
    None.
    """
    model = loadModel(modelName)
    plt.figure()
    sim = compileModel(model)
    for i in range(len(factorVars)):
        sim.reset()
        sim = editSim(sim, factorVars[i], factorVals[i])
        sim, tmax = addProtocol(sim, protocol, protocolName)
        current = simulate(sim, tmax, output = [output])
        plt.plot(np.arange(tmax),current[0])
    if len(legend)==0:
        legend = factorVals
    plt.legend(legend)
    if len(title)==0:
        title = modelName
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()
    if plotProtocol:
        sim, tmax = addProtocol(sim, protocol, protocolName, plotProtocol)

# Single exponential
def single(t, a, b, c):
    if c < 1:
        return np.ones(len(t)) * float('inf')
    return a + b * np.exp(-t / c)

# Double exponental
def f(t, a, b1, c1, b2, c2):
    if c1 < 1 or c2 < c1:
        return np.ones(len(t)) * float('inf')
    if b1 * b2 > 0:
        return np.ones(len(t)) * float('inf')
    return a + b1 * np.exp(-t / c1) + b2 * np.exp(-t / c2)

def summaryCurves(modelName, factorVars = [[]], factorVals = [[]], legend = []):
    """
    Plots the Peak IV curve for a model using Protocol 5 from 4 ways to fit an ion channel model. It plots a curve for the model specified by 'modelName'. This model can be editted by defining a 'factorVars' variable to change in the myokit model. This variable will be set to the values in 'factorVals' and used to plot multiple IV curves on the same graph.
    It also plots the tau-V curves similarly.
    Parameters
    ----------
    modelName : string
        The name of the myokit model to be used. Example, 'larsen2010', or 'sale2008-1a1b'
    factorVars : list of lists of strings
        A list of factorVars (a single factorVar is a list of names of variables in the model to change)
    factorVals : list of lists of floats
        A list of factorVals (a single factorVal is a list of values of variables specified in factorVars in the model to change)
    legend : list of strings
        Legend names to use in the plot. Default is factorVals.
    Returns
    -------
    None.
    """
    protocol = 'Clerx 2019\Traditional\Pr5'
    voltageKeys = getProtocolKeys(protocol)
    if len(legend)==0:
        legend = [str(a) for a in factorVals]
    plt.figure(1)
    plt.xlabel('Voltage (mV)')
    plt.ylabel('Peak current')
    plt.title(modelName)
    
    currents = np.zeros((8662-2662,len(voltageKeys),len(factorVals))) #Times are hardcoded for Pr5
    model = loadModel(modelName)
    sim = compileModel(model)
    for i in range(len(factorVals)):
        sim = editSim(sim, factorVars[i],factorVals[i])
        
        peaks = np.zeros((len(voltageKeys),1))
        
        for j in range(len(voltageKeys)):
            sim.reset()
            sim, tmax = addProtocol(sim, protocol, protocolName = voltageKeys[j])
            
            tmp = simulate(sim, tmax)[0]
            peaks[j] = np.max(tmp[2662:8662]) #Hardcoded range for where voltage varies between protocols
            currents[:,j,i] = tmp[2662:8662]
        plt.plot(np.linspace(-40,-120,len(voltageKeys)), peaks) #Add IV curve for factorVals[i] to the plot
    
    plt.legend(legend)
    plt.show()
    
    plt.figure(2)
    plt.xlabel('Voltage (mV)')
    plt.ylabel('Time constants - act')
    plt.title(modelName)
    
    plt.figure(3)
    plt.xlabel('Voltage (mV)')
    plt.ylabel('Time constants - rec')
    plt.title(modelName)
    
    
    plotTimes = np.arange(tmax)
    
    voltages = np.linspace(-40,-120,len(voltageKeys))
    for i in range(len(factorVals)):
        # Determine time constants
        tau_rec = []
        tau_act = []
        for j in range(len(voltageKeys)):
            # Fit exponential
            t = plotTimes[2662:8662] - plotTimes[2661]
            c = currents[:,j,i]
    
            # Guess some exponential parameters
            peak = np.argmax(np.abs(c))
    
            # Deactivation pre-fit, for guess
            guess = 200 if voltages[i] < -60 else 2000
            a2, b2, c2 = 0, c[peak], guess
            popt, pcov = curve_fit(single, t[peak:], c[peak:], [a2, b2, c2], maxfev=10000)
            a2, b2, c2 = popt
    
            # Recovery pre-fit, for guess
            guess = 5
            a1, b1, c1 = c[peak], -c[peak], guess
            if peak < 3:
                # Very fast: Only happens for simulations
                peak = 3
                a1, b1, c1 = -3, 3, 0.1
            try:
                popt, pcov = curve_fit(single, t[:peak], c[:peak], [a1, b1, c1], maxfev=10000)
                a1, b1, c1 = popt
            except RuntimeError:
                pass
    
            # Double exponential
            popt, pcov = curve_fit(f, t, c, [a1, b1, c1, b2, c2])
            tau_rec.append(popt[2])
            tau_act.append(popt[4])
            #plt.figure(5+i*len(voltageKeys)+j)
            #plt.plot(t,c,label='Data')
            #plt.plot(t,f(t,popt[0],popt[1],popt[2],popt[3],popt[4]),label='Exponential fit')
            #plt.title(str(factorVals[i])+'%1b, '+str(voltageKeys[j])+' voltage')
            #plt.legend()
            #plt.show()
        
        plt.figure(2)
        plt.plot(voltages, tau_act)
        plt.legend(legend)
        plt.figure(3)
        plt.plot(voltages, tau_rec)
        plt.legend(legend)
    
    plt.figure(2)
    plt.show()
    plt.figure(3)
    plt.show()