import matplotlib.pyplot as plt
import numpy as np
from Code.myokitGeneral import *

def simplePlot(modelName, protocol, protocolName = 'voltage', factorVars = [], factorVals = [], plotProtocol=False, xlabel = 'Time (ms)', ylabel = 'Current (nA)', title = ''):
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
    Returns
    -------
    None.
    """
    model = loadModel(modelName)
    sim = compileModel(model, factorVars, factorVals)
    sim, tmax = addProtocol(sim, protocol, protocolName, plotProtocol)
    current = simulate(sim, tmax)[0]
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
    plt.legend = legend
    if len(title)==0:
        title = modelName
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()

def multiSettingPlot(modelName, protocol, protocolName, factorVars = [], factorVals = [], legend = [], plotProtocol = False, xlabel = 'Time (ms)', ylabel = 'Current (nA)', title = ''):
    """
    Plots the output of a model on a single figure for a protocol while varying model settings. 
    Parameters
    ----------
    modelName : string
        model filename
    protocol : string
        Filename/path to a csv protocol file. 
    factorVars : list of strings
        A list of names of varaibles in the model to change
    factorVals : list of floats
        A list of factorVals (a single factorVal is a list of values of varaibles specified in factorVars in the model to change)
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
    sim = compileModel(model, factorVars, factorVals)
    for i in range(len(factorVars)):
        sim.reset()
        sim = editSim(sim, [factorVars[i]], [factorVals[i]])
        if i==0:
            sim, tmax = addProtocol(sim, protocol, protocolName, plotProtocol)
        else:
            sim, tmax = addProtocol(sim, protocol, protocolName)
        current = simulate(sim, tmax)
        plt.plot(np.arange(tmax),current[0])
    if len(legend)==0:
        legend = factorVals
    plt.legend = legend
    if len(title)==0:
        title = modelName
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()