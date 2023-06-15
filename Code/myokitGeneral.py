import matplotlib.pyplot as plt
import myokit
import numpy as np

def loadModels(modelNames):
    """
    Finds the Myokit models given by the filenames in modelNames and returns them
    Parameters
    ----------
    modelNames : list of strings
        The filenames of the models to load.
    Returns
    -------
    models : list of myokit model objects
    """
    models = []
    for i in range(len(modelNames)):
        models.append(myokit.load_model('.\\Models\\'+modelNames[i]+'.mmt'))
    return models

def loadModel(modelName):
    """
    Finds the Myokit model given by the filename in modelName and returns it
    Parameters
    ----------
    modelName : strings
        The filename of the model to load.
    Returns
    -------
    model : myokit model object
    """
    model = myokit.load_model('.\\Models\\'+modelName+'.mmt')
    return model

def compileModel(model, factorVars = [], factorVals = []):
    """
    Generates a Myokit sim object from model after changing the factorVars variables/constants to the values given in factorVals
    Parameters
    ----------
    model : Myokit model object
        The Myokit model object to compile
    factorVars : list of strings, optional
        A list of names of varaibles to change in the model file before compiling. Default is []
    factorVals : list of floats, optional
        A list of floats corresponding to the variables in factorVars. Default is []
    Returns
    -------
    sim : Myokit simulation object
    """
    for i in range(len(factorVals)):
        model.set_value(factorVars[i],factorVals[i])
    sim = myokit.Simulation(model)
    return sim

def editSim(sim, factorVars, factorVals):
    """
    Changes the constants in the sim Myokit simulation object defined by factorsVars to the values defined in factorVals
    Parameters
    ----------
    model : Myokit model object
        The Myokit model object to compile
    factorVars : list of strings
        A list of names of varaibles to change in the sim object
    factorVals : list of floats
        A list of floats corresponding to the variables in factorVars
    Returns
    -------
    sim : Myokit simulation object
    """
    for i in range(len(factorVals)):
        sim.set_constant(factorVars[i], factorVals[i])
    return sim

def addProtocol(sim, protocol, protocolName = 'voltage', plotProtocol = False):
    """
    Adds a csv protocol to a sim object.
    Parameters
    ----------
    sim : Myokit simulation object
        The sim object to add the protocol to.
    protocol : string
        The filename/path for a csv protocol.
    protocolName : string, optional
        The voltage key in the protocol's csv file. Either 'voltage' for single voltage protocols or numbered from 1 upwards for multi-voltage protocols.
    plotProtocol : bool, optional
        Whether or not to plot the protocol in an additional plot. Default is False
    Returns
    -------
    sim : Myokit simulation object
        The sim object with the requested voltage attached
    tmax : float
        The end time for the protocol
    """
    log = myokit.DataLog.load_csv('.\\Protocols\\'+protocol+'.csv')
    times = log.time()
    voltages = log[protocolName]
    tmax = times[-1]
    if plotProtocol:
        plt.plot(times,voltages)
    sim.set_fixed_form_protocol(times, voltages)
    return sim, tmax

def getProtocolKeys(protocol):
    """
    Get the dictionary keys for the voltage
    Parameters
    ----------
    protocol : string
        The filename/path for a csv protocol.
    Returns
    -------
    voltageKeys : list of strings
        The dictionary keys for the voltage, as defined in the protocol csv
    """
    log = myokit.DataLog.load_csv('.\\Protocols\\'+protocol+'.csv')
    voltageKeys = list(log.keys())[1:]
    return voltageKeys

def simulate(sim, tmax, plotTimes = [], output = ['environment.IKr']):
    """
    Runs the inputted Myokit simulation object for the specified tmax end time, returns the requested output variables from the simulation, recorded at the specified plotTimes
    Parameters
    ----------
    sim : Myokit simulation object
        Simulation object to run. Should already have the attached protocol and any requsted changes.
    tmax : float
        The end time of the simulation/protocol
    plotTimes : numpy array, optional
        The times to record outputs, default is 1ms intervals between 0 ms and "tmax" ms
    output : list of strings, optional
        A list of strings for output variables. The default is ['environment.IKr'] ie only current.
    Returns
    -------
    log[output] : array
        Array of output data. The order and variables are defined by output
    """
    if len(plotTimes) == 0:
        plotTimes = np.arange(0,tmax)
    log = sim.run(tmax, log_times=plotTimes)
    return [log[a] for a in output]
