# Peak IV curve
#%% Peak IV curve from 4 ways to fit, Pr5
import matplotlib.pyplot as plt
import myokit
import numpy as np

from scipy.optimize import curve_fit

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

def summaryCurves(modelName, factorVar = '', factorVal = []):
    """
    Plots the Peak IV curve for a model using Protocol 5 from 4 ways to fit an ion channel model. It plots a curve for the model specified by 'modelName'. This model can be editted by defining a 'factorVar' variable to change in the myokit model. This variable will be set to the values in 'factorVal' and used to plot multiple IV curves on the same graph.
    It also plots the tau-V curves similarly.

    Parameters
    ----------
    modelName : string
        The name of the myokit model to be used. Example, 'larsen2010', or 'sale2008-1a1b'
    factorVar : string, optional
        The variable in the myokit model that is to be editted. Example, 'iKr_Markov.hERG1b'. Default is '' for no changes to the default model.
    factorVal : list of values, optional
        If the aim is to explore the IV curves when varying a parameter of the model, then this is the list of parameter values to explore. The default is [] for no changes to the default model.

    Returns
    -------
    None.
    """
    
    #Load Protocol
    log = myokit.DataLog.load_csv(r'.\Protocols\Clerx 2019\Traditional\Pr5.csv')
    times = log.time()
    voltageKeys = list(log.keys())[1:]
    
    plt.figure(1)
    plt.xlabel('Voltage (mV)')
    plt.ylabel('Peak current')
    plt.title(modelName)
    
    currents = np.zeros((8662-2662,len(voltageKeys),len(modelName)))
    
    for i in range(len(factorVal)):
        #Load model
        model = myokit.load_model('.\\Models\\'+modelName+'.mmt')
        
        if factorVar != '':
            model.set_value(factorVar,factorVal[i])
        
        #Set up simulation
        sim = myokit.Simulation(model)
        
        # plot variables
        tmax = (times[-1] + 1)
        plotTimes = np.arange(0,tmax)
        
        peaks = np.zeros((len(voltageKeys),1))
        
        for j in range(len(voltageKeys)):
            sim.reset()
            sim.set_fixed_form_protocol(times, log[voltageKeys[j]])
            log2 = sim.run(tmax,log_times=plotTimes)
            peaks[j] = np.max(log2['environment.IKr'][2662:8662]) #Hardcoded range for where voltage varies between protocols
            currents[:,j,i] = log2['environment.IKr'][2662:8662]
        
        plt.plot(np.linspace(-40,-120,len(voltageKeys)), peaks, label = str(factorVal[i])) #Add IV curve for factorVal[i] to the plot
    
    plt.legend()
    plt.show()
    
    #Plot voltage potentional
    plt.figure(2)
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (mV)')
    for i in range(len(voltageKeys)):
        plt.plot(times,log[voltageKeys[i]])
    
    plt.figure(3)
    plt.xlabel('Voltage (mV)')
    plt.ylabel('Time constants - act')
    plt.title(modelName)
    
    plt.figure(4)
    plt.xlabel('Voltage (mV)')
    plt.ylabel('Time constants - rec')
    plt.title(modelName)
    
    voltages = np.linspace(-40,-120,len(voltageKeys))
    for i in range(len(factorVal)):
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
            #plt.title(str(factorVal[i])+'%1b, '+str(voltageKeys[j])+' voltage')
            #plt.legend()
            #plt.show()
        
        plt.figure(3)
        print('--------')
        print(voltages)
        print(tau_act)
        plt.plot(voltages, tau_act)
        plt.figure(4)
        plt.plot(voltages, tau_rec)
    
    plt.figure(3)
    plt.show()
    plt.figure(4)
    plt.show()

modelName = 'larsen2010'
factorVar = 'iKr_Markov.hERG1b'
factorVal = [0, 20, 40, 60, 80, 100]

summaryCurves(modelName, factorVar, factorVal)
