import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import csv
if os.getcwd()[-9:] == 'Modelling':
    os.chdir('Code')
import myokitGeneral, myokitFitting
if os.getcwd()[-4:] == 'Code':
    os.chdir('..')

#Options
plotting = False
filenameOut = 'test-002'
maxIter = 2500
iterCount = 5
parallel = False

#Define function to estimate leak current from staircase protocol current
def fit_leak_lr(staircase_protocol, current, V_win=[-115, -85], V_full=[-120, -80],
        ramp_start=0.3, ramp_end=0.7, dt=2e-4):
    # Fitting leak during the first ramp in staircaseramp prt
    #
    # staircase_protocol: full staircase ramp protocol
    # current: corresponding current for the staircase ramp protocol
    # V_win: Voltage window for fitting (in the direction of time)
    # V_full: Full voltage range during the ramp (in the direction of time)
    # ramp_start: starting time of the ramp that matches the input protocol
    # ramp_end: ending time of the ramp that matches the input protocol
    # dt: duration of each time step to work out the index in the input protocol
    from scipy import stats
    rampi, rampf = int(ramp_start / dt), int(ramp_end / dt)
    n_samples = rampf - rampi
    idxi = int(np.abs(V_win[0] - V_full[0])\
            / np.abs(V_full[1] - V_full[0])\
            * n_samples)
    idxf = int(np.abs(V_win[1] - V_full[0])\
            / np.abs(V_full[1] - V_full[0])\
            * n_samples)
    # Assumed V_win, V_full where given correctly!!
    x = staircase_protocol[rampi:rampf][idxi:idxf]
    y = current[rampi:rampf][idxi:idxf]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return slope, -1 * intercept / slope  # g_leak, E_leak

np.random.seed(1)
print(datetime.now())

#Define artefact variables
voffset = np.random.normal(0,1.5)
gLeak = np.random.lognormal(-2,0.3)
gKr = np.random.lognormal(2,0.83255)

print('Voltage Offset (mV): '+str(voffset))
print('Leak conductance (nS): '+str(gLeak))
print('gKr (nS): '+str(gKr))

#Save artefact variables
with open(filenameOut+'-trueParams.csv', 'w', newline = '') as csvfile:
    writer = csv.writer(csvfile, delimiter = ',')
    writer.writerow(['voffset', 'gLeak', 'gKr'])
    writer.writerow([voffset, gLeak, gKr])

#Simulate artefact model
model = myokitGeneral.loadModel('fink2008-artefactModel')
sim = myokitGeneral.compileModel(model, factorVars = ['voltageclamp.voffset_eff', 'voltageclamp.gLeak', 'iKr_Markov.p17'], factorVals = [voffset, gLeak, gKr])
sim, tmax = myokitGeneral.addProtocol(sim, 'Additional Protocols/staircase-ramp')

output = myokitGeneral.simulate(sim, tmax, output = ['environment.Imeasured', 'environment.IKr', 'membrane.V','voltageclamp.Vc'])
noisyOutput = output
noisyOutput[0] = output[0] + np.random.normal(size=len(output[0]))

if plotting:
    plt.plot(output[0])
    plt.title('Measured Current')
    plt.show()
    plt.plot(noisyOutput[0])
    plt.title('Measured Current (noisy)')
    plt.show()
    plt.plot(output[1])
    plt.title('IKr current')
    plt.show()
    plt.plot(output[2])
    plt.title('Membrane Voltage')
    plt.show()

#Find leak current
gLeak_est, ELeak_est = fit_leak_lr(noisyOutput[3],noisyOutput[0],dt=1e-3)
print('Estimated Leak Conductance: '+str(gLeak_est))
print('Estimated Leak Reversal Potential: '+str(ELeak_est))
sim.reset()
sim = myokitGeneral.editSim(sim, factorVars = ['voltageclamp.voffset_eff','voltageclamp.gLeak','iKr_Markov.p17', 'voltageclamp.gLeak_est','voltageclamp.ELeak_est'], factorVals = [voffset, gLeak, gKr, 0.9*gLeak_est, ELeak_est])

#Find leak subtracted current
outputArtefact = myokitGeneral.simulate(sim, tmax, output = ['environment.Imeasured', 'environment.IKr', 'membrane.V','voltageclamp.Vc','voltageclamp.ILeak_est'])
leakSubtractedCurrent = np.subtract(outputArtefact[0],outputArtefact[4])
leakSubtractedCurrent = leakSubtractedCurrent + np.random.normal(size=len(leakSubtractedCurrent))

if plotting:
    plt.plot(leakSubtractedCurrent)
    plt.title('Measured Current: Leak Subtracted')
    plt.show()

#Artefact free data
model = myokitGeneral.loadModel('fink2008-artefactFreeModel')
sim = myokitGeneral.compileModel(model, factorVars = ['iKr_Markov.p17'], factorVals = [gKr])
sim, tmax = myokitGeneral.addProtocol(sim, 'Additional Protocols/staircase-ramp')

outputArtefactFree = myokitGeneral.simulate(sim, tmax, output = ['environment.Imeasured', 'membrane.V'])
outputArtefactFree[0] = outputArtefactFree[0] + np.random.normal(size=len(outputArtefactFree[0]))
#Save artefact free data

if plotting:
    plt.plot(outputArtefactFree[0])
    plt.title('Measured Current: Artefact Free Data')
    plt.show()

with open(filenameOut+'-dataCurrents.csv', 'w', newline = '') as csvfile:
    writer = csv.writer(csvfile, delimiter = ',')
    writer.writerow(['time', 'artefactCurrent', 'artefactFreeCurrent'])
    writer.writerows(tuple(zip(np.arange(tmax), leakSubtractedCurrent, outputArtefactFree[0])))

#%% Fitting
def fittingArtefactModels(modelName, values, paramNames, parameters, localBounds, kCombinations, logTransforms, filename):
    protocol = 'Additional Protocols/staircase-ramp'
    outputName = 'environment.Imeasured'
    model = myokitFitting.Model(modelName, paramNames, protocol, outputName)
    
    times = np.arange(0,tmax)
    
    xbest, fbest, xbests, fbests = myokitFitting.fitting(model, times, values, parameters, iterCount, localBounds, kCombinations, logTransforms, returnAll = True, plotting = plotting, maxIter = maxIter, parallel = parallel)
    
    with open(filename+'-fittedParams.csv', 'w', newline = '') as csvfile:
        writer = csv.writer(csvfile, delimiter = ',')
        writer.writerow(['paramNames','initialParams', 'bestParams-f-'+str(fbest)])
        writer.writerows(tuple(zip(paramNames, parameters, xbest)))
    
    with open(filename+'-posteriors.csv', 'w', newline = '') as csvfile:
        writer = csv.writer(csvfile, delimiter = ',')
        writer.writerow(['paramNames'] + ['params-f-'+str(fbests[i]) for i in range(iterCount)])
        a = [list(i) for i in zip(*xbests)]
        writer.writerows([[paramNames[i]] + a[i] for i in range(len(paramNames))])
    
    return xbest, xbests, fbest, fbests

print('First fitting: Artefact Data and Artefact Model')
modelName = 'fink2008-artefactModel'
values = leakSubtractedCurrent
filenameOutFull = filenameOut+'-fit1'
paramNames = ['iKr_Markov.p1', 'iKr_Markov.p2', 'iKr_Markov.p3', 'iKr_Markov.p5', 'iKr_Markov.p6', 'iKr_Markov.p7', 'iKr_Markov.p8', 'iKr_Markov.p9', 'iKr_Markov.p10', 'iKr_Markov.p11', 'iKr_Markov.p13', 'iKr_Markov.p14', 'iKr_Markov.p15', 'iKr_Markov.p16', 'iKr_Markov.p17', 'voltageclamp.voffset_eff', 'voltageclamp.gLeak']
parameters = [0.20618, 0.0112, 0.04209, 0.02202, 0.0365, 0.41811, 0.0223, 0.13279, -0.0603, 0.08094, 0.00023, -0.0399, 0.04150, -0.0312, 0.1524*1e3, 0, 5]
parameters = np.multiply(parameters,np.random.uniform(0.5,1.5,len(parameters)))
localBounds = [[0,1e-7,1e3], [1,1e-7,0.4], [2,1e-7,1e3], [3,1e-7,1e3], [4,1e-7,0.4], [5,1e-7,1e3], [6,1e-7,0.4], [7,1e-7,1e3], [8,-0.4,-1e-7], [9,1e-7,1e3], [10,1e-7,1e3], [11,-0.4,-1e-7], [12,1e-7,1e3], [13,-0.4,-1e-7]]
kCombinations = [[0,1], [3,4], [5,6], [7,8], [10,11], [12,13]]
logTransforms = [0, 2, 3, 5, 7, 9, 10, 12]
xbest, fbest, xbests, fbests = fittingArtefactModels(modelName, values, paramNames, parameters, localBounds, kCombinations, logTransforms, filenameOutFull)
print('Finished first fitting')
print(datetime.now())

print('Second fitting: Artefact Data and Artefact Free Model')
modelName = 'fink2008-artefactFreeModel'
values = leakSubtractedCurrent
filenameOutFull = filenameOut+'-fit2'
paramNames = ['iKr_Markov.p1', 'iKr_Markov.p2', 'iKr_Markov.p3', 'iKr_Markov.p5', 'iKr_Markov.p6', 'iKr_Markov.p7', 'iKr_Markov.p8', 'iKr_Markov.p9', 'iKr_Markov.p10', 'iKr_Markov.p11', 'iKr_Markov.p13', 'iKr_Markov.p14', 'iKr_Markov.p15', 'iKr_Markov.p16', 'iKr_Markov.p17']
parameters = [0.20618, 0.0112, 0.04209, 0.02202, 0.0365, 0.41811, 0.0223, 0.13279, -0.0603, 0.08094, 0.00023, -0.0399, 0.04150, -0.0312, 0.1524*1e3]
parameters = np.multiply(parameters,np.random.uniform(0.5,1.5,len(parameters)))
localBounds = [[0,1e-7,1e3], [1,1e-7,0.4], [2,1e-7,1e3], [3,1e-7,1e3], [4,1e-7,0.4], [5,1e-7,1e3], [6,1e-7,0.4], [7,1e-7,1e3], [8,-0.4,-1e-7], [9,1e-7,1e3], [10,1e-7,1e3], [11,-0.4,-1e-7], [12,1e-7,1e3], [13,-0.4,-1e-7]]
kCombinations = [[0,1], [3,4], [5,6], [7,8], [10,11], [12,13]]
logTransforms = [0, 2, 3, 5, 7, 9, 10, 12]
xbest, fbest, xbests, fbests = fittingArtefactModels(modelName, values, paramNames, parameters, localBounds, kCombinations, logTransforms, filenameOutFull)
print('Finished second fitting')
print(datetime.now())

print('Third fitting: Artefact Free Data and Artefact Free Model')
modelName = 'fink2008-artefactFreeModel'
values = outputArtefactFree[0]
filenameOutFull = filenameOut+'-fit3'
parameters = [0.20618, 0.0112, 0.04209, 0.02202, 0.0365, 0.41811, 0.0223, 0.13279, -0.0603, 0.08094, 0.00023, -0.0399, 0.04150, -0.0312, 0.1524*1e3]
parameters = np.multiply(parameters,np.random.uniform(0.5,1.5,len(parameters)))
xbest, fbest, xbests, fbests = fittingArtefactModels(modelName, values, paramNames, parameters, localBounds, kCombinations, logTransforms, filenameOutFull)
print('Finished third fitting')
print(datetime.now())
