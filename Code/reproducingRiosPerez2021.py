## Reproduce figure 6 of Rios Perez 2021 using the Sale 2009 1a1b model
#%%
# Imports
import matplotlib.pyplot as plt
import myokit
import numpy as np
import scipy


#Function definitions
def f(t,a1,a2,tau1,tau2,c):
    return a1*np.exp(-t/tau1) + a2*np.exp(-t/tau2) + c



#Load model
model1a1b = myokit.load_model('.\Models\sale2008-1a1b.mmt')
model1a1b.set_value('environment.Ko',5.4)
model1a1b.set_value('environment.Ki',150)
model1a1b.set_value('environment.T',273 + 35)

#%% Figure 6A
#Load protocol excel file
log = myokit.DataLog.load_csv(r'.\Protocols\Rios-Perez 2021\riosPerez6ACurrents.csv')

print(log.keys())

protocolTimes = log.time()
tmax = protocolTimes[-1] + 1
plotTimes = np.arange(0,tmax)
voltageKeys = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13']

plt.figure(1)
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.figure(2)
plt.xlabel('Time (ms)')
plt.ylabel('Current')
plt.title('Figure 6A using Sale model')

outputCurrents = np.zeros((len(voltageKeys),len(plotTimes)))
sim = myokit.Simulation(model1a1b)

for i in range(len(voltageKeys)):
    print('------------')
    print('Running')
    print(voltageKeys[i])
    voltages = log[voltageKeys[i]]
    

    sim.set_fixed_form_protocol(protocolTimes, voltages)

    
    log2 = sim.run(tmax,log_times=plotTimes)
    print('Simulation finished')
    plt.figure(1)
    plt.plot(plotTimes, log2['stimulus.IStim'])
    plt.figure(2)
    plt.plot(log2['environment.t'],log2['environment.IKr'])
    outputCurrents[i,:] = log2['environment.IKr']
    sim.reset()

plt.figure(1)
plt.show()
plt.figure(2)
plt.show()

#%% Figures 6D and 6E
deactivationRates = np.zeros((len(voltageKeys),2))
for i in range(len(voltageKeys)):
    #plt.figure(i+3)
    #plt.title(voltageKeys[i])
    xdata = plotTimes[5001:5200]-5000
    ydata = outputCurrents[i,5001:5200]
    print('Identifying deactivation time constants')
    print(voltageKeys[i])
    popt,v,inf,msg,ir=scipy.optimize.curve_fit(f, xdata, ydata, p0 = [-0.08, 0.1, 2, 24, 2e-3],\
                                               maxfev = 1000000, \
                                                   bounds = ([-np.inf,-np.inf,0,0,-1], \
                                                             [np.inf,np.inf,np.inf,np.inf,1]), \
                                                       full_output = True, ftol=1e-14,gtol=-1)
    print('Optimised Parameters')
    print(popt)
    print(msg)
    #plt.plot(xdata,ydata,label='Data')
    #plt.plot(xdata,f(xdata,popt[0],popt[1],popt[2],popt[3],popt[4]),label='Exponential fit')
    #plt.legend()
    #plt.show()
    deactivationRates[i,:] = [max(popt[2:4]), min(popt[2:4])]

plt.boxplot(np.delete(deactivationRates,2,axis=0)[:,0]) #Drop the -60mV to -60mV decay from the plot
plt.title('Slow deactivation rate')
plt.show()

plt.boxplot(np.delete(deactivationRates,2,axis=0)[:,1])
plt.title('Fast deactivation rate')
plt.show()

#%% Figure 6B
#Load protocol excel file
log = myokit.DataLog.load_csv(r'.\Protocols\Rios-Perez 2021\riosPerez6BSteady.csv')

print(log.keys())

protocolTimes = log.time()
tmax = protocolTimes[-1] + 1
plotTimes = np.arange(0,tmax)
voltageKeys = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13']

plt.figure(1)
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.figure(2)
plt.xlabel('Time (ms)')
plt.ylabel('Current')


outputCurrents = np.zeros((len(voltageKeys),len(plotTimes)))
sim.reset()

for i in range(len(voltageKeys)):
    print('------------')
    print('Running')
    print(voltageKeys[i])
    voltages = log[voltageKeys[i]]
    

    sim.set_fixed_form_protocol(protocolTimes, voltages)

    
    log2 = sim.run(tmax,log_times=plotTimes)
    print('Simulation finished')
    plt.figure(1)
    plt.plot(plotTimes, log2['stimulus.IStim'])
    plt.figure(2)
    plt.plot(log2['environment.t'],log2['environment.IKr'])
    outputCurrents[i,:] = log2['environment.IKr']
    sim.reset()

plt.figure(1)
plt.show()
plt.figure(2)
plt.show()

plt.plot(np.arange(-80,50,10), outputCurrents[:,-1])
plt.title('Figure 6B using Sale model')
plt.xlabel('Voltage (mV)')
plt.ylabel('Steady State Current')
plt.show()



#%% Figure 6C
#Load protocol excel file
log = myokit.DataLog.load_csv(r'.\Protocols\Rios-Perez 2021\riosPerez6CTail.csv')

print(log.keys())

protocolTimes = log.time()
tmax = protocolTimes[-1] + 1
plotTimes = np.arange(0,tmax)
voltageKeys = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11']

plt.figure(1)
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.figure(2)
plt.xlabel('Time (ms)')
plt.ylabel('Current')


outputCurrents = np.zeros((len(voltageKeys),len(plotTimes)))
sim.reset()

for i in range(len(voltageKeys)):
    print('------------')
    print('Running')
    print(voltageKeys[i])
    voltages = log[voltageKeys[i]]
    

    sim.set_fixed_form_protocol(protocolTimes, voltages)

    
    log2 = sim.run(tmax,log_times=plotTimes)
    print('Simulation finished')
    plt.figure(1)
    plt.plot(plotTimes, log2['stimulus.IStim'])
    plt.figure(2)
    plt.plot(log2['environment.t'],log2['environment.IKr'])
    outputCurrents[i,:] = log2['environment.IKr']
    sim.reset()

plt.figure(1)
plt.show()
plt.figure(2)
plt.show()

plt.plot(np.arange(-60,50,10), np.amax(outputCurrents[:,4101:],axis=1))
plt.title('Figure 6C using Sale model')
plt.xlabel('Voltage (mV)')
plt.ylabel('Max Tail Current')
plt.show()

