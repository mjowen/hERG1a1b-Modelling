import matplotlib.pyplot as plt
import myokit
import numpy as np

model1a1b = myokit.load_model('.\Models\sale2008-1a1b.mmt')
model1a = myokit.load_model('.\Models\sale2008-1a.mmt')

log = myokit.DataLog.load_csv('.\Protocols\Sale 2008\AP Protocol.csv')

print(log.keys())

times = log['time']
voltages = log['voltage']

sim = myokit.Simulation(model1a1b)
sim2 = myokit.Simulation(model1a)


sim.set_fixed_form_protocol(times, voltages)
sim2.set_fixed_form_protocol(times, voltages)

tmax = times[-1] + 1
log2 = sim.run(tmax,log_times=times)
log3 = sim2.run(tmax,log_times=times)


voltages = np.subtract(voltages, min(voltages))

plt.figure()
plt.xlabel('Time (ms)')
plt.ylabel('Normalised voltage and current')
plt.plot(times, np.divide(voltages,max(voltages)),label='Voltage')
plt.plot(log2['environment.t'],np.divide(log2['environment.IKr'],max(log2['environment.IKr'])))
plt.plot(log3['environment.t'],np.divide(log3['environment.IKr'],max(log2['environment.IKr'])))

#%% State occupancy for Sale
from Code import myokitPlot
modelName = 'sale2008-1a1b'
protocol = 'Perissinotti 2018\deactivation'
factorVar = ['environment.e4031']
factorVal = [0]

myokitPlot.simplePlot(modelName = modelName, protocol=protocol, protocolName = '1', factorVars = factorVar, factorVals = factorVal, plotProtocol = True)
states = ['states.c3', 'states.c2', 'states.c1', 'states.o', 'states.i', 'states.b']
for i in range(len(states)):
    myokitPlot.simplePlot(modelName = modelName, protocol=protocol, protocolName = '1', factorVars = factorVar, factorVals = factorVal, plotProtocol = False, ylabel = states[i], output = states[i])
