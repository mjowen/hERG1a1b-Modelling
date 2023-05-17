import matplotlib.pyplot as plt
import myokit
import numpy as np

"""
Minor discrepencies in validation. Curves match well. Voltage at maximum current should be -48.4mV for hERG1a and -28.4mV for hERG1b (this code says -48.56mV for hERG1a and -28.28mV for hERG1b) only off by 1 timestep (-48.44 is 1ms before, -28.4 is 1 ms after), probably just a solver detail, maybe in rounding of parameters
""" 

model = myokit.load_model('.\Models\larsen2010.mmt')

log = myokit.DataLog.load_csv('.\Protocols\Larsen 2010\larsenRamp.csv')

print(log.keys())

times = log.time()
voltages = log['voltage']

sim = myokit.Simulation(model)

plt.figure(1)
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.plot(times, voltages,label='Voltage')
plt.show()

hERG1bRatio = [0, 20, 40, 60, 80, 100]
print('-----------')
for i in range(len(hERG1bRatio)):
    sim.set_constant('iKr_Markov.hERG1b',hERG1bRatio[i])
    sim.set_fixed_form_protocol(times, voltages)
    
    tmax = times[-1] + 1
    plotTimes = np.arange(0,tmax,1)
    log2 = sim.run(tmax,log_times=plotTimes)
    
    plt.xlabel('Time (ms)')
    plt.ylabel('Normalised Current (mA)')
    plt.plot(log2['cell.V'],np.divide(log2['IKr.i_Kr'],max(log2['IKr.i_Kr'])))
    plt.xlim((40,-80))
    print('Voltage at max current:')
    print(log2['cell.V'][np.argmax(log2['IKr.i_Kr'])])
    sim.reset()
plt.show()
