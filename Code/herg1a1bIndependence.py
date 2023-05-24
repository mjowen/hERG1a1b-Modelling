#hERG 1a/1b independence model - Fitting to data and test cases
#%%
# Imports
import matplotlib.pyplot as plt
import myokit
import numpy as np
import pints

class Model(pints.ForwardModel):
    """A two-state Hodgkin-Huxley ion current simulation."""
    def __init__(self):
        # Load model and protocol
        self.model = myokit.load_model(r'.\Models\hERG1a1bIndependence.mmt')
        log = myokit.DataLog.load_csv(r'.\Protocols\Larsen 2010\larsenRamp.csv')
        times = log.time()
        voltages = log['voltage']

        # Create a simulation
        self.sim = myokit.Simulation(self.model)
        self.sim.set_fixed_form_protocol(times, voltages)
        
    def n_parameters(self):
        return 8
    
    def simulate(self, parameters, times):
        # Reset the simulation
        self.sim.reset()
        
        # Update the parameters
        for i, p in enumerate(parameters):
            self.sim.set_constant('aSub.p' + str(1 + i), p)
        
        # Run a simulation
        tmax = times[-1] + 0.1
        log = self.sim.run(tmax, log_times=times, log=['environment.IKr'])
        
        return np.divide(log['environment.IKr'],max(log['environment.IKr']))

# Seed numpy's random generator, just to make sure we get the same example every time
np.random.seed(1)

# Create a model
model = Model()

# Define a parameter vector
parameters = np.array([5e-4, 0.0699, 2e-4, 0.05462, 0.2, 8.91e-3, 0.031, 0.03158])

# Load data csv
log = myokit.DataLog.load_csv(r'.\Data\Larsen 2010\herg1b-0-larsen.csv')
times = log['time'] # Skips evaluation for the first 200ms then tracks 1000ms
values = log['current']

# Set up a problem, and define an error measure
problem = pints.SingleOutputProblem(model, times, values)
error = pints.MeanSquaredError(problem)

# Choose a starting point
x0 = parameters

# Show the intial guess
plt.figure(figsize=(16, 5))
plt.xlabel('Time (ms)')
plt.ylabel('Current (nA)')
plt.plot(times, values, label='Data')
plt.plot(times, problem.evaluate(x0), label='Default parameters')
plt.legend(loc='upper left')
plt.show()

#%% Advanced boundaries
class AdvancedBoundaries(pints.Boundaries):
    
    def __init__(self):
        self.a_min = 1e-7
        self.a_max = 1e3
        self.b_min = 1e-7
        self.b_max = 0.4
        self.km_min = 1.67e-5
        self.km_max = 1e3
        self.v_low = -100
        self.v_high = 40
        
        # Univariate paramater bounds
        self.lower_bounds = np.array([
            self.a_min, self.b_min,
            self.a_min, self.b_min,
            self.a_min, self.b_min,
            self.a_min, self.b_min,
        ])          
        self.upper_bounds = np.array([
            self.a_max, self.b_max,
            self.a_max, self.b_max,
            self.a_max, self.b_max,
            self.a_max, self.b_max,
        ])
        
    def n_parameters(self):
        return 8
    
    def check(self, parameters):
        
        # Check parameter boundaries
        if np.any(parameters <= self.lower_bounds) or np.any(parameters >= self.upper_bounds):
            return False
        
        # Check rate boundaries
        k1m = parameters[0] * np.exp(parameters[1] * self.v_high)
        if k1m <= self.km_min or k1m >= self.km_max:
            return False
        k2m = parameters[2] * np.exp(-parameters[3] * self.v_low)
        if k2m <= self.km_min or k2m >= self.km_max:
            return False
        k3m = parameters[4] * np.exp(parameters[5] * self.v_high)
        if k3m <= self.km_min or k3m >= self.km_max:
            return False
        k4m = parameters[6] * np.exp(-parameters[7] * self.v_low)
        if k4m <= self.km_min or k4m >= self.km_max:
            return False
        
        # All tests passed!
        return True

boundaries = AdvancedBoundaries()

#%% Setup and run optimisiation
np.random.seed(2)

iterCount = 50

transformation = pints.ComposedTransformation(
    pints.LogTransformation(n_parameters=1),       # p1 (a-type)
    pints.IdentityTransformation(n_parameters=1),  # p2 (b-type)
    pints.LogTransformation(n_parameters=1),       # p3 (a-type)
    pints.IdentityTransformation(n_parameters=1),  # p4 (b-type)
    pints.LogTransformation(n_parameters=1),       # p5 (a-type)
    pints.IdentityTransformation(n_parameters=1),  # p6 (b-typ)
    pints.LogTransformation(n_parameters=1),       # p7 (a-type)
    pints.IdentityTransformation(n_parameters=1)   # p8 (b-type)
)

xbests = np.zeros((iterCount,8))
fbests = np.zeros((iterCount,1))

for i in range(iterCount):
    x0 = parameters * 2**np.random.normal(0, 0.5, len(parameters))
    # Create an optimisation controller
    opt = pints.OptimisationController(error, x0, boundaries=boundaries, transformation=transformation, method=pints.CMAES)
    
    # Run the optimisation
    xbest, fbest = opt.run()
    xbests[i,:] = xbest
    fbests[i] = fbest

i = np.argmin(fbests)
xbest = xbests[i,:]
fbest = fbests[i]

print('Best rates found. f value of ')
print(fbest)
print('found at parameters')
print(xbest)
#%% Plot best case

# Show the final result
plt.figure(figsize=(16, 5))
plt.xlabel('Time (ms)')
plt.ylabel('Current (nA)')
plt.plot(times, values, label='Noisy data')
plt.plot(times, problem.evaluate(xbest), label='Best fit')
plt.legend(loc='upper left')
plt.show()

#%% Explore cost function around minimums
# Import PINTS plotting module
import pints.plot
import itertools

for (i,j) in itertools.combinations(range(iterCount),2):
    try:
        # Call "function_between_points"
        fig, axes = pints.plot.function_between_points(
    error, point_1=xbests[i,:], point_2=xbests[j,:], padding=0.5, evaluations=100)
    except:
        print('Failed, moving on')

plt.show()

#%% Plots minimised parameters to check bounds
# Define the boundaries on the parameters and maximum rate coefficients
a_min, a_max = 1e-7, 1e3
b_min, b_max = 1e-7, 0.4
km_min, km_max = 1.67e-5, 1e3

def boundary_plot():

    # Define a range on which to plot the rate coefficient boundaries
    # We use a range that's linear in the log-transformed space
    px = np.exp(np.linspace(np.log(a_min), np.log(a_max), 200))

    # Calculate the lower and upper boundaries on p2 and p4 (which are the same as those on p6 and p8)
    p2_min = (np.log(km_min) - np.log(px)) / 40
    p2_max = (np.log(km_max) - np.log(px)) / 40
    p4_min = (np.log(km_min) - np.log(px)) / 100
    p4_max = (np.log(km_max) - np.log(px)) / 100

    # But p2 and p4 are also bounded by the parameter boundaries, so add that in too:
    p2_min = np.maximum(p2_min, b_min)
    p4_min = np.maximum(p4_min, b_min)

    # Create a figure
    fig = plt.figure(figsize=(10, 4.2))
    fig.subplots_adjust(wspace=0.5)

    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_xlabel('p1 and p5')
    ax1.set_ylabel('p2 and p6')
    ax1.set_xscale('log')
    ax1.set_xlim(3e-8, 3e3)
    ax1.set_ylim(-0.02, 0.42)
    ax1.axvline(a_min, color='#bbbbbb')
    ax1.axvline(a_max, color='#bbbbbb')
    ax1.axhline(b_min, color='#bbbbbb')
    ax1.axhline(b_max, color='#bbbbbb')
    ax1.plot(px, p2_min)
    ax1.plot(px, p2_max)
    ax1.fill_between(px, p2_min, p2_max, color='#dddddd')

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_xlabel('p3 and p7')
    ax2.set_ylabel('p4 and p8')
    ax2.set_xscale('log')
    ax2.set_xlim(3e-8, 3e3)
    ax2.set_ylim(-0.02, 0.42)
    ax2.axvline(a_min, color='#bbbbbb')
    ax2.axvline(a_max, color='#bbbbbb')
    ax2.axhline(b_min, color='#bbbbbb')
    ax2.axhline(b_max, color='#bbbbbb')
    ax2.plot(px, p4_min)
    ax2.plot(px, p4_max)
    ax2.fill_between(px, p4_min, p4_max, color='#dddddd')
    
    return ax1, ax2

# Create a plot of the boundaries
ax1, ax2 = boundary_plot()

# Select a point within the boundaries
p_good = np.array([1e-2, 0.15, 1e-2, 0.05, 1e-2, 0.15, 1e-2, 0.05])

# Create a method that draws a point in p1/p2 space
def check_p1p2(p1, p2, foract):
    p = np.copy(p_good)
    p[0] = p1
    p[1] = p2
    if foract:
        ax1.plot(p1, p2, 'o', color='tab:blue')
    else:
        ax1.plot(p1, p2, 'x', color='tab:red')

# Create a method that draws a point in p3/p4 space
def check_p3p4(p3, p4, foract):
    p = np.copy(p_good)
    p[2] = p3
    p[3] = p4
    if foract:
        ax2.plot(p3, p4, 'o', color='tab:blue')
    else:
        ax2.plot(p3, p4, 'x', color='tab:red')

for i in range(iterCount):
    check_p1p2(xbests[i,0],xbests[i,1], True)
    check_p3p4(xbests[i,2],xbests[i,3], True)
    check_p1p2(xbests[i,4],xbests[i,5], False)
    check_p3p4(xbests[i,6],xbests[i,7], False)

plt.show()

#%% Posterior distributions
fig = plt.figure(figsize=(10, 4.2))
fig.subplots_adjust(wspace=0.5)

for i in range(8):
    ax = fig.add_subplot(2, 4, i+1)
    ax.hist(xbests[:,i])
    ax.axvline(xbest[i],color='black')
    ax.axvline(parameters[i],color='red')
    

