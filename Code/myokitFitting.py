#%%
# Imports
import os
import matplotlib.pyplot as plt
import numpy as np
import pints
if os.getcwd()[-9:] == 'Modelling':
    os.chdir('Code')
import myokitGeneral
if os.getcwd()[-4:] == 'Code':
    os.chdir('..')

class Model(pints.ForwardModel):
    def __init__(self, modelName, paramNames, protocol, outputName):
        self.model = myokitGeneral.loadModel(modelName)
        self.sim = myokitGeneral.compileModel(self.model)
        self.sim, self.tmax = myokitGeneral.addProtocol(self.sim, protocol)
        
        self.paramNames = paramNames
        self.outputName = outputName
        
        self.sim.pre(200)
        
    def n_parameters(self):
        return len(self.paramNames)
    
    def simulate(self, parameters, times):
        # Reset the simulation
        self.sim.reset()
        
        # Update the parameters
        self.sim = myokitGeneral.editSim(self.sim, factorVars = self.paramNames, factorVals = parameters)
        
        # Run a simulation
        output = myokitGeneral.simulate(self.sim, self.tmax, output = [self.outputName])
        
        return output[0]

class AdvancedBoundaries(pints.Boundaries):
    def __init__(self, paramCount, localBounds, kCombinations, vHigh = 40, vLow = -120):
        # localBounds = [[0,1e-7,1e3],[3,1e-3,1e5]] sets the bounds for parameter index 0 to be [1e-7,1e3] and index 3 to be [1e-3,1e5]
        # kCombinations = [[0,1],[4,5]] means param[0]*exp(param[1]*V) and param[4]*exp(param[5]*V) satisfy bounds
        # self.a_min = 1e-7
        # self.a_max = 1e3
        # self.b_min = 1e-7
        # self.b_max = 0.4
        self.km_min = 1.67e-5
        self.km_max = 1e3
        self.vLow = vLow
        self.vHigh = vHigh
        self.paramCount = paramCount
        self.kCombinations = kCombinations
        
        self.lowerBounds = [-np.inf for i in range(paramCount)]
        self.upperBounds = [np.inf for i in range(paramCount)]
        for bound in localBounds:
            self.lowerBounds[bound[0]] = bound[1]
            self.upperBounds[bound[0]] = bound[2]
    
    def n_parameters(self):
        return self.paramCount
    
    def check(self, parameters):
        
        # Check parameter boundaries
        if np.any(parameters <= self.lowerBounds) or np.any(parameters >= self.upperBounds):
            return False
        
        for comb in self.kCombinations:
            kLow = parameters[comb[0]] * np.exp(parameters[comb[1]] * self.vLow)
            kHigh = parameters[comb[0]] * np.exp(parameters[comb[1]] * self.vHigh)
            if kLow < self.km_min or kHigh < self.km_min or kLow > self.km_max or kHigh > self.km_max:
                return False
        
        # All tests passed!
        return True

def fitting(model, times, values, parameters, iterCount, localBounds, kCombinations, logTransforms, maxIter = 10000, returnAll = False, plotting = False):
    problem = pints.SingleOutputProblem(model, times, values)
    error = pints.MeanSquaredError(problem)
    boundaries = AdvancedBoundaries(paramCount = model.n_parameters(), localBounds = localBounds, kCombinations = kCombinations)
    
    for i in range(len(parameters)):
        if i == 0:
            if i in logTransforms:
                transformation = pints.LogTransformation(n_parameters=1)
            else:
                transformation = pints.IdentityTransformation(n_parameters=1)
        else:
            if i in logTransforms:
                transformation = pints.ComposedTransformation(transformation, pints.LogTransformation(n_parameters=1))
            else:
                transformation = pints.ComposedTransformation(transformation, pints.IdentityTransformation(n_parameters=1))
    if plotting:
        plt.figure(figsize=(16, 5))
        plt.xlabel('Time (ms)')
        plt.ylabel('Current (nA)')
        plt.plot(times, values, label='Data')
        plt.plot(times, problem.evaluate(parameters), label='Default parameters')
        plt.legend(loc='upper left')
        plt.show()
    
    xbests = np.zeros((iterCount,model.n_parameters()))
    fbests = np.zeros((iterCount))
    for i in range(iterCount):
        x0 = parameters * 2**np.random.normal(0, 0.5, len(parameters))
        while not boundaries.check(x0):
            x0 = parameters * 2**np.random.normal(0, 0.5, len(parameters))
        # Create an optimisation controller
        opt = pints.OptimisationController(error, x0, boundaries=boundaries, transformation=transformation, method=pints.CMAES)
        opt.set_max_iterations(maxIter)
        # Run the optimisation
        xbest, fbest = opt.run()
        xbests[i,:] = xbest
        fbests[i] = fbest
        
        if plotting:
            # Show the final result
            plt.figure(figsize=(16, 5))
            plt.xlabel('Time (ms)')
            plt.ylabel('Current (nA)')
            plt.plot(times, values, label='Noisy data')
            plt.plot(times, problem.evaluate(xbest), label='Fit')
            plt.legend(loc='upper left')
            plt.show()

    i = np.argmin(fbests)
    xbest = xbests[i,:]
    fbest = fbests[i]

    print('Best rates found. f value of ')
    print(fbest)
    print('found at parameters')
    print(xbest)
    
    if plotting:
        # Show the final result
        plt.figure(figsize=(16, 5))
        plt.xlabel('Time (ms)')
        plt.ylabel('Current (nA)')
        plt.plot(times, values, label='Noisy data')
        plt.plot(times, problem.evaluate(xbest), label='Best fit')
        plt.legend(loc='upper left')
        plt.show()
    if returnAll:
        return xbest, fbest, xbests, fbests
    else:
        return xbest, fbest
