from scipy.integrate import odeint
from scipy.optimize import differential_evolution
from scipy.optimize import Bounds
from joblib import Parallel
from joblib import delayed
import numpy as np
import pandas as pd

#sub-Microbial model
from subMFermGroups import subMFermGroups as model

#Import data to confront with
#Y = pd.read_csv("/mnt/9C64E5F864E5D554/Documents and Settings/capekp00/Documents/data_statistika/Junior/Inkubace/IncubationExperiment/Python/pysourcefile/Y.csv")
Y = pd.read_csv("/home/pcapek/IncubationExperiment/Python/pysource/Y.csv")
#Import initial variables (Fungal PLFA, bacterial PLFA, Cflush and pH)
#incond = pd.read_csv("/mnt/9C64E5F864E5D554/Documents and Settings/capekp00/Documents/data_statistika/Junior/Inkubace/IncubationExperiment/Python/pysourcefile/InitialConditions.csv")
incond = pd.read_csv("/home/pcapek/IncubationExperiment/Python/pysource/InitialConditions.csv")
incond = np.array(incond)
#Import H+ ions derived from CO2 dissolution
#HExtra = pd.read_csv("/mnt/9C64E5F864E5D554/Documents and Settings/capekp00/Documents/data_statistika/Junior/Inkubace/IncubationExperiment/Python/pysourcefile/H.csv")
HExtra = pd.read_csv("/home/pcapek/IncubationExperiment/Python/pysource/H.csv")
HExtra = np.array(HExtra)

#Importing data frame with rows indicating parameters and their possible combinations, which may differ between fungi and bacteria 
vpars = pd.read_csv("/home/pcapek/IncubationExperiment/Python/pysource/vpars3free.csv",
                    usecols=range(1, 9))
'''
vpars = pd.read_csv("/mnt/9C64E5F864E5D554/Documents and Settings/capekp00/Documents/data_statistika/Junior/Inkubace/IncubationExperiment/Python/pysourcefile/vpars3free.csv",
                    usecols=range(1, 9))
'''

#Those are all model parameters - lower and upper bounds are defined
ParmsFermGroups = pd.DataFrame({'first': [6, 50, 50, 3, 1.6, 2, 0.5, 0.001, 0.2, 0.4],
                                'lower': [0.1, 1, 1, 0.5, 0.8, 0.05, 0.1, 1e-8, 1e-3, 1e-3],
                                'upper': [20, 500, 500, 10, 10, 3, 5, 1, 1, 1]})

#===================Defining lower and upper bounds of model parameters, which will be optimized - different length of parameter space will be applied so correct reordering and sorting of parameters is essential 
def getBounds(ID):
    #==========Lower bound
    x0L = np.array(ParmsFermGroups.loc[:, 'lower'])
    x0L = pd.Series(x0L, index=['Im', 'KmG', 'KmA', 'Gm', 'edemand', 'psi', 'eta', 'k', 'nb', 'ng'])
    #===============================
    #parameters nb, ng and a are always fixed so they are removed here but will be attached later
    lastxL = x0L[['nb', 'ng']]
    firstxL = x0L.drop(lastxL.index)
    #===============================
    #Those are parameters, which will be assumed to differ between fungi and bacteria
    xLfree = x0L[vpars.loc[ID, :].dropna()]
    #The rest of parameters is, thus fixed
    fixedL = firstxL.index.difference(xLfree.index, sort=False)
    xLfixed = x0L[fixedL]
    #now, the lower bound of parameters could be added together - 2times free parameters appended by fixed parameters and last parameters
    LowerBound = np.concatenate([xLfree, xLfree, xLfixed, lastxL]) #
    #==========Upper bound
    x0U = np.array(ParmsFermGroups.loc[:, 'upper'])
    x0U = pd.Series(x0U, index=['Im', 'KmG', 'KmA', 'Gm', 'edemand', 'psi', 'eta', 'k', 'nb', 'ng'])
    #===============================
    #parameters nb, ng and a are always fixed so they are removed here but will be attached later
    lastxU = x0U[['nb', 'ng']]
    firstxU = x0U.drop(lastxU.index)
    #===============================
    #Those are parameters, which will be assumed to differ between fungi and bacteria
    xUfree = x0U[vpars.loc[ID, :].dropna()]
    #The rest oU parameters is, thus fixed
    fixedU = firstxU.index.difference(xUfree.index, sort=False)
    xUfixed = x0U[fixedU]
    #now, the lower bound of parameters could be added together - 2times free parameters appended by fixed parameters and last parameters
    UpperBound = np.concatenate([xUfree, xUfree, xUfixed, lastxU]) #
    return Bounds(lb=LowerBound, ub=UpperBound)

#===================Defining cost function
def cost(x, Y, incond, HExtra, ID):
    n0 = pd.Series(np.array([1,1,1,1,1,1,1,1]), index = ['Im', 'KmG', 'KmA', 'Gm', 'edemand', 'psi', 'eta', 'k'])
    #Conversion factors
    conversions = x[(len(x)-2):len(x)]
    #intial fungal abundance
    a = incond[4]
    free0 = vpars.loc[ID, :].dropna()
    #Fungal-specific parameters
    fungalP = pd.Series(x[0:len(free0)], index=free0)
    #Bacteria-specific parameters
    bacterialP = pd.Series(x[len(free0):(len(free0)*2)], index=free0)
    #Fixed parameters
    fixedI = n0.index.difference(fungalP.index, sort=False)
    fixedP = pd.Series(x[(len(free0)*2):(len(x)-2)], index=fixedI)
    #===============================================
    fungalAll = pd.concat([fungalP, fixedP])
    fungalAll = fungalAll.reindex(['Im', 'KmG', 'KmA', 'Gm', 'edemand', 'psi', 'eta', 'k'])
    bacterialAll = pd.concat([bacterialP, fixedP])
    bacterialAll = bacterialAll.reindex(['Im', 'KmG', 'KmA', 'Gm', 'edemand', 'psi', 'eta', 'k'])
    #===============================================
    #input to odeint - first all fungal, then all bacterial. !!!fixed parameters are duplicated because they are independent of microbial groups
    p0 = [np.concatenate([fungalAll, bacterialAll])]
    #initial concentration of H+ ions
    H0 = 10**(-incond[3])
    #initial condition
    Bn0 = incond[2]/conversions[0]
    Bf0 = Bn0*a
    Bb0 = Bn0*(1 - a)
    y0 = [500, 0, 0, 0, 0, 0, 0, 0, Bf0.item(), Bb0.item(), 0]
    #additional conversion factors for fungal and bacterial PLFA 
    conversions = np.concatenate([conversions, incond[0]/Bf0, incond[1]/Bb0])

    #time steps
    t = pd.to_numeric(Y.Time, errors='coerce')

    #Solving the model
    Yhat0 = odeint(model, y0=y0, args=tuple(p0), t=t)

    #Reformating simulations
    Yhat = np.concatenate((Yhat0[:, 0].reshape(1, len(t)), #Glucose
                       Yhat0[:, 7].reshape(1, len(t)),#labelled CO2
                       (Yhat0[:, 3]*conversions[2]).reshape(1, len(t)),
                       (Yhat0[:, 6]*conversions[3]).reshape(1, len(t)),
                       ((Yhat0[:, 3] + Yhat0[:, 6])*conversions[0] + 
                        ((Yhat0[:, 3] + Yhat0[:, 8])*(Yhat0[:, 1] + Yhat0[:, 2])+
                         (Yhat0[:, 6] + Yhat0[:, 9])*(Yhat0[:, 4] + Yhat0[:, 5]))*conversions[1]).reshape(1, len(t)),
                         (-np.log10(H0 + HExtra[:, 0] + ((Yhat0[:, 10]/2*0.28/(1-0.28)/1e3)*1.75e-5/(HExtra[:, 1] + 1.75e-5)))).reshape(1, len(t))), axis = 0)

    #Re-formating observations
    Y = np.array(Y.drop('Time', axis=1)).T
    #Weights of residuals 
    Wmin = np.array([])
    Wmax = np.array([])
    for i in range(Y.shape[0]):
        Wmin = np.append(Wmin, min(Y[i, ~np.isnan(Y[i, :])]))
        Wmax = np.append(Wmax, max(Y[i, ~np.isnan(Y[i, :])]))

    Wmin = np.repeat(Wmin, Y.shape[1]).reshape(Y.shape)
    Wmax = np.repeat(Wmax, Y.shape[1]).reshape(Y.shape)
    W = Wmax - Wmin

    #Residuals
    residuals = (((Yhat - Y)/W)**2).reshape(1, (Y.shape[0]*Y.shape[1]))
    #sum of residuals
    out = sum(residuals[~np.isnan(residuals)])
    return out

'''
#this part of script can be used to do single calculations  
bounds = getBounds(ID=91)
result = differential_evolution(cost, bounds, args=(Y, incond, HExtra, 91))
np.array(result.x).tofile("/mnt/9C64E5F864E5D554/Documents and Settings/capekp00/Documents/data_statistika/Junior/Inkubace/IncubationExperiment/Python/pyresults/optpars92.csv", sep = ";")
'''
'''
#this part of script is used for parallel calculation across rows of vpars
def MinimizeParallel(ID):
    bounds = getBounds(ID)
    result = differential_evolution(cost, bounds, args=(Y, incond, HExtra, ID))
    fileout = "/home/pcapek/IncubationExperiment/Python/pyresults/" + "optpars" + str((ID+1)) + ".csv"
    #fileout = "/mnt/9C64E5F864E5D554/Documents and Settings/capekp00/Documents/data_statistika/Junior/Inkubace/IncubationExperiment/Python/pyresults/" + "optpars" + str((ID+1)) + ".csv"
    result.x.tofile(fileout, sep = ';')

#Parallel(n_jobs=2)(delayed(MinimizeParallel)(i) for i in [0, 1])

Parallel(n_jobs=vpars.shape[0])(delayed(MinimizeParallel)(i) for i in range(0, vpars.shape[0]))
'''

#This part of script is used to calculate best parameters 100 times to know the variability in parameters estimates
def MinimizeParallel(ID):
    bounds = getBounds(6)
    result = differential_evolution(cost, bounds, args=(Y, incond, HExtra, 6))
    fileout = "/home/pcapek/IncubationExperiment/Python/pyresults/" + "PLFinalpars" + str(ID) + ".csv"
    #fileout = "/mnt/9C64E5F864E5D554/Documents and Settings/capekp00/Documents/data_statistika/Junior/Inkubace/IncubationExperiment/Python/pyresults/" + "optpars" + str((ID+1)) + ".csv"
    result.x.tofile(fileout, sep = ';')

Parallel(n_jobs=100)(delayed(MinimizeParallel)(i) for i in range(1, 101))


'''!!!NOT TO RUN
#======================Lets run 2 examples sequentially and in parallel
Parallel(n_jobs=2)(delayed(MinimizeParallel)(i) for i in [0, 1])


for i in [0, 1]:
    MinimizeParallel(i)
'''
'''
#initialize pandas dictionary
ncols = range(1, 2)
out = {V: [] for V in ncols}

bounds = getBounds(ID=0)
x = bounds.lb

x
print(fileout)
print(vpars.loc[89, :])
'''