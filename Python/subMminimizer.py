from scipy.integrate import odeint
from scipy.optimize import differential_evolution
from scipy.optimize import Bounds
import numpy as np
import pandas as pd

#sub-Microbial model
from subMFermGroups import subMFermGroups as model

#Import data to confront with
Y = pd.read_csv("/mnt/9C64E5F864E5D554/Documents and Settings/capekp00/Documents/data_statistika/Junior/Inkubace/IncubationExperiment/Python/pysourcefile/Y.csv")
#Y = pd.read_csv("/home/pcapek/IncubationExperiment/Python/pysource/Y.csv")
#Import initial variables (Fungal PLFA, bacterial PLFA, Cflush and pH)
incond = pd.read_csv("/mnt/9C64E5F864E5D554/Documents and Settings/capekp00/Documents/data_statistika/Junior/Inkubace/IncubationExperiment/Python/pysourcefile/InitialConditions.csv")
#incond = pd.read_csv("/home/pcapek/IncubationExperiment/Python/pysource/InitialConditions.csv")
incond = np.array(incond)
#Import H+ ions derived from CO2 dissolution
HExtra = pd.read_csv("/mnt/9C64E5F864E5D554/Documents and Settings/capekp00/Documents/data_statistika/Junior/Inkubace/IncubationExperiment/Python/pysourcefile/H.csv")
#HExtra = pd.read_csv("/home/pcapek/IncubationExperiment/Python/pysource/H.csv")
HExtra = np.array(HExtra)
#model parameters
ParmsFermGroups = pd.DataFrame({'first': [6, 50, 50, 3, 1.6, 2, 0.5, 0.001, 0.2, 0.4],
                                'lower': [0.1, 1, 1, 0.5, 0.8, 0.05, 0.1, 1e-8, 1e-3, 1e-3],
                                'upper': [20, 500, 500, 10, 10, 3, 5, 1, 1, 1]})

#cost function
def cost(x, Y, incond, HExtra):
    #input to odeint
    p0 = [np.concatenate([x[0:8], x[0:8]])]
    #Conversion factors
    conversions = x[8:10]
    #intial fungal abundance
    a = incond[4]
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

bounds = Bounds(lb=np.array(ParmsFermGroups.loc[:, 'lower']),
                ub=np.array(ParmsFermGroups.loc[:, 'upper']))
result = differential_evolution(cost, bounds, args=(Y, incond, HExtra), updating='deferred', workers=2)

np.array(result.x).tofile("/mnt/9C64E5F864E5D554/Documents and Settings/capekp00/Documents/data_statistika/Junior/Inkubace/IncubationExperiment/Python/pyresults/optparsPL0.csv", sep = ";")
#np.array(result.x).tofile("/home/pcapek/IncubationExperiment/Python/pyresults/optpars.csv", sep = ";")