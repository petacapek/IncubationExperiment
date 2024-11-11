import numpy as np
import pandas as pd
from scipy.integrate import odeint

def subMSolverFerm (model, pars, t, y0):
    #model parameters
    ##Im, Km, yA, Gm, k, emax
    pars_model=pars[0:6]
    #conversion factors
    ##ng, nb
    conversions=pars[6:8]
    
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate labelled biomass (Bl), kec factor, and 14CFC Flush
    Bl  = y[:, 2] + (y[:, 2] + y[:, 5])*y[:, 1]
    kec = (conversions[1]*y[:, 2] + conversions[0]*(y[:, 2] + y[:, 5])*y[:, 1])/Bl
    Flush = Bl*kec
                    
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(len(t), 1), #Glucose
                           y[:, 3].reshape(len(t), 1),#labelled CO2
                           #kec.reshape(len(t), 1), 
                           Flush.reshape(len(t), 1),
                           y[:, 6].reshape(len(t), 1)), 
                           axis = 1)
    
    return yhat
