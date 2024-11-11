import numpy as np
import pandas as pd
from scipy.integrate import odeint

def subMSolver0 (model, pars, t, y0):
    #model parameters
    ##Im, Km, yA, Em, g, k
    pars_model=pars[0:6]
    #conversion factors
    ##nG, nB
    conversions=pars[6:8]
    
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #labelled biomass C (Bl) including everything death
    Bl  = y[:, 2] + (y[:, 2] + y[:, 5])*y[:, 1]
    #kec defines the extractable fraction - all metabolites can be extracted
    kec = (y[:, 2]*conversions[1] + (y[:, 2] + y[:, 5])*y[:, 1]*conversions[0])/Bl
    #labelled chloroform flush
    Cflush = Bl*kec
    #kPLFA defines conversion factor between labelled PLFA and Bl
    #kPLFA = y[:, 2]*conversions[2]/Bl
    #Total PLFA
    #PLFA = Bl*kPLFA
                        
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(len(t), 1), #Glucose
                           y[:, 3].reshape(len(t), 1),#labelled CO2
                           Cflush.reshape(len(t), 1)), 
                           axis = 1)
    
    return yhat
