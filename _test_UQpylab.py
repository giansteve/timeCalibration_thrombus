# -*- coding: utf-8 -*-
'''
Created on Thu Nov 11 14:17:59 2021

@author: pc
'''
#%% packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from uqpylab import sessions
#%% preamble
myToken = 'b81b324bac9a3ebe0768e69efb5b98cd91b62279' # The user's token to access the UQCloud API
UQCloud_instance = 'https://beta.uq-cloud.io' # The UQCloud instance to use

# Start the session
# mySession = sessions.cloud(host=UQCloud_instance, token=myToken)
mySession = sessions.cloud() # since I have saved my credentials already
# (Optional) Get a convenient handle to the command line interface
uq = mySession.cli

# Reset the session
mySession.reset()

#%% Create INPUT space
#%% from Matlab
# kbp = [8.000e-12, 8.000e-8]
# cbpt = [10, 250e4]
Ns = 2000;

InputOpts = {
    'Marginals': [
        {'Type': 'Uniform','Parameters': [1.00e-10, 1.00e-6]},      # Dc
        {'Type': 'Uniform','Parameters': [200e2, 200e4]},           # kc 
        {'Type': 'Uniform','Parameters': [8.000e-12, 8.000e-8]},    # kBP
        {'Type': 'Uniform','Parameters': [10.00e2, 10.00e4]},       # Ct
        {'Type': 'Uniform','Parameters': [10, 250e4]},              # cBPt
        {'Type': 'Uniform','Parameters': [1.00e-01, 3.00e+00]},     # TRt
        {'Type': 'Uniform','Parameters': [100, 10.00e4]},           # kc,wall
        {'Type': 'Uniform','Parameters': [1.5e14/20, 4.5e14/20]},   # cAP
        {'Type': 'Uniform','Parameters': [100, 250e3]},             # cBPbt
        {'Type': 'Uniform','Parameters': [0.1, 50]}                 # gamma_t
    ]
}

# Create the bivariate normal random vector
myInput = uq.createInput(InputOpts)
# Draw 10 samples from the bivariate normal distribution
# theSamples = uq.getSample(myInput, Ns)

#%% Read INPUT and OUTPUT from excel file
exp_design = pd.read_excel('fromMATLABtoPYTHON.xlsx', usecols='A:J', sheet_name='input', index_col=None, header=None)
OUT_hs = pd.read_excel('fromMATLABtoPYTHON.xlsx', usecols='A:BVL', sheet_name='OUT_HS', index_col=None, header=None)
OUT_ls = pd.read_excel('fromMATLABtoPYTHON.xlsx', usecols='A:BVL', sheet_name='OUT_LS', index_col=None, header=None)
OUT_time = pd.read_excel('fromMATLABtoPYTHON.xlsx', usecols='A', sheet_name='OUT_time', index_col=None, header=None)

# transform DataFrames into lists
ed_list = exp_design.values.tolist()
outputHS = OUT_hs.values.tolist()
outputLS = OUT_ls.values.tolist()

#%% test
NED = [50, 100, 150, 200, 500]
for n in enumerate(NED):
    # print(i)
    print(n)

#%% Create metamodel of the thrombus model
MRI_time_index = list(range(0,70,10))

for i,n in enumerate(MRI_time_index):
    # print(i)
    # print(n)

    MetaOpts= {
        'Type': 'Metamodel',
        'MetaType': 'PCE',
        'Method': 'LARS',
        'Display': 'verbose',
        'Degree': np.arange(8,10).tolist(),
        'DegreeEarlyStop': False,
    }
    
    MetaOpts['TruncOptions'] = {'qNorm': 0.95}
    MetaOpts['ExpDesign'] = {
        'X': ed_list,
        'Y': outputHS[n]
        }
    
    myPCE_LARS = uq.createModel(MetaOpts)
    uq.print(myPCE_LARS)

# code from Matlab
# metamodel.Type = 'Metamodel';
# metamodel.MetaType = 'PCE';
# metamodel.Display = 'verbose';
# metamodel.Method = 'lars';
# metamodel.Degree = 3:15;
# metamodel.TruncOptions.qNorm = 0.95;
# metamodel.DegreeEarlyStop = false;
# metamodel.Input = INPUT;
# metamodel.ExpDesign.NSamples = Ns;
# metamodel.ExpDesign.X = exp_design;
# % PCE H/S
# metamodel.ExpDesign.Y = OUT_HS(MRI_time_index,:)';
# PCE_HS = uq_createModel(metamodel);
# % PCE L/S
# metamodel.ExpDesign.Y = OUT_LS(MRI_time_index,:)';
# PCE_LS = uq_createModel(metamodel);