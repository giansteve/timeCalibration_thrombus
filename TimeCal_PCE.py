# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 14:17:59 2021

@author: pc
"""
#%% packages
# import uqpylab
from uqpylab import sessions
#%% preamble
myToken = 'b81b324bac9a3ebe0768e69efb5b98cd91b62279' # The user's token to access the UQCloud API
UQCloud_instance = 'https://beta.uq-cloud.io' # The UQCloud instance to use

# Start the session
#mySession = sessions.cloud(host=UQCloud_instance, token=myToken)
mySession = sessions.cloud() # since I have saved my credentials already
# (Optional) Get a convenient handle to the command line interface
uq = mySession.cli

# Reset the session
mySession.reset()

#%% code start
# Specify the options for a bivariate normal random vector

InputOpts = {
    'Marginals': [
        {
            'Type': 'Gaussian',
            'Parameters': [0,1]
        },
        {
            'Type': 'Gaussian',
            'Parameters': [0,1]
        }
    ]
}

# Create the bivariate normal random vector
myInput = uq.createInput(InputOpts)

# Draw 10 samples from the bivariate normal distribution
theSamples = uq.getSample(myInput, 10)