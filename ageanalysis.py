# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 14:21:02 2016

@author: kiritani
"""

import os
from pickle import load
import matplotlib.pyplot as plt

os.chdir("C:\\Users\\kiritani\\Documents\\data\\analysis\\free whisking\\tdTomatoSOM201609131454\\neurons")

plt.figure(figsize=[9,4.5])
for file in os.listdir("C:\\Users\\kiritani\\Documents\\data\\analysis\\free whisking\\tdTomatoSOM201609131454\\neurons"):
    if file.endswith(".pickle"):
        neuron = load(open(file,'rb'))
        if 100 < neuron['depth'] < 300:
            delta = neuron['whisk']['mean Vm (mV)'] - neuron['quiet']['mean Vm (mV)']
            plt.plot(neuron['age'], delta,'o',color='orange')
plt.xlabel('postnatal day')
plt.ylabel('$\Delta$ Vm (whisk - quiet)')