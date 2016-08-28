# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 14:27:52 2016

@author: kiritani
"""

from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import pickle
from tkinter.filedialog import askopenfilename
from itertools import zip_longest
from pandas import Series as series

neurons=pickle.load(open(askopenfilename(),'rb'))

for __, n in neurons.items():
    for qw in ['quiet', 'whisk']:
        segs=[]
        for i, r in enumerate(n['recordings']):
            if 'seg' in r['stimxsgs'][qw].keys():
                segs=segs+r['stimxsgs'][qw]['seg']
        segs = list(map(lambda x: x.values, segs))
        if len(segs)>3:
            med =np.nanmedian(list(zip_longest(segs, fillvalue=np.nan)),axis=0)
            sr=r['samplerate']
            med = series(med[0],(np.arange(len(med[0]))-1)/sr-0.25)
            pre=med[(med.index>-0.1)*(med.index<0)].median()
            post=med[(med.index>0.1)*(med.index<0.2)].median()
            n[qw]={'Rm (Mohms)':(pre-post)/100*1000} #in Mega ohm
        else:
            n[qw]={'Rm (Mohms)':np.nan}
analyses = n[qw].keys()
fig = plt.figure()
nrows=len(analyses)
ncols=3
groupdata = defaultdict(dict)
for ii, analysis in enumerate(analyses):
    for __, n in neurons.items():
        layers=[(550,'L5'),(450, 'L4'),(0,'L23')]
        for l in layers:
            if n['depth'] > l[0]:
                layer = l[1]
                break

        if not(np.isnan(n['quiet'][analysis])) and not(np.isnan(n['whisk'][analysis])):
            groupdata[layer].setdefault(analysis,[]).append([n['quiet'][analysis], n['whisk'][analysis]])
    for col, L in enumerate(['L23','L4']):
        ax = plt.subplot(nrows, ncols, nrows*ii+col+1)        
        ax.plot(np.transpose(groupdata[L][analysis]),color='blue',alpha=0.5)
        ax.errorbar([0,1], np.mean(groupdata[L][analysis],axis=0),yerr=np.std(groupdata[L][analysis],axis=0), color='black')
        ax.set_xticks([0,1])
        ax.set_xticklabels(['quiet','whisk'])
        ax.set_xlim([-.3, 1.3])
    plt.subplot(nrows, ncols, nrows*ii+1).set_ylabel(analysis)

plt.subplot(nrows,ncols,1).set_title('=< 450 um')
plt.subplot(nrows, ncols,2).set_title('> 450 um')
#pickle.dump(neurons,open(asksaveasfilename(defaultextension=".pickle", initialdir='C:/User/kiritani/Desktop'),'wb'))