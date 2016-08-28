# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 10:53:12 2016

@author: kiritani
"""

from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import pickle
from tkinter.filedialog import askopenfilename
from itertools import zip_longest
from pandas import Series as series
from os.path import join, dirname

filename = askopenfilename()
neurons=pickle.load(open(filename,'rb'))

for k, n in neurons.items():
    for qw in ['quiet', 'whisk','all']:
        segs=[]
        for i, r in enumerate(n['recordings']):
            if 'seg' in r['stimxsgs'][qw].keys():
                segs=segs+r['stimxsgs'][qw]['seg']
        segs = list(map(lambda x: x.values, segs))
        fig = plt.figure()
        if len(segs)>3:
            med =np.nanmedian(list(zip_longest(segs, fillvalue=np.nan)),axis=0)
            sr=r['samplerate']
            n[qw]['med'] = series(med[0],(np.arange(len(med[0]))-1)/sr-0.25)
            ax = plt.subplot(1, 3, i)
            ax.plot(med)
            
        else:
            n[qw]={'Rm (Mohms)':np.nan}
        fig.savefig(join(dirname(filename),'cell_'+str(k)+'.png'))
analyses = n[qw].keys()
fig = plt.figure()
nrows=len(analyses)
ncols=3
groupdata = defaultdict(dict)