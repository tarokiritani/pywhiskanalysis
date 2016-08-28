# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 16:07:30 2016

@author: kiritani
"""
from pandas import Series as series
from scipy.io import loadmat
import preprocessfp
import numpy as np
from itertools import groupby
from scipy.signal import hilbert
f='C:/Users/kiritani/Documents/data/analysis/whiskerTracking/TK435/AAAA0001whisk.mat'

mat = loadmat(f,struct_as_record=False,squeeze_me=True)
a=mat['angleArray']
tseries=series(a,(np.arange(len(a))+1)/500)
threshold=0.8
whisking=(tseries.diff().abs()>threshold).astype(int)
ons= whisking.loc[whisking.diff()==1].index
offs = whisking.loc[whisking.diff()==-1].index
whiskseg = []
for off in offs:
    if ons.max() > off:
        on = ons[ons>off][0]
        if (on-off) < 0.5:
            whisking.loc[off:on]=1
bl = whisking.reindex(tseries.index,method='nearest').values
idx = 0
for tf,g in groupby(bl):
    dur = len(list(g))
    if tf:
        whiskseg.append(tseries.iloc[idx:idx+dur-1])
    idx = idx + dur
    
for seg in whiskseg:
    hil=hilbert(seg-seg.mean())
    re = np.real(hil)
    img = np.imag(hil)
    phase = np.arctan2(img,re)
    figure()
    plot(phase)
    amp = (re**2+img**2)**0.5
    figure()
    plot(amp)
    figure()    
    plot(amp*sin(arctan2(re,img)))
    