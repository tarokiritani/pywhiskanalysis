# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:31:08 2016

@author: kiritani
"""

from collections import defaultdict
from pandas import concat
import numpy as np
import matplotlib.pyplot as plt
import pickle
import scipy.ndimage as ni
from tkinter.filedialog import askopenfilename, asksaveasfilename
from scipy.fftpack import fft, fftfreq
from scipy.interpolate import interp1d
from scipy.stats import wilcoxon
from functools import reduce
from os.path import join, dirname
import itertools
from itertools import groupby, zip_longest
from scipy.signal import correlate

#neurons=pickle.load(open(fn,'rb'))

def calc_fft(neuron):
    for qw in ['quiet', 'whisk']:
        neuron[qw]=defaultdict(dict)
        ffts=[]
        fft_lens = []
        for r in neuron['recordings']:
            for s in r[qw]['seg']:
                if (s.index[-1]-s.index[0]) > 1:
                    yf = fft(s.values)
                    xf = fftfreq(s.size,s.index[1]-s.index[0])
                    fftint = interp1d(xf[0:(len(xf)+1)/2],2/s.size*np.abs(yf[0:(len(xf)+1)/2]))
                    ffts.append(fftint(np.linspace(1.0001,5,num=50)))
                    fft_lens.append(s.index[-1]-s.index[0])
                    print(len(fft_lens))
        ffts=np.average(ffts,axis=0,weights=fft_lens)
        neuron[qw]['FFT1-5 (mV)']=np.trapz(ffts,dx=4/50)
    return neuron
    
def calc_ap_rates(n):
    for i, qw in enumerate(['quiet', 'whisk']):
        n[qw]=defaultdict(dict)
        aps=[]
        seglens=[]
        for r in n['recordings']:
            aps.append(sum(r['whisking'].reindex(r['spkpeaks'],method='nearest')==i))
            seglens.append(sum([rs.index[-1] - rs.index[0] for rs in r[qw]['seg']]))
        n[qw]['AP rate (Hz)'] = sum(aps)/sum(seglens)
    return n
    
def calc_std_vm(n):
    for i, qw in enumerate(['quiet', 'whisk']):
        n[qw]['std Vm (mV)'] = concat([concat(r[qw]['seg']) for r in n['recordings']]).std()
    return n
    
def calc_mean_vm(n):
    for qw in ['quiet', 'whisk']:
        n[qw]['mean Vm (mV)'] = concat([concat(r[qw]['seg']) for r in n['recordings']]).mean()
    return n

def calc_onset(n):
    segs=[r['whiskonxsgs']['all']['seg'] for r in n['recordings']]
    segs=list(itertools.chain.from_iterable(segs))
    segs=list(zip_longest(*segs, fillvalue=np.nan))
    if np.ndim(segs)==2:
        segs=np.nanmean(segs,axis=1)
    n['onset_mean']=segs
    return n

def calc_crosscorr(n):
    for r in n['recordings']:
        xsgsegs=r['whisk']['seg']
        whisksegs = [r['whisktrace'].reindex(s.index,method='nearest') for s in xsgsegs]
    for segs in zip(xsgsegs, whisksegs):
        print(correlate(segs[0],segs[1],'full'))
        


    
def freewhiskanalysis(neurons):
    neurons=neurons.values()
    neurons = reduce(lambda a, x: map(x,a),[calc_ap_rates, calc_mean_vm, calc_std_vm,calc_onset], neurons)
    
    groupdata = defaultdict(dict)
    analyses = ['mean Vm (mV)','std Vm (mV)','AP rate (Hz)']
    for n in neurons:
        #analyses = n['quiet'].keys()
        print(analyses)
        for ii, analysis in enumerate(analyses):
            if n['depth']>550:
                layer = 'L5'
            elif n['depth']>=450:
                layer = 'L4'
            else:
                layer = 'L23'
            groupdata[layer].setdefault(analysis,[]).append([n['quiet'][analysis], n['whisk'][analysis]])
        groupdata[layer].setdefault('depth',[]).append(n['depth'])
        groupdata[layer].setdefault('onsets',[]).append(n['onset_mean'])
        
    #figure()
    #for onset in groupdata['L23']['onsets']:
    #    plt.plot(onset)
    
    fig, axes = plt.subplots(len(analyses),3,figsize=(16, 12))
    for ii, analysis in enumerate(analyses):
        for col, L in enumerate(['L23','L4']):
            gmat = np.array(groupdata[L][analysis])
            axes[ii, col].plot(gmat.transpose(),color='blue',alpha=0.5)        
            axes[ii, col].errorbar([0,1], np.nanmean(gmat,axis=0),yerr=np.nanstd(gmat,axis=0), color='black')
            [T, p]=wilcoxon(gmat[:,0], gmat[:,1])
            axes[ii, col].text(0.5, axes[ii, col].get_ylim()[1], 'p={:0.2f}'.format(p),ha='center',va='top')
            axes[ii, col].set_xticks([0,1])
            axes[ii, col].set_xticklabels(['quiet','whisk'])
            axes[ii, col].set_xlim([-.3, 1.3])
            if analysis == 'mean Vm (mV)':
                gdiff = gmat[:,1]-gmat[:,0]
                title='$\Delta$ Vm'
                y=0
            else:
                gdiff = gmat[:,1]/gmat[:,0]
                title='whisk/quiet'
                y=1
            axes[ii,2].scatter(groupdata[L]['depth'], gdiff)
            axes[ii,2].axhline(y, color='gray')
            axes[ii,2].set_ylabel(title)
        axes[ii,0].set_ylabel(analysis)
    axes[ii,col+1].set_xlabel('$\mu$m')
    axes[1,2].plot()
    axes[0,0].set_title('=< 450 um')
    axes[0,1].set_title('> 450 um')
    fn=askopenfilename()
    fig.savefig(join(dirname(fn),'groupanalysis.png'))
    pickle.dump(neurons,open(asksaveasfilename(defaultextension=".pickle", initialdir='C:/Users/kiritani/Desktop'),'wb'))