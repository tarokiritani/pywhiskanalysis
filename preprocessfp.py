# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 14:18:17 2016

@author: kiritani
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 23:19:19 2016

@author: kiritani
"""

import sqlite3
from collections import defaultdict
from pandas import Series as series
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
from os import chdir, makedirs
import pickle
import time
from functools import reduce
from itertools import groupby
from scipy.io import loadmat
from scipy.signal import correlate
from scipy.fftpack import fft, fftfreq
from scipy.interpolate import interp1d
from scipy.stats import wilcoxon
from pandas import concat, read_sql_query
import itertools
from scipy.ndimage.filters import maximum_filter, minimum_filter, median_filter
import pdb
from statistics import median_low, StatisticsError


def sql2neurons(exptype,celltype,snames):#coil stim 0p5Hz, free whisking, intrinsic current injection, intrinsic200ms
    c = sqlite3.connect('C:/Users/kiritani/Documents/data/GitHub/experiments/db/development.sqlite3')
    sql = "SELECT * FROM analyses INNER JOIN cells ON analyses.cell_id = cells.id INNER JOIN mice ON cells.mouse_id = mice.id"
    df = read_sql_query(sql,c)
    c.close()
    df = df[df['species_strain'].str.contains('|'.join(snames))]
    df = df[df['analysis_type'] == exptype]
    df = df[df['cell_type'] == celltype]
    neurons = []
    for ci in df['cell_id'].unique():
        neuron=df[df['cell_id']==ci]
        files=neuron['file'].values
        recs = [{'files':f} for f in files]
        neurons.append({'id':ci, 'recordings':recs, 'depth':neuron['depth'].values[0]})
    return neurons

def getfiles(neuron):
    for r in neuron['recordings']:
        mat={}
        for f in r['files'].split(';'):
            if f.endswith(('.mat','.xsg','.signal')):
                mat.update(loadmat(f,struct_as_record=False,squeeze_me=True))
        r['samplerate']=mat['header'].ephys.ephys.sampleRate
        r['ephys']=series(mat['data'].ephys.trace_1)
        r['ephys'].index += 1
        r['ephys'].index /= r['samplerate']
        dur=r['ephys'].index[-1]

        if mat['header'].stimulator.stimulator.pulseNameArray[1] == 'coil0p5Hz_2':
            r['stimtimes']=np.arange(1,dur,2)
        elif 'neg100' in mat['header'].ephys.ephys.pulseNameArray[0]:
            delay=mat['header'].ephys.ephys.pulseParameters[0].squarePulseTrainDelay
            width=mat['header'].ephys.ephys.pulseParameters[0].squarePulseTrainWidth
            r['stimtimes']=np.arange(delay,dur,width)
        
        r['whisker']=series(mat['angleArray'])
        r['whisker'].index += 1
        r['whisker'].index /= 500
        if 'peakTiming' in mat:
            r['spkpeaks']=np.reshape(mat['peakTiming'],(-1))
        if 'onOffTiming' in mat:
            r['onOff']=mat['onOffTiming']
            
    return neuron

def whiskclass(neuron, mindur=0.5, threshold=0.8):
    for i, r in enumerate(neuron['recordings']):
        whisking = (r['whisker'].diff().abs() > threshold).astype(int)
        idx = 0       
        for movement, g in groupby(whisking):
            tpoints = whisking.iloc[idx:idx + len(list(g))].index
            dur = tpoints[-1] - tpoints[0]
            if not(movement) and dur < 0.5:
                whisking[tpoints] = 1
            idx += len(tpoints)
        
        idx = 0
        for movement, g in groupby(whisking):
            tpoints = whisking.iloc[idx:idx + len(list(g))].index
            dur = tpoints[-1] - tpoints[0]
            if movement and dur < mindur:
                whisking[tpoints] = 0
            idx += len(tpoints)

        r['whisking'] = whisking
        r['whiskons'] = whisking[whisking.diff()==1].index
        r['whiskoffs'] = whisking[whisking.diff()==-1].index
        r['whisk'] = {'seg':[]}
        r['quiet'] = {'seg':[]}
        idx = 0
        whisking_high_res = r['whisking'].reindex(r['ephys'].index,method='nearest')
        for tf,g in groupby(whisking_high_res):
            dur = len(list(g))
            if tf:
                r['whisk']['seg'].append(r['ephys'].iloc[idx:idx+dur])
            else:
                r['quiet']['seg'].append(r['ephys'].iloc[idx:idx+dur])
            idx += dur
        
        fig, axes = plt.subplots(2,1)
        for p in zip(axes,['ephys','whisker'],['k','g']):
            p[0].plot(r[p[1]],p[2])
            p[0].fill_between(r['whisking'].index, *p[0].get_ylim(), r['whisking'], facecolor='red',alpha=.5)
        axes[0].set_title(neuron['depth'])
        if 'protractions' in r:
            axes[1].plot(r['protractions'],'+r')
            axes[1].plot(r['retractions'],'+g')
        if 'onOff' in r:
            for contact in r['onOff']:
                axes[0].axvspan(*contact,color='red')
        fname = 'cell_%s_%s' % (neuron['id'], i)
        fig.savefig(fname+'.png')
        fig.set_canvas(fig.canvas)
        pickle.dump(fig, open(fname+'.pickle', 'wb'))
        plt.close()
    neuron['quiet']=defaultdict(dict)
    neuron['whisk']=defaultdict(dict)
    return neuron

def onsetanalysis(neuron, exptype):
    keys=['t','pre','post','xsg','spk']
#    whiskon = ['whiskons',0.5,1,'whiskonxsgs','whiskonspks']
#    whiskoff = ['whiskoffs',0.5,1,'whiskoffxsgs','whiskoffspks']
    spkpeaks = ['spkpeaks',0.05,0.1,'spkvxsgs','spkspks']
    stims = ['stimtimes',0.1,0.5,'stimxsgs','stimspks']
    intrinsic = ['stimtimes',0.25,0.75,'stimxsgs','stimspks']
    #intrinsic200 = ['stimtimes',0.125,0.375,'stimxsgs','stimspks']
    rectypes = {'free whisking': [spkpeaks],'coil stim 0p5Hz': [spkpeaks, stims],'active touch': [stims], 'intrinsic current injection': [spkpeaks, intrinsic]}[exptype]
    analyses = list(map(lambda x: dict(zip(keys,x)), rectypes))
    for r in neuron['recordings']:
        for kk, a in enumerate(analyses):
            if a['t'] in r.keys():
                t0 = a['pre']
                t1 = a['post']
                whiskstates=r['whisking'].reindex(r[a['t']],method='nearest')
                if exptype=='intrinsic current injection':
                    for t in whiskstates.index:
                        wr = r['whisking'][(r['whisking'].index>t-0.1)*(r['whisking'].index<t+0.2)].values
                        if not all([x==wr[0] for x in wr]):
                            print('drop!')
                            whiskstates=whiskstates.drop(t)
            r[a['xsg']]=defaultdict(dict)
            r[a['spk']]=defaultdict(dict)
            for ii, (state, boo) in enumerate(zip(['whisk','quiet','all'],[[1],[0],[0,1]])):
                sr = r['samplerate']
                tvec=np.arange(-t0*sr,t1*sr)/sr
                r[a['xsg']][state]['seg']=list(map(lambda x:r['ephys'].reindex(tvec+x,method='nearest',tolerance=1/sr),whiskstates[whiskstates.isin(boo)].index))
                r[a['spk']][state]['seg']=list(map(lambda x:r['spkpeaks'][(r['spkpeaks']>x-t0)&(r['spkpeaks']<x+t1)]-x, whiskstates[whiskstates.isin(boo)].index))
                tracesegs=r[a['xsg']][state]['seg']
                if np.ndim(tracesegs)==2:
                    for fname, fcn in zip(['mean','median'], [np.nanmean, np.nanmedian]):
                        reducedtr=fcn(tracesegs,axis=0)
                        r[a['xsg']][state][fname]=series(reducedtr, tvec)
    return neuron

def calc_fft(neuron):
    for qw in ['quiet', 'whisk']:
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
        try:
            ffts=np.average(ffts,axis=0,weights=fft_lens)
            neuron[qw]['FFT1-5 (mV)']=np.trapz(ffts,dx=4/50)
        except ZeroDivisionError:
            neuron['quiet']['FFT1-5 (mV)']=np.nan
            neuron['whisk']['FFT1-5 (mV)']=np.nan
            break
    return neuron

def calc_ap_rates(n):
    for i, qw in enumerate(['quiet', 'whisk']):
        aps=[]
        seglens=[]
        for r in n['recordings']:
            aps.append(sum(r['whisking'].reindex(r['spkpeaks'],method='nearest')==i))
            seglens.append(sum([rs.index[-1] - rs.index[0] for rs in r[qw]['seg']]))
        n[qw]['AP rate (Hz)'] = sum(aps)/sum(seglens)
    return n

def calc_std_vm(n):
    for qw in ['quiet', 'whisk']:
        values = []
        for r in n['recordings']:
            for seg in r[qw]['seg']:
                values.extend(seg.values)
        n[qw]['std Vm (mV)'] = concat([concat(r[qw]['seg']) for r in n['recordings']]).std()
    return n

def calc_mean_vm(n):
    for qw in ['quiet', 'whisk']:
        values = []        
        for r in n['recordings']:
            for seg in r[qw]['seg']:
                values.extend(seg.values)
        n[qw]['mean Vm (mV)'] = np.mean(values)
    return n

def calc_crosscorr(n):
    for r in n['recordings']:
        xsgsegs=r['whisk']['seg']
        whisksegs = [r['whisker'].reindex(s.index,method='nearest') for s in xsgsegs]
    for segs in zip(xsgsegs, whisksegs):
        print(correlate(segs[0],segs[1],'full'))
    return n

def getwhiskpeaks(n, windowsize=30):
    protractions=[]
    retractions=[]
    for r in n['recordings']:
        tr = r['whisker']
        maxfiltered = maximum_filter(tr, size=windowsize/2) # 500 fps
        minfiltered = minimum_filter(tr, size=windowsize/2)
        propeaks = tr[(tr==maxfiltered) & (tr - minfiltered > 5)]
        propeaks_index = propeaks.index
        propeaks_index = np.split(propeaks_index, np.where(np.diff(propeaks_index)>0.01)[0]+1)
        try:
            r['protractions']=propeaks[map(median_low, propeaks_index)]
        except StatisticsError:
            r['protractions']=DataFrame()
        repeaks = tr[(tr==minfiltered) & (maxfiltered - tr > 5)]
        repeaks_index = repeaks.index
        repeaks_index = np.split(repeaks_index, np.where(np.diff(repeaks_index)>0.01)[0]+1)
        try:
            r['retractions']=repeaks[map(median_low, repeaks_index)]        
        except StatisticsError:
            r['retractions']=DataFrame()
        
#        protractions.append(cut_up(r['ephys'], r['protractions'].index, 0.1, 0.1, r['samplerate']))
#        retractions.append(cut_up(r['ephys'], r['retractions'].index, 0.1, 0.1, r['samplerate']))
#    n['ephysprotractions']=concat(protractions, axis=1)
#    n['ephysretractions']=concat(retractions, axis=1)
    return n

def hilbert_transform(n):
    return n

def intrinsics(n):
    dfall = []
    dfwhisk = []
    for r in n['recordings']:
        dfall.append(cut_up(r['ephys'], r['stimtimes'],0.2, 0.2, r['samplerate']))
        
        dfwhisk.append()
    
    return n

def whiskOnOff(n):
    dfons=[]
    dfoffs=[]
    for r in n['recordings']:
        dfons.append(cut_up(r['ephys'], r['whiskons'], 0.2,0.2, r['samplerate']))
        dfoffs.append(cut_up(r['ephys'], r['whiskoffs'], 0.2,0.2, r['samplerate']))
    n['whiskons']=concat(dfons, axis=1)
    n['whiskoffs']=concat(dfoffs, axis=1)
    return n

def touchOnOff(n):
    dfons = []
    dfoffs = []
    for r in n['recordings']:
        ons = r['onOff'][:,0]
        offs = r['onOff'][:,1]
        dfons.append(cut_up(r['ephys'], ons, 0.2,0.2, r['samplerate']))
        dfoffs.append(cut_up(r['ephys'], offs, 0.2,0.2, r['samplerate']))
    n['touchons']=concat(dfons, axis=1)
    n['touchoffs']=concat(dfoffs, axis=1)
    return n

def cluster_touches(n, gap=0.3):
    for r in n['recordings']:
        intervals = r['onOff'][1:,0]-r['onOff'][:-1,1]
        cut_here = np.where(intervals > gap)[0] + 1
        r['bouts'] = np.split(r['onOff'],cut_here)
    return n

def cut_up(trace, times, pre, post, samplerate):
    timevec = np.arange(-pre*samplerate, post*samplerate)/samplerate
    df=DataFrame()
    for t in times:
        ts = trace.reindex(timevec + t, limit = 1, method = 'pad')
        ts.index = timevec
        df[t]=ts
    return df
    
def remove_aps(n, rolling_size=2):
    """
    This created a trace without action potentials using median filter. The filter size in ms is given as an argument.
    """
    for r in n['recordings']:
        rolling_size_in_points = int(r['samplerate'] * rolling_size / 1000)
        no_ap = series(median_filter(r['ephys'], size=rolling_size_in_points))
        no_ap.index += 1
        no_ap.index /= r['samplerate']
        r['ephys_no_ap'] = no_ap
    return n



def main():
    exptype = 'free whisking'
    celltype = 'tdTomato'
    snames = ['SOM']
    dirpath='C:/Users/kiritani/Documents/data/analysis/'+exptype+'/'+celltype+ ''.join(snames)+time.strftime('%Y%m%d%H%M')
    makedirs(dirpath+'/neurons')
    chdir(dirpath)
    if exptype == 'free whisking':
        fs = [getfiles,remove_aps,getwhiskpeaks,whiskclass,calc_ap_rates, calc_mean_vm, calc_std_vm,#calc_fft,
              ]
    elif exptype == 'coil stim 0p5Hz':
        fs = [getfiles,whiskclass,lambda x: onsetanalysis(x,exptype),]
    elif exptype == 'active touch':
        fs = [getfiles,whiskclass,touchOnOff,cluster_touches]
    neurons = reduce(lambda a, x: map(x,a),fs, sql2neurons(exptype,celltype,snames))

    groupdata = defaultdict(dict)
    for n in neurons:
        layer=(n['depth']>550 and 'L5')or(n['depth']>432 and 'L4')or('L23')
        for i, r in enumerate(n['recordings']):
            if exptype == 'free whisking':
                spikefig, spikeaxes = plt.subplots(2,2)
                on_off_fig, on_off_axes = plt.subplots(3,2)
                peakfig, peakaxes = plt.subplots(3,2)
                plts = (
                [spikeaxes[0,0],'spike average',r['ephys'], r['spkpeaks'], .1, .1, r['samplerate']],
                [on_off_axes[0,0],'whisk onsets',r['ephys'], r['whiskons'], .1, .2, r['samplerate']],
                [on_off_axes[0,1],'whisk offsets',r['ephys'], r['whiskoffs'], .1, .2, r['samplerate']],
                [on_off_axes[2,0],'',r['whisker'], r['whiskons'], .1, .2, 500],
                [on_off_axes[2,1],'',r['whisker'], r['whiskoffs'], .1, .2, 500],
                [peakaxes[0,0],'protractions',r['ephys'], r['protractions'].index, .1, .1, r['samplerate']],
                [peakaxes[0,1],'retractions',r['ephys'], r['retractions'].index, .1, .1, r['samplerate']],
                [peakaxes[2,0],'',r['whisker'], r['protractions'].index, .1, .1, 500],
                [peakaxes[2,1],'',r['whisker'], r['retractions'].index, .1, .1, 500],
                )
                rasters = (
                    [on_off_axes[1,0],r['whiskons'],.1,.2],
                    [on_off_axes[1,1],r['whiskoffs'],.1,.2],
                    [peakaxes[1,0],r['protractions'].index,.1,.1],
                    [peakaxes[1,1],r['retractions'].index,.1,.1],
                           )
                figs = zip(['spike','whiskonoff','peaks'], [spikefig, on_off_fig, peakfig])
            elif exptype == 'intrinsic current':
                rn, cl = 4,4
                subaxes=(
                [(0,0), r['stimxsgs']['whisk'], 'whisking'],
                [(0,1), r['stimxsgs']['quiet'], 'quiet'])
            elif exptype == 'coil ':
                rn, cl = 4,4
                subaxes=(
                )
            elif exptype == 'active touch':
                on_off_fig, on_off_axes = plt.subplots(2,2)
                plts = (
                [on_off_axes[0,0],'touch onsets', r['ephys'], r['onOff'][:,0], 1, 1, r['samplerate']],
                [on_off_axes[0,1],'touch offsets', r['ephys'], r['onOff'][:,1], 1, 1, r['samplerate']]
                )
                rasters = (
                [on_off_axes[1,0],r['onOff'][:,0], 1, 1],
                [on_off_axes[1,1],r['onOff'][:,1], 1, 1],  
                )
                figs = zip(['touchonoff',], [ on_off_fig, ])
            for p in plts:
                df = cut_up(*p[2:])
                try:
                    p[0].plot(df.mean(1),'r', zorder=3)
                except:
                    print(df)
                try:
                    p[0].plot(df,'b')
                except:
                    print(df)
                p[0].set_xlim(-p[4],p[5])
                p[0].set_title(p[1])

            for raster in rasters:
                segs = list(map(lambda x:r['spkpeaks'][(r['spkpeaks']>x-raster[2])&(r['spkpeaks']<x+raster[3])]-x, raster[1]))
                for ii, trial in enumerate(segs):
                    raster[0].vlines(trial, ii+0.5, ii+1.5)
                raster[0].set_ylim(.5, ii+.5)
                raster[0].set_xlim(-raster[2],raster[3])

            for name, fig in figs:
                fname = 'cell_%s_%s' % (n['id'], i)
                fig.savefig(fname+name+'.png')
                pickle.dump(fig, open(fname+name+'.pickle', 'wb'))
                plt.close()

        for ii, analysis in enumerate(n['quiet'].keys()):
            groupdata[layer].setdefault(analysis,[]).append([n['quiet'][analysis], n['whisk'][analysis]])
        for analysis in set(n.keys())-{'recordings','quiet','whisk'}:
            groupdata[layer].setdefault(analysis,[]).append(n[analysis])
        pickle.dump(n, open(dirpath+'/neurons/'+fname+'.pickle','wb'))

    if exptype == 'free whisking':
        fig, axes = plt.subplots(4,4,figsize=(16, 12))
        for ii, analysis in enumerate(['AP rate (Hz)','mean Vm (mV)','std Vm (mV)','FFT1-5 (mV)']):
            for col, L in enumerate(['L23','L4']):
                gmat = np.array(groupdata[L][analysis])
                axes[ii, col].plot(gmat.transpose(),color='blue',alpha=0.5)
                axes[ii, col].errorbar([0,1], np.nanmean(gmat,axis=0),yerr=np.nanstd(gmat,axis=0), color='black')
                [T, p]=wilcoxon(gmat[:,0], gmat[:,1])
                axes[ii, col].text(0.5, axes[ii, col].get_ylim()[1], 'p=%s' % float('%.2g' %p),ha='center',va='top')
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
                axes[ii,2].axvline(432, color='gray')
                axes[ii,2].set_ylabel(title)
                axes[ii,3].scatter(groupdata[L]['depth'], gmat[:,0])
            axes[ii,0].set_ylabel(analysis)
        axes[ii,col+1].set_xlabel('$\mu$m')
        axes[1,2].plot()
        axes[0,0].set_title('=< 432 um')
        axes[0,1].set_title('> 432 um')
        fig.savefig('group.png')
        pickle.dump(fig, open('group.pickle', 'wb'))
        
        fig, axes = plt.subplots(4,4,figsize=(16, 12))
        
        
    elif exptype == 'active touch':
        pass
        
if __name__ == "__main__":
    main()