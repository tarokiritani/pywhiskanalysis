# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 15:48:28 2016

@author: kiritani
"""
from preprocessfp import *
from os import chdir, makedirs
import pickle
import time
import numpy as np
import matplotlib.pyplot as plt
import pdb

def main():
    celltype = 'tdTomato'
    snames = ['SOM']
    dirpath='C:/Users/kiritani/Documents/data/analysis/freeVsActive/'+celltype+ ''.join(snames)+time.strftime('%Y%m%d%H%M')
    makedirs(dirpath+'/neurons')
    chdir(dirpath)
    freeneurons = sql2neurons('free whisking', celltype, snames)
    activeneurons = sql2neurons('active touch', celltype, snames)
    neuron_ids = {n['id'] for n in freeneurons} & {n['id'] for n in activeneurons}
    pdb.set_trace()
    neuron_ids = list(neuron_ids)
    for neuron_id in neuron_ids:
        print(neuron_id)
        fneuron = next(fn for fn in freeneurons if fn['id']==neuron_id)
        aneuron = next(an for an in activeneurons if an['id']==neuron_id)
        fneuron = getfiles(fneuron)
        aneuron = getfiles(aneuron)
        aneuron = touchOnOff(aneuron)
        fneuron = whiskclass(fneuron)
        f, axes = plt.subplots(2,2)
        for r in aneuron['recordings']:
            df = cut_up(r['ephys'], r['onOff'][:,0],0.1,0.1,r['samplerate'])
            axes[0,0].plot(df.mean(1),'r',zorder=3)
            axes[0,0].plot(df,'b')
            dfwhisk = cut_up(r['whisker'], r['onOff'][:,0],0.1,0.1,500)
            axes[1,0].plot(dfwhisk.mean(1),'r',zorder=3)
            axes[1,0].plot(dfwhisk,'b')
            threshold = r['whisker'].reindex(r['onOff'][:,0],method='nearest').mean()
        for r in fneuron['recordings']:
            crossing = r['whisker'] > threshold
            crossing = crossing.values
            crossing[1:][crossing[:-1] & crossing[1:]] = False
            d = np.gradient(r['whisker']) > 0
            crosses = r['whisker'][crossing & d]
            df = cut_up(r['ephys'],crosses.index,0.1,0.1,r['samplerate'])
            axes[0,1].plot(df.mean(1),'r', zorder=3)
            axes[0,1].plot(df,'b')
            dfwhisk = cut_up(r['whisker'], crosses.index,0.1,0.1,500)
            axes[1,1].plot(dfwhisk.mean(1),'r',zorder=3)
            axes[1,1].plot(dfwhisk,'b')
        axes[0,0].set_title('active touch')
        axes[0,1].set_title('free whisking')
        f.show()
        fname = 'cell_%s' % neuron_id
        f.savefig(fname+'.png')
        f.set_canvas(f.canvas)
        pickle.dump(f, open(fname+'.pickle', 'wb'))

if __name__ == "__main__":
    main()