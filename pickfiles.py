# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 20:27:34 2016

@author: kiritani
"""
def pickfiles(filename):
    mouse = str(int(filename[2:6]))
    cell = filename[6:10]
    recording = filename[10:14]
    print(mouse)
    print(cell)
    print(recording)
    xsgpath = "C:\\Users\\kiritani\\Documents\\data\\cells\\TK"+mouse+"\\Recording\\"+cell+"\\"+filename
    spikepath = xsgpath.replace(".xsg","spike.mat")
    whiskpath = "C:\\Users\\kiritani\\Documents\\data\\analysis\\whiskerTracking\\TK"+mouse+"\\"+ cell + recording + "whisk.mat"
    activetouchpath = whiskpath.replace("whisk.","active.")
    paths = ";".join([xsgpath, spikepath, whiskpath, activetouchpath])
    print(paths)