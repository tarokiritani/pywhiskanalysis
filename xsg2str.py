# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 12:57:56 2016

@author: kiritani
"""
import sqlite3
from scipy.io import loadmat
import os
import re
import pdb

dfolder = 'C:/Users/kiritani/Documents/data/'

def xsg2str(xsgpath):
    mat = loadmat(xsgpath,struct_as_record=False,squeeze_me=True)
    trace = mat['data'].ephys.trace_1[::4]
    #global dfolder    
    outputpath = dfolder + 'amazons3/'+os.path.basename(xsgpath)
    outputpath = outputpath.replace('.xsg','.txt')
    with open(outputpath,'w') as text_file:
        text_file.write("sec,vm\\n")
        for i, t in enumerate(trace):
            text_file.write("%.4f,%.2f\\n" % (i/10000, t))

def mat2str(matpaths, xsgpath):
    for path in matpaths:
        mat = loadmat(path, struct_as_record=False,squeeze_me=True)
        outputpath = dfolder + 'amazons3/'+os.path.basename(xsgpath)
        if 'angleArray' in mat:
            trace=mat['angleArray']
            with open(outputpath.replace('.xsg','.whisk'),'w') as text_file:
                text_file.write("sec,degree\\n")
                for i, t in enumerate(trace):
                    text_file.write("%.3f,%.1f\\n" % (i/500, t))
        if 'onOffTiming' in mat:
            onoffs = mat['onOffTiming']
            if onoffs.ndim == 1:
                onoffs = onoffs.reshape(-1,2)
            ons = ','.join(str(e) for e in onoffs[:,0]) + '\n'
            offs = ','.join(str(e) for e in onoffs[:,1])
            with open(outputpath.replace('.xsg','.contact'),'w') as text_file:
                text_file.write(ons)
                text_file.write(offs)

def main():
    tables = ['cells','analyses','mice']
    conn = sqlite3.connect(dfolder+'GitHub/experiments/db/development.sqlite3')    
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    conn_new = sqlite3.connect(dfolder +'GitHub/barrel_vm/db/development.sqlite3')
    conn_new.row_factory = sqlite3.Row
    c_new = conn_new.cursor()
    for table in tables:
        c.execute('SELECT * FROM ' + table)
        rows = c.fetchall()
        for row in rows:
            keys = row.keys()
            try:
                keys.remove('user_id')
            except ValueError:
                pass
            columns = ','.join(keys)    
            qmarks = '?,'*len(keys)
            qmarks = qmarks[0:-1]
            c_new.execute("INSERT OR REPLACE INTO " + table + " ("+columns+") VALUES ("+qmarks+");",tuple(map(lambda x: row[x], keys)))
            print(row['id'])
            conn_new.commit()
    conn.close()
    c_new.execute("SELECT * FROM analyses")
    rows = c_new.fetchall()

    for row in rows:
        try:
            xsgpath = re.search('\w[^;]*.xsg',row['file'])
            matpaths = re.findall('\w[^;]*.mat',row['file'])
            
            if xsgpath:
                xsgpath = xsgpath.group(0)
                xsg2str(xsgpath)
                try:
                    mat2str(matpaths,xsgpath)
                except FileNotFoundError:
                    pass
                file = os.path.basename(xsgpath).replace('.xsg','')
                print(file)
                c_new.execute("UPDATE analyses SET file = '" + file + "' WHERE id = %s" % row['id'])
            else:
                c_new.execute("DELETE FROM analyses WHERE id = %s" % row['id'])
            conn_new.commit()
        except AttributeError:
            pass
    conn_new.close()
    
if __name__ == "__main__":
    main()