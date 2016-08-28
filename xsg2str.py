# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 12:57:56 2016

@author: kiritani
"""
import sqlite3
from scipy.io import loadmat
import numpy as np
import os
import re
import pdb
def xsg2str(xsgfile):
    mat = loadmat(xsgfile,struct_as_record=False,squeeze_me=True)
    trace = mat['data'].ephys.trace_1[::4]
    sec = np.arange(10000*60)/10000
    outputfile = 'C:/Users/kiritani/Documents/data/amazons3/'+os.path.basename(xsgfile)
    outputfile = outputfile.replace('.xsg','.txt')
    with open(outputfile,'w') as text_file:
        text_file.write("sec,vm\\n")
        for s, t in zip(sec,trace):
            text_file.write("%.4f,%.2f\\n" % (s, t))
    return outputfile

def mat2str(matpaths,xsgfilename):
    for path in matpaths:
        mat = loadmat(path, struct_as_record=False,squeeze_me=True)
        if 'angleArray' in mat:
            trace=mat['angleArray']
            sec=np.arange(500*60)/10000
            outputfile = 'C/Users/kiritani/Documents/data/amazons3/whisk'+xsgfilename
            with open(outputfile,'w') as text_file:
                text_file.write("sec,degree\\n")
                for s, t in zip(sec,trace):
       	           text_file.write("%.3f,%.1f\\n" % (s, t))
    return outputfile

def main():
    tables = ['cells','analyses','mice']
    conn = sqlite3.connect('C:/Users/kiritani/Documents/data/GitHub/experiments/db/development.sqlite3')
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    conn_new = sqlite3.connect('C:/Users/kiritani/Documents/data/GitHub/barrel_vm/db/development.sqlite3')
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
                xsgfilename = xsg2str(xsgpath)
                mat2str(matpaths,xsgfilename)
                file = os.path.basename(xsgpath).replace('.xsg','.txt')
                c_new.execute("UPDATE analyses SET file = '" + file + "' WHERE id = %s" % row['id'])
            else:
                c_new.execute("DELETE FROM analyses WHERE id = %s" % row['id'])
            conn_new.commit()
            
        except AttributeError:
            pass
    conn_new.close()
    
if __name__ == "__main__":
    main()