# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 14:30:02 2014

@author: dibakarsigdel
"""



import numpy as np
import matplotlib.pyplot as plt






class wnr(object):
    
    def __init__ (self,filename,ver):
                self.filename = filename
                self.ver = ver
        

    def writer(self):
                 DataOut = np.column_stack(self.ver)    
                 np.savetxt(self.filename,DataOut)
                 return


    def reader(self):
                colno = len(self.ver)
                for k in range (colno):
                    self.ver[k]  = np.loadtxt(self.filename, unpack = True, usecols = [k])
                return self.ver











#===========================================================================
x = [0.0 for k in range(30)]
y = [0.0 for k in range(30)]
er = [0.0 for k in range(30)]
ver = [x,y,er]
wnr('su201.dat',ver).reader()
plt.figure(102)
plt.grid()
plt.xlabel('beta')
plt.ylabel('average plaquette')
plt.errorbar(ver[0],ver[1],yerr = ver[2],fmt = '8')
plt.show()




