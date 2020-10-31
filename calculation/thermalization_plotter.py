# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 09:18:39 2014

@author: dibakarsigdel
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 11:27:13 2014

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
x = [0.0 for k in range(100)]
y = [0.0 for k in range(100)]
ver = [x,y]
for k in range(3,8):
    st = str(k)
    wnr('su'+st+'t.dat',ver).reader()
    plt.figure(102)
    plt.grid()
    plt.scatter(ver[0],ver[1])



    plt.text(102,0.15,'SU(3)')
    plt.text(102,0.26,'SU(4)')
    plt.text(102,0.65,'SU(5)')
    plt.text(102,0.88, 'SU(6)')
    plt.text(102,0.93, 'SU(7)')







    plt.title('Thermalization')
    plt.xlabel('Iterations')
    plt.ylabel('Average plaquette')
    plt.show()




