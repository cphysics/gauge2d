import numpy as np


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
#x = [0.0 for k in range(100)]
#y = [0.0 for k in range(100)]
#for k in range(100):
#    x[k] = random.random()
#    y[k] = random.random()
#wnr('data.dat',[x,y]).writer()
#---------------------------------------------------------------------------
#[x,y] = [[0.0 for k in range (100)]for l in range (2)]
#[xx,yy] = wnr('data.dat',[x,y]).reader()
#=====================================================================
