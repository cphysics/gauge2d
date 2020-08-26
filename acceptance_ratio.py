# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 12:41:56 2014

@author: dibakarsigdel
"""



import matplotlib.pyplot as plt
import random as random
import numpy as np
import math as math




class Normalized_Crutz_Distribution(object):
    
        def  __init__ (self,beta):
                  self.beta = beta
    
        def funct(self,t):
                   ft = (math.sqrt(1 - (t)**2))*(math.exp(self.beta*(t))) 
                   return ft
                   
        def nrn(self,npt):
                dx = 2.0/(2.0*float(npt))
                I = [0.0 for k in range(-npt,npt+1)]
                x = [k*dx for k in range(-npt,npt+1)]
                for k in range(-npt,npt):
                    I[k] = self.funct(x[k])*dx
                II = sum(I)
                return II           
                   
        def normalized_funct(self,t):   
              nfunct = self.funct(t)/self.nrn(50)
              return nfunct
              
    
        def plot_function(self):
                    x = np.arange(-1.0,1.0,0.001)
                    ct = len(x)
                    y = []
                    for k in range(ct):
                        tt = self.normalized_funct(x[k])
                        y.append(tt)
                    return x,y 
                    
                    
                    
                    
class Bounded_Crutz_Distribution(object): 
    
        def  __init__ (self,beta):
                  self.beta = beta      
        
        def find_max(self,ptn,fpt,dx):
                    dt = dx/float(ptn)
                    x = [0.0 for k in range(ptn)]
                    fx = [0.0 for k in range(ptn)]
                    x = [(fpt + dt*k) for k in range(ptn)]
                    NCD = Normalized_Crutz_Distribution(self.beta)
                    for k in range(ptn):
                        fx[k] = [NCD.normalized_funct(x[k])]
                    fmax = max(fx)
                    mx = (fpt + (fx.index(fmax))*dt)
                    xo = mx - dt
                    return xo,fmax
            
            
               
        def optimize(self,slice):
                    ftp = -1.0
                    dx= [0.0 for k in range (slice)]
                    for k in range(slice):
                        dx[k]  = 2.0/float(10**(k))
                    for k in range(slice):
                            ftp,fmax = self.find_max(20,ftp,dx[k])
                    return  fmax           
         
                 
         
        def bounded_funct(self,a0):
                    den = self.optimize(10)
                    NCD = Normalized_Crutz_Distribution(self.beta)
                    ft = NCD.normalized_funct(a0)/den[0]
                    return ft 
                    
                    
                    
                    
class Check_Crutz_Distribution(object):   
      
        def  __init__ (self,beta):
                  self.beta = beta        
                    
        def  crutz_generation(self):
                   count = 0
                   failed = 0
                   BCD = Bounded_Crutz_Distribution(self.beta)
                   while (count <1):
                       r = [random.uniform(0.0,1.0),\
                           random.uniform(0.0,1.0)]
                      
                       a0 = (2*r[0] - 1)
                       ff = BCD.bounded_funct(a0)
                       b0 = r[1] 
                       if (b0 < ff):
                         count = count +1 
                         return a0,failed
                       else:
                        count = 0
                        failed =  failed+1
                           
       
        def collect_x(self,LN): 
                    ll = 0
                    x = []
                    while ll < LN+1:
                        a0 = self.crutz_generation()
                        x.append(a0)
                        ll= ll+1
                    return x
                    
                    

def subplot_plotter(LN,SN):
                    nt = 1
                    while nt < 4+1:
                        beta = 1.0*nt
                        NCD = Normalized_Crutz_Distribution(beta)
                        CCD = Check_Crutz_Distribution(beta)
                        xx,yy = NCD.plot_function()
                        x = CCD.collect_x(LN)    
                        num_bins = SN
                        st = 'beta = '+ str(beta)
                        plt.figure(22)
                        plt.subplot(2,2,nt)
                        plt.scatter(xx,yy)
                        plt.xlabel(' eigen value eta at:'+ st)
                        plt.ylabel('normalized probability' )
                        plt.hist(x,num_bins, normed= 1.0, facecolor='green',alpha = 0.5)
                        nt = nt+1
           
           
#subplot_plotter(5000,50)
beta = 16.0
ll = 0
sfailed = 0
ttln = 2000
while ll < ttln :         
     a0,failed = Check_Crutz_Distribution(beta).crutz_generation()
     sfailed = sfailed+failed
     print a0
     ll = ll+1
print 'beta = ',beta,',  ','acceptance ratio = ',ttln/float(ttln+sfailed)            
          
          
          
          
          
          
          
          
          
          
          
          
          





