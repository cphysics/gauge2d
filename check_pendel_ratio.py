# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 13:19:49 2014

@author: dibakarsigdel
"""


import matplotlib.pyplot as plt
import random as random
import numpy as np
import math as math




class Check_Pendelton_Distribution(object):
    
        def  __init__ (self,beta):
                  self.beta = beta
    
        def funct(self,t):
                   ft = (math.sqrt(1 - (t)**2))*(math.exp(self.beta*(t))) 
                   return ft
                   
        def nrn(self):
                dx = 0.002
                I = [0.0 for k in range(-500,501)]
                x = [k*dx for k in range(-500,501)]
                for k in range(-500,500):
                    I[k] = self.funct(x[k])*dx
                II = sum(I)
                return II           
                   
        def normalized_funct(self,t):   
              nfunct = self.funct(t)/self.nrn()
              return nfunct
             
        def plot_function(self):
                    x = np.arange(-1.0,1.0,0.001)
                    ct = len(x)
                    y = []
                    for k in range(ct):
                        tt = self.normalized_funct(x[k])
                        y.append(tt)
                    return x,y 
                    
        def pendelton_generation(self):
                    count = 0
                    failed = 0
                    while (count < 1):
                         r = [ random.uniform(0.0,1.0),\
                         random.uniform(0.0,1.0),\
                         random.uniform(0.0,1.0),\
                         random.uniform(0.0,1.0)]
                         dtm = 1.0
                         x  = [-(math.log(r[0])/(dtm*self.beta)),\
                               -(math.log(r[1])/(dtm*self.beta))]
                 
                         C = (math.cos(2.0*math.pi*r[2]))**2
                         A = x[0]*C
                         delta = x[1]+A
           
                         if (r[3]**2) < (1-(0.5*delta)):
                             a0 = (1- delta)
                             count = count+1
                             return a0,failed
                         else: 
                            count = 0
                            failed = failed+1
       
        def collect_x(self,LN): 
                    ll = 0
                    x = []
                    while ll < LN+1:
                        a0 = self.pendelton_generation()
                        x.append(a0)
                        ll= ll+1
                    return x
                    
                    

def subplot_plotter(LN,SN):
                    nt = 1
                    while nt < 4+1:
                        beta = 1.0*nt
                        xx,yy = Check_Pendelton_Distribution(beta).plot_function()
                        x = Check_Pendelton_Distribution(beta).collect_x(LN)    
                        num_bins = SN
                        st = 'beta = '+ str(beta)
                        plt.figure(22)
                        plt.subplot(2,2,nt)
                        plt.scatter(xx,yy)
                        plt.xlabel(' eigen value eta at:'+ st)
                        plt.ylabel('normalized probability' )
                        plt.hist(x,num_bins, normed= 1.0, facecolor='green',alpha = 0.5)
                        nt = nt+1
           
           
#subplot_plotter(10000,50)   
#subplot_plotter(5000,50)
beta = 0.10
ll = 0
sfailed = 0
ttln = 2000
while ll < ttln :         
     a0,failed = Check_Pendelton_Distribution(beta).pendelton_generation()
     sfailed = sfailed+failed
     print a0
     ll = ll+1
print 'beta = ',beta,',  ','acceptance ratio = ',ttln/float(ttln+sfailed)            