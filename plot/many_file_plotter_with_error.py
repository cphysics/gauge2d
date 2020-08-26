# -*- coding: utf-8 -*-
'''
############################################################################
## This program PLOTES many set of [x,y] data( ex:set 6) in Subplot manner.#
## Data should be available as "su" + n + ".dat" format for variable n     #
##                                                                         #
##-------------------------------------------------------------------------#
## Date: 2014-OCT-1        Sub: Lattice gauge theory                       #
## By Dibakar Sigdel        Collection: Python-1-10-02                     #
############################################################################
'''





import WNR
import matplotlib.pyplot as plt


            
def fx(x): 
      if x > 2 or x == 2:
            fnx = 1/x
      elif x < 2:
            fnx = 1-x/4.0
      return fnx
      

def largeN_limit(N,dlambda):
    
        x = [0.0 for k in range(N)]
        y = [0.0 for k in range(N)]

        for k in range (N):
                xt = dlambda*k
                yt = fx(xt)
                x[k] = xt
                y[k] = yt
        return  x,y
                                
def plot_me(pn,xx,yy,ye):
        st = 'SU('+str(pn+3)+')'
        xt,yt = largeN_limit(300,0.025)                              
        plt.subplot(3,2,pn)
        plt.text(2.5, 0.9, st)
        plt.text(2.5, 0.8,'L=30')
        plt.xlabel('lambda = g^2*N')
        plt.ylabel('wilson_loop')
        plt.plot(xt,yt,'-')    
        plt.errorbar(xx[pn],yy[pn],yerr = ye[pn],fmt = '8')
        plt.show()    
        return
        
        
   
##########################################################################
ND = 28
NP = 5
 
x = [0.0 for k in range(ND)]
y = [0.0 for k in range(ND)] 
er =  [0.0 for k in range(ND)] 
xx = [[0.0 for k in range(ND)]for l in range(NP)] 
yy = [[0.0 for k in range(ND)] for l in range(NP)] 
ye = [[0.0 for k in range(ND)] for l in range(NP)] 

data_dir = "plot-lab/"

for k in range(NP):    
   [xx[k],yy[k],ye[k]] = WNR.wnr( data_dir + 'su'+ str(k+3)+'.dat',[x,y,er]).reader()   
      
         
for k in range(NP):
       plot_me(k,xx,yy,ye)               
