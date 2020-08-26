




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
                                
def plot_me(pn,xx,yy):
        st = 'SU('+str(pn+2)+')'
        xt,yt = largeN_limit(100,0.04)                              
        plt.subplot(3,2,pn)
        plt.text(2.5, 0.9, st)
        plt.text(2.5, 0.8,'L=30')
        plt.xlabel('lambda = g^2*N')
        plt.ylabel('wilson_loop')
        plt.plot(xt,yt,'-')    
        plt.scatter(xx[pn],yy[pn])
        plt.show()    
        return
        
        
        
##########################################################################

x = [[0.0 for k in range(16)]for l in range(6)] 
y = [[0.0 for k in range(16)] for l in range(6)] 
xx = [[0.0 for k in range(16)]for l in range(6)] 
yy = [[0.0 for k in range(16)] for l in range(6)] 

data_dir = "/Users/dibakarsigdel/Documents/Python/Lab/gross_witten/"

for k in range(6):    
   [xx[k] ,yy[k] ] = WNR.wnr( data_dir + 'su'+ str(k+2)+'.dat',[x[0],y[0]]).reader()   
      
         
for k in range(6):
       plot_me(k,xx,yy)               
       
    
