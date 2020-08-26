# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 18:56:33 2014

@author: dibakarsigdel
"""

import matplotlib.pyplot as plt
import math as math
import numpy as np
import WNR
#import CSV
import random as random
import cmath as cmath
from scipy import linalg 
#from scipy.special import iv



def fun(s):
      if s ==0:
          fn = 1
      else:
          fn = 0
      return fn 



class Start(object):
    def __init__(self,L):
        self.L = L

    def cold_start(self): 
        I = np.matrix(np.identity(2))
        UU = [[[I for x in range(self.L)]for y in range(self.L)]for r in range(2)]
        return UU 


    def hot_start(self):
         ai =  complex(0,1)
         I = np.matrix(np.identity(2))
         U = [[[I for x in range(self.L)]for y in range(self.L)]for z in range(2)]
         
         for i in range (2):
             for j in range(self.L):
                 for k in range(self.L):
                     
                      t = [random.random(),random.random(),\
                           random.random()]                                     
                      
                      xi = (math.pi*(2*t[0]-1))
                      theta =0.5*(math.acos(2*t[1]-1))
                      phi = (math.pi*(2*t[2]-1))
                      
                      a = [0.0 for l in range(2)]
                      a = [math.cos(theta)*(cmath.exp(ai*phi)),\
                           math.sin(theta)*(cmath.exp(ai*xi))]
                           
                      SU2 = []
                      SU2 = np.matrix([[a[0],a[1]],[-a[1].conjugate(),a[0].conjugate()]])
                      
              
                      
                      U[i][j][k] = SU2
                     
         return U
         
         
         
class PendlCrutz(object):
    
            def __init__(self,alpha,k):
                    self.alpha = alpha
                    self.k = k
        

            def su2(self,a0): 
                    aa = math.sqrt(1 - (a0**2))
    
                    t = [random.uniform(0.0,1.0),\
                        random.uniform(0.0,1.0)]
           
                    theta = math.acos(2.0*t[0]-1)
                    xi = 2.0*math.pi*(t[1])
    
                    a =  [ a0,\
                    aa*math.sin(theta)*math.cos(xi),\
                    aa*math.sin(theta)*math.sin(xi),\
                    aa*math.cos(theta)]

                    XX  = np.matrix([[complex(a[0],a[3]),complex(a[2],a[1])] \
                            ,[complex(-a[2],a[1]),complex(a[0],-a[3])]])
           
                    return XX 
         

               
            def pendlgnr(self):
                     count = 0
                     while (count < 1):
                         r = [ random.uniform(0.0,1.0),\
                         random.uniform(0.0,1.0),\
                         random.uniform(0.0,1.0),\
                         random.uniform(0.0,1.0)]
                 
                         x  = [-(math.log(r[0])/(self.k*self.alpha)),\
                               -(math.log(r[1])/(self.k*self.alpha))]
                 
                         C = (math.cos(2.0*math.pi*r[2]))**2
                         A = x[0]*C
                         delta = x[1]+A
           
                         if (r[3]**2) < (1-(0.5*delta)):
                             a0 = (1- delta)
                             count = count+1
                             XX = self.su2(a0)
                             return XX
                         else: 
                            count = 0
       
       
            def funct(self,a0):
                    return (math.sqrt(1.0 - (a0)**2.0))*(math.exp(self.alpha*self.k*a0))
             
             
    
               
            def find_max(self,ptn,fpt,dx):
                    dt = dx/float(ptn)
                    x = [0.0 for k in range(ptn)]
                    fx = [0.0 for k in range(ptn)]
            
                    x = [(fpt + dt*k) for k in range(ptn)]
            
                    for k in range(ptn):
                        fx[k] = [self.funct(x[k])]
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
         
                 
         
            def nfn(self,a0):
                    den = self.optimize(10)
                    ft = self.funct(a0)/den[0]
                    return ft
             
    
            
       
            def  crutzgnr(self):
                   count = 0
                   failed = 0
                   while (count <1):
                       r = [random.uniform(0.0,1.0),\
                           random.uniform(0.0,1.0)]
                      
                       a0 = (2*r[0] - 1)
                       ff = self.nfn(a0)
                       b0 = r[1] 
                       if (b0 < ff):
                         count = count +1 
                         XX = self.su2(a0)
                         return XX 
                       else:
                        count = 0
                        failed =  failed+1
                           
        
                  
class Update(object):
    
            def __init__(self,U,L):
                     self.U = U
                     self.L = L
                    
  
                                                  
            def staple(self,r,t,s):
                    if r ==0:
                        Q = [1,0,1,1,0,1]
                    elif r==1:
                        Q = [0,1,0,0,1,0]
                    #LK = np.matrix(self.UU[r][t][s])
                    D = [ np.matrix(self.U[Q[0]][(s+1)%self.L][(t-1) + (self.L*fun(t))]).getH(),\
                          np.matrix(self.U[Q[1]][(t-1) + (self.L*fun(t))][s]).getH(),\
                          np.matrix(self.U[Q[2]][s][(t-1) + (self.L*fun(t))]),\
                          np.matrix(self.U[Q[3]][(s+1)%self.L][t]),\
                          np.matrix(self.U[Q[4]][(t+1)%self.L][s]).getH(),\
                          np.matrix(self.U[Q[5]][s][t]).getH()]
           
       
                    W = np.dot(D[0],np.dot(D[1],D[2])) \
                     + np.dot(D[3],np.dot(D[4],D[5]))
                    k = math.sqrt(linalg.det(W).real)
                    A = W/(k)
                    return k,A 
    
                                                       
            def link(self,r,t,s,alpha,flip):
                    k,A = self.staple(r,t,s)
                    AD = A.getH()
                    if alpha > flip:
                        XX = PendlCrutz(alpha,k).pendlgnr()
                    else:
                        XX = PendlCrutz(alpha,k).crutzgnr()
                    NU = np.dot(XX,AD)
                    self.U[r][t][s] = NU
                    return self.U 
 





class Calculate(object):
            def __init__(self,U,L):
                    self.U = U
                    self.L = L
                   

            def plqt(self,s,t):
                    D  = [ np.matrix(self.U[0][s][t]),\
                           np.matrix(self.U[1][(t+1)%self.L][s]),\
                           np.matrix(self.U[0][(s+1)%self.L][t]).getH(),\
                           np.matrix(self.U[1][t][s]).getH()]
                    return D         


            def  avplqt(self):
                    sum_trp = 0.0  
                    for s  in range(self.L):
                        for t in range(self.L):
                            D = self.plqt(s,t)
                            UP = np.dot(D[0],np.dot(D[1],np.dot(D[2],D[3])))
                            trup = (1.0 - ((1.0/float(self.N))*np.trace(UP).real))
                            sum_trp = sum_trp + (trup/float(self.L*self.L))
                    return sum_trp  
 

            def wloop11(self,s,t):
                    D = self.plqt(s,t)       
                    UP = np.dot(D[0],np.dot(D[1],np.dot(D[2],D[3])))
                    wtr =  (1.0/float(2))*np.trace(UP).real
                    return wtr 





def Mean_Error(stor_w11):
        nt = len(stor_w11)
        ver = [0.0 for k in range(nt)] 
        mw11 = 0.0
         
        for k in range (nt):
            mw11 = mw11+stor_w11[k]/float(nt)
        for l in range (nt):
            ver[l] = (stor_w11[l]-mw11)**2
        #print ver
        s_error = math.sqrt(sum(ver)/float(nt*nt))
        return  mw11, s_error
        
        



def exportstorw(l,titr,alpha,flip):
            sitr=50
            storw = [0.0 for k in range(titr - sitr)]
            ll = 1
            U = Start(l).cold_start()
            while (ll < titr+1): 
                    
                    for s in range(l):
                        for t in range(l):
                            for r in range(2):
                                U = Update(U,l).link(r,s,t,alpha,flip)
                    w11 = Calculate(U,l).wloop11(5,5)
                    print alpha,w11
                    if ll > sitr:
                        storw[ll-sitr-1] = w11
                    ll = ll+1
             
            return  storw


def errorbar(N,l,titr,flip):
            plot_dir = "/Users/dibakarsigdel/Dropbox/Plots/"  
            data_dir = "/Users/dibakarsigdel/Dropbox/Data/" 
            Nmax = 28
            Nvalue = 1
            dlamda = 0.25
            x = [0.0 for k in range(Nmax)]
            y = [0.0 for k in range(Nmax)]
            y_error = [0.0 for k in range(Nmax)]
            while Nvalue < Nmax+1:
                    lamda =  dlamda*Nvalue
                    alpha = (2.0*N)/lamda
                    x[Nvalue-1] =  (2.0*N)/alpha
                    storw = exportstorw(l,titr,alpha,flip) 
                    y[Nvalue-1],y_error[Nvalue-1] = Mean_Error(storw) 
                    st = str(N)
                    plt.figure(104)
                    plt.xlabel('lambda')
                    plt.ylabel('W11')
                    plt.grid()
                    plt.errorbar(x,y, yerr = y_error, fmt='8')
                    plt.savefig(plot_dir + 'plotsu'+st+'.png')
                    plt.show()
                    WNR.wnr('su'+st+'.dat',[x,y,y_error]).writer()
                    WNR.wnr(data_dir +'su'+st+'.dat',[x,y,y_error]).writer()
                    Nvalue = Nvalue+1
            
            return  
                                                    
                    
############################################################################           
#Declerations------------------------
N = 2                  
l = 10
titr = 1000
flip = 0.5

#------------------------------------------
errorbar(N,l,titr,flip)