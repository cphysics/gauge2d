
# -*- coding: utf-8 -*-

##########################################################################
## This program calculates the thermal average of Wilson loop "W(K1,K2)" #
## in case of SU(2) gauge theory in 2D with periodic boundary conditions #
## in different representations using "Heat bath algorithm".             #
##                                                                       #
##-----------------------------------------------------------------------#
## Date: 2014-Jun-25        Sub: Lattice gauge theory                    #
## By Dibakar Sigdel        Collection: Python-25-06-02                  #
##########################################################################


import numpy as np
import Pendelton
import cmath as cmath
import Selector


def fun(s):
      if s ==0:
          fn = 1
      else:
          fn = 0
      return fn 
      
      
      
                                    
class update(object):
    
            def __init__(self,UU,L,N):
                     self.UU = UU
                     self.L = L
                     self.N = N
  
                                                  
            def staple(self,r,t,s):
                    if r ==0:
                        Q = [1,0,1,1,0,1]
                    elif r==1:
                        Q = [0,1,0,0,1,0]
                    LK = np.matrix(self.UU[r][t][s])
                    D = [ np.matrix(self.UU[Q[0]][(s+1)%self.L][(t-1) + (self.L*fun(t))]).getH(),\
                          np.matrix(self.UU[Q[1]][(t-1) + (self.L*fun(t))][s]).getH(),\
                          np.matrix(self.UU[Q[2]][s][(t-1) + (self.L*fun(t))]),\
                          np.matrix(self.UU[Q[3]][(s+1)%self.L][t]),\
                          np.matrix(self.UU[Q[4]][(t+1)%self.L][s]).getH(),\
                          np.matrix(self.UU[Q[5]][s][t]).getH()]
           
       
                    W = np.dot(D[0],np.dot(D[1],D[2])) \
                     + np.dot(D[3],np.dot(D[4],D[5]))
       
                    WW = np.dot(LK,W)
                    return WW
       

            def  findZk(self,W,ct):
                    Nn = Selector.selector(self.N)
                    WD = Nn.extractw(W,ct)
                    X = WD[0,0] + (WD[1,1]).conjugate()
                    Y = (WD[0,1]).conjugate() - WD[1,0]
                    k = cmath.sqrt(abs(X)**2 + abs(Y)**2).real
                    x = X/k
                    y = Y/k
                    Z = np.matrix([[(x).conjugate(), - (y).conjugate()] ,[y,x]])
                    return k,Z
      
                                                       
            def link(self,r,t,s,alpha):
                    LK =  np.matrix(self.UU[r][t][s])
                    W =  self.staple(r,t,s)
                    Nn = Selector.selector(self.N)
                    Cn = Selector
                    V = [Cn.CI(self.N) for lt in range(Nn.count())]
                    ct = 0
                    while ct < Nn.count():
                        k,Z = self.findZk(W,ct)
                        XX = Pendelton.Pendelton(alpha,k).pendlgnr()
                        VD = np.dot(XX,Z)
                        V[ct] = Nn.expandv(VD,ct)
                        W = np.dot(V[ct],W)
                        ct = ct+1
          
                    NU = Cn.CI(self.N)
                    for q in range(Nn.count()):   
                       NU = np.dot(NU,V[q])
                    NNU = np.dot(NU,LK)
       
                    self.UU[r][t][s] = NNU
       
                    return self.UU 
 


      

   
