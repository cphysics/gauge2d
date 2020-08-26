# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 09:12:45 2014

@author: dibakarsigdel
"""


# -*- coding: utf-8 -*-

##########################################################################
## This MODULE generates the initial configuration of SU(2) matrix for   #
##  L by L lattice in two dimension.                                     #
##                                                                       #
##-----------------------------------------------------------------------#
## Date: Mon Sep 29 15:07:50 2014        Sub: Lattice gauge theory       #
## By Dibakar Sigdel        Collection: Python-014-09-12                 #
##########################################################################





import matplotlib.pyplot as plt
import random as random
import numpy as np
import math as math
import cmath as cmath
#import grosswitt




def fun(s):
              if s ==0:
                  fn = 1
              else:
                  fn = 0
              return fn 
      

def  CI(NB):
            IN = np.matrix([[complex(0,0) for k in range(NB)]for l in range(NB)])
            for k in range(NB):
                for l in range(NB):
                    if k == l:
                        IN[k,l] = complex(1,0)
            return IN  




class Start(object):
            def __init__(self,L,N):
                    self.L = L
                    self.N = N

            def cold_start(self): 
                I = np.matrix(np.identity(self.N))
                UU = [[[I for x in range(self.L)]for y in range(self.L)]for r in range(2)]
                return UU 



            def SU2(self):
                   ai = complex(0,1)
                   r = [random.random(),random.random(),\
                        random.random()] 
                   xi = (math.pi*(2*r[0]-1))
                   theta =0.5*(math.acos(2*r[1]-1))
                   phi = (math.pi*(2*r[2]-1))
                   a = [0.0 for l in range(2)]
                   a = [math.cos(theta)*(cmath.exp(ai*phi)),\
                   math.sin(theta)*(cmath.exp(ai*xi))]
                   su2 = []
                   su2 = np.matrix([[a[0],a[1]],\
                       [-a[1].conjugate(),a[0].conjugate()]])
                   return su2

       
            def su2tosun(self,s,t):
        
                    SUN = CI(self.N)
                    SU2 = self.SU2()
                    SUN[s,s] = SU2[0,0]
                    SUN[s,t] = SU2[0,1]
                    SUN[t,s] = SU2[1,0]
                    SUN[t,t] = SU2[1,1]
                    return SUN

                            
            def sun_gnr(self):
                    SUNM = CI(self.N)
                    s = 1
                    while s < self.N:
                        t = s+1
                        while t < self.N+1:
                            SUN = self.su2tosun(s-1,t-1)
                            SUNM = np.dot(SUNM,SUN)
                            t = t+1
                        s = s+1
                        ZSUN = SUNM
                        return ZSUN
                      
                    
            def hot_start(self):
                    I = np.matrix(np.identity(self.N))
                    UU = [[[I for x in range(self.L)]for y in range(self.L)]for z in range(2)]
         
                    for i in range (2):
                         for j in range(self.L):
                             for k in range(self.L):
                                 SUN = self.sun_gnr()     
                                 UU[i][j][k] = SUN
                    return UU
         



class Selector(object):
    
            def __init__(self,N):
                    self.N = N
        
            def count(self):
                    ct =  0
                    for k in range(self.N):
                        ct = ct + k
                    return ct


            def  select(self):
                    ct = 0
                    s = 0
                    P = [0 for k in range(self.count())]
                    Q = [0 for k in range(self.count())]
                    while s < self.N-1:
                        t = s+1
                        while t < self.N:
                            P[ct] = s
                            Q[ct] = t
                            ct = ct+1
                            t = t+1
                        s = s+1
                    return P,Q

  
            def extractw(self,W,ct):
                    WD =  CI(2)
                    P,Q = self.select()
                    s = P[ct]
                    t = Q[ct]
                    WD[0,0] = W[s,s]
                    WD[0,1] = W[s,t]
                    WD[1,0] = W[t,s]
                    WD[1,1] = W[t,t]
                    return WD
       

      
            def expandv(self,VD,ct):
                    V = CI(self.N)
                    P,Q = self.select()
                    s = P[ct]
                    t = Q[ct]
                    V[s,s] = VD[0,0]
                    V[s,t] = VD[0,1]
                    V[t,s] = VD[1,0]
                    V[t,t] = VD[1,1]
                    return V       


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
                     failed = 0  
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
                             return failed , XX
                         else: 
                            count = 0
                            failed = failed+1
        


    
       
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
                        # fm = self.optimize2(a0,2000)
                        #print fm
                 
                        b0 = r[1] 
                        if (b0 < ff):
                            count = count +1 
                            XX = self.su2(a0)
                            return failed, XX 
                        else:
                            count = 0
                            failed =  failed+1
                            #print "failed"
                             
         
                                   
class Update(object):
    
            def __init__(self,U,L,N):
                     self.U = U
                     self.L = L
                     self.N = N
  

                                     
            def staple(self,r,t,s):
                
                    if r ==0:
                        Q = [1,0,1,1,0,1]
                    elif r==1:
                        Q = [0,1,0,0,1,0]
                  
                    #LK = np.matrix(self.UU[r][t][s])
                    D = [ (self.U[Q[0]][(s+1)% self.L][(t-1) + (self.L*fun(t))]).getH(),\
                         (self.U[Q[1]][(t-1) + (self.L*fun(t))][s]).getH(),\
                         (self.U[Q[2]][s][(t-1) + (self.L*fun(t))]),\
                         (self.U[Q[3]][(s+1)%self.L][t]),\
                         (self.U[Q[4]][(t+1)%self.L][s]).getH(),\
                         (self.U[Q[5]][s][t]).getH()]
           
                    W = np.dot(D[0],np.dot(D[1],D[2])) \
                        + np.dot(D[3],np.dot(D[4],D[5]))
                        
                    LK = self.U[r][t][s]
                    WW = np.dot(LK,W)
                    return WW
       

            def  findZk(self,W,ct):
                    Nn = Selector(self.N)
                    WD = Nn.extractw(W,ct)
                    X = WD[0,0] + (WD[1,1]).conjugate()
                    Y = (WD[0,1]).conjugate() - WD[1,0]
                    k = cmath.sqrt(abs(X)**2 + abs(Y)**2).real
                    x = X/k
                    y = Y/k
                    Z = np.matrix([[(x).conjugate(), - (y).conjugate()] ,[y,x]])
                    return k,Z
                    
                     
                                  
            def link(self,r,t,s,alpha,flip):
                    LK =  self.U[r][t][s]
                    W =  self.staple(r,t,s)
                    Nn = Selector(self.N)
                   
                    V = [CI(self.N) for lt in range(Nn.count())]
                    ct = 0
                    while ct < Nn.count():
                        k,Z = self.findZk(W,ct)
                        if alpha > flip :     
                           failed, XX = PendlCrutz(alpha,k).pendlgnr()
                        else:
                           failed, XX = PendlCrutz(alpha,k).crutzgnr()
                        VD = np.dot(XX,Z)
                        V[ct] = Nn.expandv(VD,ct)
                        W = np.dot(V[ct],W)
                        ct = ct+1
          
                    NU = CI(self.N)
                    
                    for q in range(Nn.count()):   
                       NU = np.dot(NU,V[q])
                    NNU = np.dot(NU,LK)
       
                    self.U[r][t][s] = NNU
       
                    return failed, self.U
                    
                    
                    
                    
 
class Calculate(object):
            def __init__(self,U,L,N):
                    self.U = U
                    self.L = L
                    self.N = N

            def plqt(self,s,t):
                    D  =  [(self.U[0][s][t]),\
                           (self.U[1][(t+1)%self.L][s]),\
                           (self.U[0][(s+1)%self.L][t]).getH(),\
                           (self.U[1][t][s]).getH()]
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
                    wtr =  (1.0/float(self.N))*np.trace(UP).real
                    return wtr 
 
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


def thermalization(N,l,titr,alpha,flip):
            ll = 1
            U = Start(l,N).cold_start()
            lx = []
            savp =[]
            while (ll < titr+1):
                    for s in range(l):
                        for t in range(l):
                            for r in range(2):
                                  failed, U = Update(U,l,N).link(r,s,t,alpha,flip)
                    avp = Calculate(U,l,N).avplqt()
                    print ll, avp
                    lx.append(ll)
                    savp.append(avp)
                    plt.figure(100)
                    plt.scatter(0.0,0.0)
                    plt.scatter(0.0,1,0)
                    plt.scatter(ll,avp)
                    plt.show()
                    ll = ll+1
            return  lx,savp




           

                                     
               
                                    
                                                    

###################################################################           
#Declerations------------------------
N = 3
l = 20
titr = 100
sitr = 20
flip = 1.0
alpha = 10.0

#------------------------------------------
#w_checker(N,l,titr,flip)
#erorbar(N,l,titr,flip)
lx,savp = thermalization(N,l,titr,alpha,flip)
ver = [lx,savp]
wnr('su3t.dat',ver).writer()













         
