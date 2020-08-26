# -*- coding: utf-8 -*-

##########################################################################
## This MODULE generates the initial configuration of SU(2) matrix for   #
##  L by L lattice in two dimension.                                     #
##                                                                       #
##-----------------------------------------------------------------------#
## Date: 2014-Sept-12        Sub: Lattice gauge theory                   #
## By Dibakar Sigdel        Collection: Python-014-09-12                 #
##########################################################################






import random as random
import numpy as np
import math as math
import cmath as cmath







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
         


#------------------------------------------------------------------------         
#A = Start(2,5)
#B = A.hot_start()
#print B         
#         
         
         
