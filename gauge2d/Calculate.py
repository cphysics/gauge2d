
# -*- coding: utf-8 -*-

##########################################################################
## This program calculates the thermal average of Wilson loop "W(K1,K2)" #
## in case of SU(2) gauge theory in 2D with periodic boundary conditions #
## in different representations using "Heat bath algorithm".             #
##                                                                       #
##-----------------------------------------------------------------------#
## Date: 2014-sept-24        Sub: Lattice gauge theory                    #
## By Dibakar Sigdel        Collection: Python-25-06-02                  #
##########################################################################

import numpy as np



class Calculate(object):
            def __init__(self,UU,L,N):
                    self.UU = UU
                    self.L = L
                    self.N = N

            def plqt(self,s,t):
                    D  = [ np.matrix(self.UU[0][s][t]),\
                           np.matrix(self.UU[1][(t+1)%self.L][s]),\
                           np.matrix(self.UU[0][(s+1)%self.L][t]).getH(),\
                           np.matrix(self.UU[1][t][s]).getH()]
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
 

            def wilsonlp(self,K):
                    I = np.matrix(np.identity(self.N))
                    WKK = np.matrix(np.identity( self.N))
                    PD = [[I for k in range(K)]for p in range(4)]
                    DD = [I for k in range(4)]
                    for s  in range(K):
                        PD[0][s] = np.matrix(self.UU[0][0][s])
                    for s in range(K):
                        PD[1][s] = np.matrix(self.UU[1][K][s])
                    for s in range(K):
                        t = K-s-1
                        PD[2][s] = np.matrix(self.UU[0][K][t]).getH()
                    for s in range(K):
                        x = K-s-1
                        PD[3][s] = np.matrix(self.UU[1][0][x]).getH()
                    for r in range(4):
                        for k in range(K):
                            DD[r] = np.dot(DD[r],PD[r][k])
                    WKK = np.dot(DD[0],np.dot(DD[1],np.dot(DD[2],DD[3])))
                    wilp =  (1.0/float(self.N))*np.trace(WKK).real    
                    return wilp    
         
 
            def polyacove(self):
                    Tx = [np.matrix(np.identity(self.N))  for i in range(self.L)]
                    Ty = [np.matrix(np.identity(self.N)) for i in range(self.L)]
       
                    T1 = np.matrix(np.identity(self.N))
                    T2 = np.matrix(np.identity(self.N))
    
                    for i in range(self.L):
                        Tx[i] = np.matrix(self.UU[0][0][i])
                        Ty[i] = np.matrix(self.UU[1][0][i])
           
                    for k in range(self.L):
                        T1 = np.dot(T1,Tx[k])
                        T2 = np.dot(T2,Ty[k])
       
                    chit1 =  (1.0/float(self.N))*np.trace(T1).real
                    chit2 =  (1.0/float(self.N))*np.trace(T2).real
       
                    return chit1,chit2   
    