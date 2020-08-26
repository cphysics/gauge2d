import numpy as np


def  CI(NB):
            IN = np.matrix([[complex(0,0) for k in range(NB)]for l in range(NB)])
            for k in range(NB):
                for l in range(NB):
                    if k == l:
                        IN[k,l] = complex(1,0)
            return IN  



class selector(object):
    
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

                       
