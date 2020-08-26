import random as random
import numpy as np
import math as math





class Pendelton(object):
    
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
        


