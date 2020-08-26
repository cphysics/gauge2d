


import matplotlib.pyplot as plt
import math as math
import numpy as np
import random as random
import cmath as cmath
from scipy import linalg 
#from scipy.special import iv





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
         
         
         
class Pendl_Crutz(object):
    
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
                       # fm = self.optimize2(a0,2000)
                       #print fm
                 
                       b0 = r[1] 
                       if (b0 < ff):
                         count = count +1 
                         XX = self.su2(a0)
                         return XX 
                       else:
                        count = 0
                        failed =  failed+1
                        #print "failed = ",failed         
def fun(s):
      if s ==0:
          fn = 1
      else:
          fn = 0
      return fn 

      
                  
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
    
                                                       
            def link(self,r,t,s,beta, flip):
                    k,A = self.staple(r,t,s)
                    AD = A.getH()
                    if beta > flip:
                        XX = Pendl_Crutz(beta,k).pendlgnr()
                    else:
                        XX = Pendl_Crutz(beta,k).crutzgnr()
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
                            trup = (1.0 - ((1.0/float(2))*np.trace(UP).real))
                            sum_trp = sum_trp + (trup/float(self.L*self.L))
                    return sum_trp  
 

            def wloop11(self,s,t):
                    D = self.plqt(s,t)       
                    UP = np.dot(D[0],np.dot(D[1],np.dot(D[2],D[3])))
                    wtr =  (1.0/float(2))*np.trace(UP).real
                    return wtr 
                    
                    
                    
            def wilsonlp(self,K1,K2):
                    I = np.matrix(np.identity(2))
                    WKK = np.matrix(np.identity(2))
                    nn = [K1,K2,K1,K2]
                    PD = [[I for k in range(nn[r])]for r in range(4)]
                    DD = [I for k in range(4)]
                    
                    for s  in range(K1):
                        PD[0][s] = np.matrix(self.U[0][0][s])
                    for s in range(K2):
                        PD[1][s] = np.matrix(self.U[1][K2][s])
                    for s in range(K1):
                        t = K1-s-1
                        PD[2][s] = np.matrix(self.U[0][K1][t]).getH()
                    for s in range(K2):
                        x = K2-s-1
                        PD[3][s] = np.matrix(self.U[1][0][x]).getH()
                        
                    for r in range(4):
                        ct = len(PD[r])
                        for k in range(ct):
                            DD[r] = np.dot(DD[r],PD[r][k])
                    WKK = np.dot(DD[0],np.dot(DD[1],np.dot(DD[2],DD[3])))
                    wilp =  (1.0/float(2))*np.trace(WKK).real    
                    return wilp    
         

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
        
        

#############################################################################

def thermalization(l,titr,beta,flip):
    
        ll = 0
        U = Start(l).cold_start()
        while (ll < titr+1):
                for s in range(l):
                    for t in range(l):
                        for r in range(2):
                           U = Update(U,l).link(r,s,t,beta,flip)
                avp = Calculate(U,l).avplqt()
                print ll, avp
                plt.figure(100)
                #plt.scatter(0,0.0)
                #plt.scatter(0,1.0)
                plt.grid()
                plt.xlabel('iterations')
                plt.ylabel('average plaquette')
                
                plt.text(101,0.175,'beta = 8.0')
                #plt.text(101,0.235,'beta = 6.0')
                #plt.text(101,0.34,'beta = 4.0')
                #plt.text(101,0.57, 'beta = 2.0')
                plt.scatter(0.0,1.0)
                plt.scatter(ll,avp)
                plt.show()
                ll = ll+1 
       


def Mean_Error(stor_w11):
        nt = len(stor_w11)
        ver = [0.0 for k in range(nt)] 
        mw11 = 0.0
         
        for k in range (nt):
            mw11 = mw11+stor_w11[k]/float(nt)
        for l in range (nt):
            ver[l] = (stor_w11[l]-mw11)**2
        s_error = math.sqrt(sum(ver)/nt**2)
        return  mw11, s_error
        
        
def export_wstor(l,ct,titr,beta,flip,K1,K2):
        sitr =10
        # titr = total no of iterations(100)
        # sitr = start counting iterations for average value(20)
        citr = titr-sitr
        wstor = [0.0 for k in range (citr)]
        ll = 0
        U = Start(l).cold_start()
        while (ll < titr+1):
                
                for s in range(l):
                    for t in range(l):
                        for r in range(2):
                           U = Update(U,l).link(r,s,t,beta,flip)
              
                if ll > sitr:
                        wk1k2 = Calculate(U,l).wilsonlp(K1,K2)
                        wk1_1k2 = Calculate(U,l).wilsonlp(K1-1,K2)
                        wk1k2_1= Calculate(U,l).wilsonlp(K1,K2-1)
                        wk1_1k2_1 = Calculate(U,l).wilsonlp(K1-1,K2-1)
                        wtr  = -math.log(abs((wk1k2*wk1_1k2_1)/(wk1k2_1*wk1_1k2)))
                        print ct,ll, wtr
                        rt = ll- sitr-1
                        wstor[rt] = wtr
                ll = ll+1 
        return wstor
                


        
def errorbar(l,titr,tdt,flip,K1,K2):
    # tdt = total no of data points(value of alpha)  taken for plot
    plot_dir = "/Users/dibakarsigdel/Dropbox/Plots/"  
    data_dir = "/Users/dibakarsigdel/Dropbox/Data/" 
    ct =1
    ax = [0.0 for k in range (tdt)]
    ay = [0.0 for k in range (tdt)]
    y_error = [0.0 for k in range (tdt)]
    st = str(2) 
    while (ct < tdt+1):
        beta = ct*0.1
        ax[ct-1] = beta
        
        wstor = export_wstor(l,ct,titr,beta,flip,K1,K2)
        ay[ct-1],y_error[ct-1] =  Mean_Error(wstor)
        plt.figure(102)
        plt.xlabel('beta')
        plt.ylabel('crutz_ratio')
        plt.grid()
        plt.errorbar(ax,ay, yerr = y_error, fmt='8')
        plt.savefig(plot_dir + 'plotsu'+st+'.png')
        plt.show()
        wnr('su'+st+'.dat',[ax,ay,y_error]).writer()
        wnr(data_dir +'su'+st+'.dat',[ax,ay,y_error]).writer()
        ct = ct+1
    return ax,ay,y_error
                    
############################################################################           
l = 30
tdt = 10
titr = 40
flip = 0.5
#------------------------------------------
errorbar(l,titr,tdt,flip,4,4)
#thermalization(l,titr,8.0,flip)
















