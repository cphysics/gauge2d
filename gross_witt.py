
import matplotlib.pyplot as plt
import csv






class CSV(object):
    
        def __init__(self,data_file):
                    self.data_file = data_file

        def csv_reader(self,nr,nc,nnr,nnc):
    
                with open( self.data_file,'rb') as csvfile:
                         reader = csv.reader(csvfile)
    
    
                         Bigdata = [[0 for l in range(nr)]for k in range (nc)]
     
                         i = 0
                         for row in reader:
                             for col in range(nc):
                                 Bigdata[col][i] = row[col]
                             i = i+1
      
                         Puredata = [[0 for l in range(nnr)]for k in range (nnc)] 
     
                         for r in range (nnr):
                             for c in range(nnc):
                                     Puredata[c][r] = Bigdata[c][r+1]
                         return Puredata
         


        def csv_writer(self,ver_to_write):
          
                          xn = len(ver_to_write)
                          yn = len(ver_to_write[0])
                          Inverse = [[0 for r in range(xn)]for c in range(yn)]
         
                          for r in range(yn):
                             for c in range(xn):
                                  Inverse[r][c] = ver_to_write[c][r] 
                            
                          with open(self.data_file, 'wb') as csvfile:
                             writer = csv.writer(csvfile)
                             writer.writerows(Inverse)
                          return
         


         


################## grosswitt plot ###################################



class Grosswittpr(object):
        # N = 16 points
        #dx = 0.25 ----> x - interval<----
    
        def __init__(self,dx,N):
            self.dx = dx
            self.N  = N
    
        
        def fx(self,x):
                if x > 2 or x == 2:
                    fnx = 1/x
                elif x < 2:
                    fnx = 1-x/4.0
                return fnx
    
      
        def plot_rec(self):
            
                x = [0.0 for k in range(self.N)]
                y = [0.0 for k in range(self.N)]

                for k in range (self.N):
                    xt = self.dx*k
                    yt = self.fx(xt)
                    x[k] = xt
                    y[k] = yt
    
    
                plt.figure(502)
                plt.plot(x,y,'-')
                plt.show()
                U = [x,y]
                CSV('wlamda.csv').csv_writer(U)
                return 

                
                                                
                
############ tabular beta ###########################################


class Beta_Tabular(object):
    
        def __init__(self,dx,dtp):
               self.dx = dx
               self.dtp = dtp

        def bet(self,N):
            
                #dtp = data points = 16
                #dx = 0.25
                
                lamda = [self.dx*k for k in range(1,self.dtp+1)]
                beta = [float(N*N)/(lamda[k]) for k in range(self.dtp)]
                return beta
                
        def betatable(self):
                lamda = [self.dx*k for k in range(1,self.dtp+1)]
                        
                U = [lamda, self.bet(2),self.bet(3),self.bet(4),self.bet(5),self.bet(6)]
                CSV('beta.csv').csv_writer(U)
                return 




###################### take action ######################################


Grosswittpr(0.25,16).plot_rec()
#-------------------------------------------
Beta_Tabular(0.25,16).betatable()

