# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 10:31:29 2014

@author: dibakarsigdel
"""

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
         
