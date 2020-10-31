# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 10:13:23 2014

@author: dibakarsigdel
"""
import CSV

def betavalue(N):
        lamda = [0.1*k for k in range(42)]
        beta = [float(N)/(lamda[k]) for k in range(1,41)]
        #print beta
        return beta
        
lamda = [0.1*k for k in range(1,41)]        
U = [lamda, betavalue(2),betavalue(3),betavalue(4),betavalue(5),betavalue(6)]

CSV.CSV('beta.csv').csv_writer(U)








