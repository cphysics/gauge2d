
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





import matplotlib.pyplot as plt
import Start
import Update
import Calculate

#Declerations------------------------
N = 3
l = 30
itr = 10
alpha= 10.0
#------------------------------------------

U = Start.Start(l,N).cold_start() 

ll = 1

while (ll < itr+1): 
     for s in range(l):
           for t in range(l):
               for r in range(2):
                   
                      U = Update.update(U,l,N).link(r,s,t,alpha)
              
              
                                  
     avp = Calculate.Calculate(U,l,N).avplqt()  
     #avp = Wilson11(U0)            
     print avp
     plt.figure(1)
     plt.scatter(ll,avp)
     plt.show()
        
     ll = ll+1 
   
   
   

 
