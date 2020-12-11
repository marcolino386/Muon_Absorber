#Program to create histogram to all data and also a histogram with the sum of all angles for each energy.

import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings("ignore")


def make_hist(A,wdt):
    Emin,Emax = A.min(),A.max()
    bs = np.arange((Emin-wdt),(Emax+wdt),wdt)
    hist = plt.hist(A,bins=list(bs))
    freq = hist[0]
    dE = hist[1]
    #plt.show()
    return freq,dE

#angles = [1.432,1.909, 2.291, 2.726,3.013,3.576, 4.399, 5.194, 5.711, 6.34, 7.125, 8.13, 8.746, 9.019, 9.462, 9.782, 9.951, 10.3, 10.89, 11.31, 11.77, 12.53, 15.12]
angles =  [2.291, 2.862, 4.764, 6.34, 7.125, 8.13, 9.462]
Es =  [1.5, 1.7, 1.9, 2.1, 2.3, 2.4, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9]
width = 0.003 # Interval width in GeV

drt = os.getcwd()   # Get current working directory
#charge = ['plus', 'minus']
charge = ['plus']


for j in charge:
   ang_index=0
   for Ei in Es:
      #all_E_L = [] # List for histogram of all energies (independent of angle)
      # Loop over angles:
      #fn = '{4}/{2}/{2}{3}/Energy{0}_Angle_{1}.dat'.format(Ei,t,pastas[t], charge[i],drt) # Data filename
      fn = '{2}/data_mu_plus2/Energy{0}_Angle_{1}.dat'.format(Ei,angles[ang_index],drt)
      ang_index += 1  
      data = np.loadtxt(fn,dtype=float)
      if np.size(data) > 7: # pois se for = ou menor a 7 quer dizer que não ha pontos de muon no detetor
         E = data[:,3]
         #all_E_L.append(list(E))

         freq,dE = make_hist(E,width)
         arq = open('{2}/histogram_Energy{0}_Angle_{1}_{3}.dat'.format(Ei,angles[ang_index],drt,j),'w')
         for i,f in enumerate(freq):
            arq.write('{0} {1} \n'.format(dE[i],int(f)))
         arq.close()
      else:
         print('Sem valores para E = {0} Gev e theta = {1}'.format(Ei,angles[ang_index])) 
       
  
            
