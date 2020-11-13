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
Es = [1,2,3,4,5,6,7,8,9,10]
pastas = [1.4,1.9,2.2,2.7,3.0,3.5,4.3,5.1,5.7,6.3,7.1,8.1,8.7,9.0,9.4,9.7,9.9,10.3,10.8,11.3,11.7,12.5,15.1]
angles = [1.432,1.909, 2.291, 2.726,3.013,3.576, 4.399, 5.194, 5.711, 6.34, 7.125]
width = 0.003 # Interval width in GeV

drt = os.getcwd()   # Get current working directory
charge = ['plus', 'minus']

for j in charge:
    for Ei in Es:
        all_E_L = [] # List for histogram of all energies (independent of angle)
        # Loop over angles:
        for t in angles:
            #fn = '{4}/{2}/{2}{3}/Energy{0}_Angle_{1}.dat'.format(Ei,t,pastas[t], charge[i],drt) # Data filename
            fn = '{2}/data_mu_plus1/Energy{0}_Angle_{1}.dat'.format(Ei,t,drt)
            #fn = '/home/lgp/absorber/data/10k/10k05T/data_mu_{2}1/Energy{0}_Angle_{1}.dat'.format(Ei,t,j)
            data = np.loadtxt(fn,dtype=float)
            if np.size(data) > 7: # pois se for = ou menor a 7 quer dizer que não ha pontos de muon no detetor
                E = data[:,3]
                all_E_L.append(list(E))

                freq,dE = make_hist(E,width)
                arq = open('{2}/histogram_Energy{0}_Angle_{1}_{3}.dat'.format(Ei,t,drt,j),'w')
                for i,f in enumerate(freq):
                    arq.write('{0} {1} \n'.format(dE[i],int(f)))
                arq.close()
            else:
                print('Sem valores para E = {0} Gev e theta = {1}'.format(Ei,t)) 
        all_E = np.array(all_E_L)
        
        if all_E_L != []:
            freq,dE = make_hist(np.hstack(all_E),width)
            arq = open('{1}/histogram_Energy{0}_Angle_all_{2}.dat'.format(Ei,drt, j),'w')
            for i,f in enumerate(freq):
                arq.write('{0} {1} \n'.format(dE[i],int(f)))
            arq.close()
  
            
