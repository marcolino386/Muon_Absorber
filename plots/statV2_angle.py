#Calculate mean, variation of all the data for every angle.

import numpy as np
import os
import glob
import matplotlib.pyplot as plt


def end_calculus(X0,Y0,E0,px0,py0,pz0):
    print("cabo")

def pt_array(px,py,pz):
    A = []
    for i in range(np.size(px)):
       pti = np.sqrt(px[i]*px[i] + py[i]*py[i])
       A.append(pti)
    return A

def theta_array(px,py,pz):
    A = []
    for i in range(np.size(px)):
       thetai = np.arccos((pz[i])/np.sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]))
       A.append(thetai)
    return A
    
def modify_array():
    X0 = X[0]; Y0 = Y[0]; Z0 = Z[0];  E0 = E[0]; px0 = px[0]; py0 = py[0];  pz0 = pz[0]
    X.remove(0); Y.remove(0); Z0.remove(0); E.remove(0); px.remove(0); py.remove(0); pz.remove(0)

def make_calculus(ang, energy, charge,path,index,data_momentum,data_position,field):
    X0,Y0,E0,px0,py0,pz0 = 0.0, 0.0, 0.0 ,0.0,0.0,0.0
    #deleta a / do argumento path
    data_name = path.replace("/", "")
    #abre dois arquivos para momentum e posição
    #entra em todos os arquivos .dat da pasta da variavel path

    #filename = open('/home/lgp/absorber/data/10k/10k{5}T/data_mu_{2}1/Energy{3}_Angle_{4}.dat'.format(path,ang,charge,energy,angles[index],field),"r")
     
    filename1 = open("{0}/data_mu_plus1/Energy{3}_Angle_{4}.dat".format(path,ang,charge,energy,angles[index]),"r")
    filename = open("{0}/data_mu_plus2/Energy{3}_Angle_{4}.dat".format(path,ang,charge,energy,angles[index]),"r")
    #transforma o arquivo em um array
    A = np.loadtxt(filename, dtype = float)
    B = np.loadtxt(filename1, dtype = float)
    
    Pxb,Pyb = B[:,4],B[:,5]
    Pxb0 = Pxb[0]; Pyb0 = Pyb[0]
    #se o numero de elementos for 7 significa que nenhum muon passou, apenas há o muon do evento sem absorber
    if np.size(A) <= 7:
       #seleciona a variavel de sua respectiva coluna da matriz(array)
       X,Y,E,px,py,pz = 0,0,0,0,0,0
       slopeX ,slopeY = 0, 0
       theta0 = 0
       if np.size(A) == 7:
          X,Y,E,px,py,pz = A[1],A[2],A[3],A[4],A[5], A[6]
          theta0 = np.arccos((pz)/np.sqrt(px*px + py*py + pz*pz))*(180.00/3.14159265)
          slopeX = px/pz
          slopeY=py/pz
          DeltaPx = px - Pxb0
          DeltaPy = py - Pyb0
                

       #Printa

       print("----------------------------------- " +str(energy)+ "GEV--------------------------------------------------------------------------------")
       print("E0, X0, Y0, Xmed, Ymed, ENmed ,dx, dy, dEn, desX, desY,desen, final0")
       print(str(E) + " " + str(X) + " " + str(Y) + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " "+ "0" +"\n")
       print("######################################")
       print("en[0], px[0], py[0],theta0, pt0 , pxmed, pymed, ptmed , thetamed , dpx, dpy, dpt, dtheta, despx, despy,slopeX0,slopeY0,slopeXmed,slopeYmed,dslopeX,dslopeY, destheta, despt, DeltaPx, DeltaPy ,final0")
       print(str(E) + " " + str(px) + " " + str(py) +  " " + str(theta0) + " " + str(np.sqrt(px*px + py*py + pz*pz)) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(slopeX) + " " + str(slopeY) + " " + str(0) + " " + str(0) + " " +  str(0) + " " + str(0) + " " +  str(0) + " " + str(0) + " " + str(DeltaPx) + " " + str(DeltaPy)+ " " + str(0) +  "\n")

       data_position.write(str(E) + " " + str(X) + " " + str(Y) + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + "\n")


       data_momentum.write(str(E) + " " + str(px) + " " + str(py) +  " " + str(theta0) + " " + str(np.sqrt(px*px + py*py + pz*pz)) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(slopeX) + " " + str(slopeY) + " " + str(0) + " " + str(0) + " " +  str(0) + " " + str(0) + " " +  str(0) + " " + str(0) + " " + str(DeltaPx) + " " + str(DeltaPy)+ " " + str(0) +  "\n")
           
    else:
       #seleciona a variação de sua respectiva coluna da matriz(array)
       X,Y,E,px,py,pz = A[:,1],A[:,2],A[:,3],A[:,4],A[:,5],A[:,6]

       #modify_array()
       X0 = X[0]; Y0 = Y[0];   E0 = E[0]; px0 = px[0]; py0 = py[0];  pz0 = pz[0]
       X = np.delete(X,0); Y = np.delete(Y,0); E = np.delete(E,0); px = np.delete(px,0); py = np.delete(py,0); pz = np.delete(pz,0);
       
  
       
        
           
       #cria as arrays de pt e theta 
       pt = pt_array(px,py,pz)
       theta = theta_array(px,py,pz)

       #calcula media e desvio padrão utilizando a bilbioteca
       Xmed = np.mean(X); Ymed = np.mean(Y); Emed = np.mean(E); pxmed = np.mean(px); pymed = np.mean(py); pzmed = np.mean(pz); thetamed = np.mean(theta)*(180.00/3.14159265);ptmed = np.mean(pt)
       desX = np.std(X); desY = np.std(Y); desE = np.std(E); despx = np.std(px); despy = np.std(py); despz = np.std(pz); despt = np.std(pt); destheta = np.std(theta)

       #calculando a variavel slope
       slopeXmed = pxmed/pzmed
       slopeYmed = pymed/pzmed
       slopeX0  = px0/pz0
       slopeY0 = py0/pz0
       dslopeX = slopeXmed - slopeX0
       dslopeY = slopeYmed - slopeY0
 
       #desvio
       dx = Xmed - X0; dy = Ymed - Y0; dE = Emed - E0; dpx = pxmed - px0; dpy = pymed - py0;  dslopeX = slopeXmed - slopeX0; dslopeY = slopeYmed - slopeY0
  
       

       #modulo do momentum e dpt
       pt0 = np.sqrt(px0*px0 + py0*py0)
       dpt = ptmed - pt0

       #angulo do momentum 
       theta0 = np.arccos((pz0)/np.sqrt(px0*px0 + py0*py0 + pz0*pz0))*(180.00/3.14159265)
       dtheta = thetamed - theta0
         
       DeltaPx = pxmed - Pxb0
       DeltaPy = pymed - Pyb0   
     
    
       #Printa dados

       # E0, X0, Y0, Xmed, Ymed, ENmed ,dx, dy, dEn, desX, des , desx, desy,desen, final0
       final = np.size(E)
       print("----------------------------------- " +str(E0)+ "GEV--------------------------------------------------------------------------------")
       print("E0, X0, Y0, Xmed, Ymed, ENmed ,dx, dy, dEn, desX, desY,desen, final0, theta0")
       print(str(E0) + " " + str(X0) + " " + str(Y0) + " " + str(Xmed) + " " + str(Ymed)  + " " + str(Emed) + " " + str(dx) + " " + str(dy) + " " + str(dE) + " " + str(desX) + " " + str(desY) + " " + str(desE) + " " + str(final) + " " + str(theta0) + "\n") 
       print("######################################")
       #en[0], px[0], py[0],theta0, pt0 , pxmed, pymed, ptmed , thetamed , dpx, dpy, dpt, dtheta, despx, despy,slopeX0,slopeY0,slopeXmed,slopeYmed,dslopeX,dslopeY, final0
       print("en[0], px[0], py[0],theta0, pt0 , pxmed, pymed, ptmed , thetamed , dpx, dpy, dpt, dtheta, despx, despy,slopeX0,slopeY0,slopeXmed,slopeYmed,dslopeX,dslopeY, final0")
       print(str(E0) + " " + str(px0) + " " + str(py0) +  " " + str(theta0) + " " + str(pt0) + " " + str(pxmed) + " " + str(pymed) + " " + str(ptmed) + " " + str(thetamed) + " " + str(dpx) + " " + str(dpy) + " " + str(dpt) + " " + str(dtheta) + " " + str(despx) + " " + str(despy) + " " + str(slopeX0) + " " + str(slopeY0) + " " + str(slopeXmed) + " " + str(slopeYmed) + " " + str(dslopeX) + " " + str(dslopeY) + " " + str(destheta) + " " + str(despt) + " " + str(DeltaPx) + " " + str(DeltaPy)+ " " + str(final) +  "\n") 
           
          


       #escreve nos arquivos de saída

       # "E0, X0, Y0, Xmed, Ymed, ENmed ,dx, dy, dEn, desX, desY,desen, final0"

       data_position.write(str(E0) + " " + str(X0) + " " + str(Y0) + " " + str(Xmed) + " " + str(Ymed) + " " + str(Emed) + " " + str(dx) + " " + str(dy) + " " + str(dE) + " " + str(desX) + " " + str(desY) + " " + str(desE) + " " + str(final) + " " + str(theta0) + "\n") 

       #"en[0], px[0], py[0],theta0, pt0 , pxmed, pymed, ptmed , thetamed , dpx, dpy, dpt, dtheta, despx, despy,slopeX0,slopeY0,slopeXmed,slopeYmed,dslopeX,dslopeY, destheta, despt, DeltaPx, DeltaPy ,final0"

       data_momentum.write(str(E0) + " " + str(px0) + " " + str(py0) +  " " + str(theta0) + " " + str(pt0) + " " + str(pxmed) + " " + str(pymed) + " " + str(ptmed) + " " + str(thetamed) + " " + str(dpx) + " " + str(dpy) + " " + str(dpt) + " " + str(dtheta) + " " + str(despx) + " " + str(despy) + " " + str(slopeX0) + " " + str(slopeY0) + " " + str(slopeXmed) + " " + str(slopeYmed)+ " " + str(dslopeX) + " " + str(dslopeY) + " " + str(destheta) + " " + str(despt) + " " + str(DeltaPx) + " " + str(DeltaPy)+ " " + str(final) +  "\n") 


    print("DONE")   
          
#charge = ['plus', 'minus']      
#angles = [1.432,1.909, 2.291, 2.726,3.013,3.576, 4.399, 5.194, 5.711, 6.34, 7.125, 8.13, 8.746, 9.019, 9.462, 9.782, 9.951, 10.3, 10.89, 11.31, 11.77, 12.53,13.09,14.04, 15.12]
angles = [ 2.291, 2.726,3.013,3.576, 4.399, 5.194, 5.711, 6.34, 7.125, 8.13, 8.746, 9.019, 9.462, 9.951]
charge = ['plus']
drt = os.getcwd()
field = "05"
  # Get current working directory
for c in charge:
    ang_index = 0
    for ang in angles:
        data_momentum = open("{0}/data_mu_{1}_momentum_{3}.dat".format(drt,c,field,ang), "w") 
        data_position = open("{0}/data_mu_{1}_position_{3}.dat".format(drt,c, field,ang), "w") 
        for i in range(1,11):
             make_calculus(ang, i,c,drt,ang_index,data_momentum, data_position,field)
        ang_index +=1
    data_momentum.close()
    data_position.close()  

