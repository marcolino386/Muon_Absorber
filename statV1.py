import numpy as np
import os
import glob
import matplotlib.pyplot as plt

def end_calculus(X0,Y0,E0,px0,py0,pz0):
    print("cabo")


def modify_array():
    X0 = X[0]; Y0 = Y[0]; Z0 = Z[0];  E0 = E[0]; px0 = px[0]; py0 = py[0];  pz0 = pz[0]
    X.remove(0); Y.remove(0); Z0.remove(0); E.remove(0); px.remove(0); py.remove(0); pz.remove(0)

def make_calculus(path):
    X0,Y0,E0,px0,py0,pz0 = 0.0, 0.0, 0.0 ,0.0,0.0,0.0
    #deleta a / do argumento path
    data_name = path.replace("/", "")
    #abre dois arquivos para momentum e posição
    data_momentum = open(data_name + "_momentum.dat", "w") 
    data_position = open(data_name + "_position.dat", "w") 
    #entra em todos os arquivos .dat da pasta da variavel path
    for filename in glob.glob(os.path.join(path, '*.dat')):
        #transforma o arquivo em um array
        A = np.loadtxt(filename, dtype = float)
        #se o numero de elementos for 7 significa que nenhum muon passou, apenas há o muon do evento sem absorber
        if np.size(A) == 7:
           #seleciona a variavel de sua respectiva coluna da matriz(array)
           X,Y,E,px,py,pz = A[1],A[2],A[3],A[4],A[5], A[6]


           #Printa

           print("----------------------------------- " +str(E)+ "GEV--------------------------------------------------------------------------------")
           print("E0, X0, Y0, Xmed, Ymed, ENmed ,dx, dy, dEn, desX, desY,desen, final0")
           print(str(E) + " " + str(X) + " " + str(Y) + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " "+ "0" +"\n")
           print("######################################")
           print("en[0], px[0], py[0],theta0, pt0 , pxmed, pymed, ptmed , thetamed , dpx, dpy, dpt, dtheta, despx, despy,slopeX0,slopeY0,slopeXmed,slopeYmed,dslopeX,dslopeY, final0")
           print(str(E) + " " + str(px) + " " + str(py) +  " " + str(np.arctan2(py,px)) + " " + str(np.sqrt(px*px + py*py + pz*pz)) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(px/pz) + " " + str(py/pz) + " " + str(0) + " " + str(0) + " " +  str(0) + " " + str(0) + " " + str(0) + " " + "\n")

           data_position.write(str(E) + " " + str(X) + " " + str(Y) + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + " " + "0" + "\n")


           data_momentum.write(str(E) + " " + str(px) + " " + str(py) +  " " + str(np.arctan2(py,px)) + " " + str(np.sqrt(px*px + py*py + pz*pz)) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(0) + " " + str(px/pz) + " " + str(py/pz) + " " + str(0) + " " + str(0) + " " +  str(0) + " " + str(0) + " " + str(0) + "\n")
           
        else:
           #seleciona a variação de sua respectiva coluna da matriz(array)
           X,Y,E,px,py,pz = A[:,1],A[:,2],A[:,3],A[:,4],A[:,5],A[:,6]
           #modify_array()
           X0 = X[0]; Y0 = Y[0];   E0 = E[0]; px0 = px[0]; py0 = py[0];  pz0 = pz[0]
           X = np.delete(X,0); Y = np.delete(Y,0); E = np.delete(E,0); px = np.delete(px,0); py = np.delete(py,0); pz = np.delete(pz,0);

           #calcula media e desvio padrão utilizando a bilbioteca
           Xmed = np.mean(X); Ymed = np.mean(Y); Emed = np.mean(E); pxmed = np.mean(px); pymed = np.mean(py); pzmed = np.mean(pz)
           desX = np.std(X); desY = np.std(Y); desE = np.std(E); despx = np.std(px); despy = np.std(py); despz = np.std(pz)  

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
           ptmed = np.sqrt(pxmed*pxmed + pymed*pymed + pzmed*pzmed)
           pt0 = np.sqrt(px0*px0 + py0*py0 + pz0*pz0)
           dpt = ptmed - pt0

           #angulo do momentum 
           theta0 = np.arctan2(py0,px0)
           thetamed = np.arctan2(pymed,pxmed)
           dtheta = thetamed - theta0
         
           
           #Printa dados

           # E0, X0, Y0, Xmed, Ymed, ENmed ,dx, dy, dEn, desX, des , desx, desy,desen, final0
           final = np.size(E)
           print("----------------------------------- " +str(E0)+ "GEV--------------------------------------------------------------------------------")
           print("E0, X0, Y0, Xmed, Ymed, ENmed ,dx, dy, dEn, desX, desY,desen, final0")
           print(str(E0) + " " + str(X0) + " " + str(Y0) + " " + str(Xmed) + " " + str(Ymed)  + " " + str(Emed) + " " + str(dx) + " " + str(dy) + " " + str(dE) + " " + str(desX) + " " + str(desY) + " " + str(desE) + " " + str(final) + "\n") 
           print("######################################")
           #en[0], px[0], py[0],theta0, pt0 , pxmed, pymed, ptmed , thetamed , dpx, dpy, dpt, dtheta, despx, despy,slopeX0,slopeY0,slopeXmed,slopeYmed,dslopeX,dslopeY, final0
           print("en[0], px[0], py[0],theta0, pt0 , pxmed, pymed, ptmed , thetamed , dpx, dpy, dpt, dtheta, despx, despy,slopeX0,slopeY0,slopeXmed,slopeYmed,dslopeX,dslopeY, final0")
           print(str(E0) + " " + str(px0) + " " + str(py0) +  " " + str(theta0) + " " + str(pt0) + " " + str(pxmed) + " " + str(pymed) + " " + str(ptmed) + " " + str(thetamed) + " " + str(dpx) + " " + str(dpy) + " " + str(dpt) + " " + str(dtheta) + " " + str(despx) + " " + str(despy) + " " + str(slopeX0) + " " + str(slopeY0) + " " + str(slopeXmed) + " " + str(slopeYmed) + " " + str(dslopeX) + " " + str(dslopeY) + " " + str(final) + "\n") 
           
          


           #escreve nos arquivos de saída

           # "E0, X0, Y0, Xmed, Ymed, ENmed ,dx, dy, dEn, desX, desY,desen, final0"

           data_position.write(str(E0) + " " + str(X0) + " " + str(Y0) + " " + str(Ymed) + " " + str(Emed) + " " + str(dx) + " " + str(dy) + " " + str(dE) + " " + str(desX) + " " + str(desY) + " " + str(desE) + " " + str(final) + "\n")

           #"en[0], px[0], py[0],theta0, pt0 , pxmed, pymed, ptmed , thetamed , dpx, dpy, dpt, dtheta, despx, despy,slopeX0,slopeY0,slopeXmed,slopeYmed,dslopeX,dslopeY, final0"

           data_momentum.write(str(E0) + " " + str(px0) + " " + str(py0) +  " " + str(theta0) + " " + str(pt0) + " " + str(pxmed) + " " + str(pymed) + " " + str(ptmed) + " " + str(thetamed) + " " + str(dpx) + " " + str(dpy) + " " + str(dpt) + " " + str(dtheta) + " " + str(despx) + " " + str(despy) + " " + str(slopeX0) + " " + str(slopeY0) + " " + str(slopeXmed) + " " + str(slopeYmed)+ " " + str(dslopeX) + " " + str(dslopeY) + " " + str(final) + "\n")


    print("DONE")   
    data_momentum.close()
    data_position.close()         
       
   
make_calculus("data_mu_plus2/")
   
make_calculus("data_mu_minus2/")
