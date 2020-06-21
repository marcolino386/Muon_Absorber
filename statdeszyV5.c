
#include <stdio.h> 
#include <dirent.h> 
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>



bool has_txt_extension(char const *name)
{
    size_t len = strlen(name);
    return len > 4 && strcmp(name + len - 4, ".dat") == 0;
}


void calc_and_write(FILE* data,FILE* res1, FILE* res2){




int n,nn,s,final0, fina;



double x[50000], y[50000], en[50000], px[50000], py[50000], pz[50000];
double xsoma, ysoma, ensoma,pxsoma,pysoma,pzsoma,xmed,ymed,enmed, pxmed, pymed, pzmed ,dx,dy,den, dpx, dpy, dpz = 0.0, final;
double desx, vasx, desy, vasy, desen,vasen, despx, vaspx, despy, vaspy, despz, vaspz, ptmed, thetamed,pt0,theta0,vaspt,vastheta,despt,destheta,dpt,dtheta,slopeX0, slopeXmed,slopeY0,slopeYmed,dslopeY,dslopeX = 0.0; 

        

	for(nn=0; !feof(data); nn++) { fscanf(data, "%i %lf %lf %lf %lf %lf %lf",&n, &x[nn], &y[nn], &en[nn], &px[nn], &py[nn], &pz[nn] ); final0=nn;}
	fclose(data);


	xsoma = 0.0, ysoma = 0.0, ensoma = 0.0, pxsoma, pysoma, pzsoma = 0.0;

	for(s=1;s<= final0;s++){
	xsoma = xsoma+x[s];
	ysoma = ysoma+y[s];
	ensoma = ensoma+en[s];
        pxsoma = pxsoma + px[s];
        pysoma = pysoma + py[s];
        pzsoma = pzsoma + pz[s];
	}

	final= (final0-1)/1.0;
        if (final0 == 0) {
           xmed, ymed, enmed,pxmed,pymed,pzmed, ptmed, thetamed, slopeXmed, slopeYmed = 0.0;
	} else {
        xmed=xsoma/final;
	ymed=ysoma/final;
	enmed=ensoma/final;
        pxmed = pxsoma/final;
        pymed = pysoma/final;
        pzmed = pzsoma/final;
        ptmed = sqrt(pxmed*pxmed + pymed*pymed);
        thetamed = atan2(pymed,pxmed);
        
	slopeXmed = pxmed/pzmed;
	
	slopeYmed = pymed/pzmed;       
	}

	

	
for(s=1;s<= final0;s++){
	vasx = (xmed-x[s])*(xmed-x[s])+ vasx;
	vasy = (ymed-y[s])*(ymed-y[s])+ vasy;
	vasen = (enmed-en[s])*(enmed-en[s])+ vasen;
        vaspx = (pxmed-px[s])*(pxmed-px[s])+ vaspx;
        vaspy = (pymed-py[s])*(pymed-py[s])+ vaspy;
        vaspz = (pzmed-pz[s])*(pzmed-pz[s])+ vaspz;
        vaspt = (ptmed - sqrt(px[s]*px[s] + py[s]*py[s]))*(ptmed - sqrt(px[s]*px[s] + py[s]*py[s])) + vaspt;
        vastheta = (thetamed -atan2(py[s],px[s]))*(thetamed -atan2(py[s],px[s]));

	}
  
      if(final0 == 0) {

	desx, vasx, desy, vasy, desen, despx, despy, despz, despt, destheta = 0.0;

} else {
        desx=sqrt(vasx)/final;
	desy=sqrt(vasy)/final;
	desen=sqrt(vasen)/final;
        despx=sqrt(vaspx)/final;
        despy=sqrt(vaspy)/final;
        despz=sqrt(vaspz)/final;
        despt = sqrt(vaspt)/final;
        destheta = sqrt(vastheta)/final;


}
	
slopeX0 = px[0]/pz[0];

slopeY0 = py[0]/pz[0];
       

dx=x[0]-xmed;
dy=y[0]-ymed;
den=en[0]-enmed;
dpx = px[0] - pxmed;
dpy = py[0] - pymed;
dpz = pz[0] - pzmed;
dslopeX = slopeXmed - slopeX0;
dslopeY = slopeYmed - slopeY0;


pt0 = sqrt(px[0]*px[0] + py[0]*py[0]);
theta0 = atan2(py[0],px[0]);

dpt = pt0 - ptmed;
dtheta = theta0 - thetamed;



printf("E0, X0, Y0, Xmed, Ymed, ENmed ,dx, dy, dEn, final \n E0], px[0], py[0], pt0 ,theta0 ,pxmed, pymed, ptmed ,thetamed ,dpx,dpy, dpt ,dtheta,despx,despy, despt, destheta,final");
printf("%f %f %f %f %f %f %f %f %f %f %f %f %i \n %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",en[0],x[0],y[0], xmed,ymed,enmed, dx,dy,den, desx, desy,desen,final0, px[0], py[0], pxmed, pymed, dpx,dpy, despx,despy, pt0, ptmed, dpt, despt, theta0, thetamed, dtheta, destheta);
printf("==================================================== \n");

fprintf(res1, "%f %f %f %f %f %f %f %f %f %f %f %f %i \n",en[0],x[0],y[0], xmed,ymed,enmed, dx,dy,den, desx, desy,desen,final0);
fprintf(res2, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n", en[0], px[0], py[0], pt0 , theta0 , pxmed, pymed, ptmed , thetamed , dpx, dpy, dpt, dtheta, despx, despy, despt, destheta, final ,slopeX0,slopeY0,slopeXmed,slopeYmed,dslopeX,dslopeY);    

printf("==================================================== \n");

//fprintf(res2, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %i \n", en[0], px[0], py[0], pt0 ,theta0 ,pxmed, pymed, ptmed ,thetamed ,dpx,dpy, dpt ,dtheta,despx,despy, despt, destheta,final0); 

}


void make_calculus(char dir_name[]) {

struct dirent *de;  // Pointer for directory entry  
    DIR *dr = opendir(dir_name);
    FILE *data; FILE *res1; FILE *res2;   
    
    if (dr == NULL)  // opendir returns NULL if couldn't open directory 
    { 
        printf("Could not open current directory" ); 
        return;
    } 
     printf("Inside %s directory\n", dir_name);
     printf(" \n", dir_name);
     while ((de = readdir(dr)) != NULL) {
          bool is_dat = has_txt_extension(de->d_name);
         // printf(de->d_name);
          if (de->d_type == DT_REG && is_dat == true) {
             char path_data[100];
             char path_write1[100];
             char path_write2[100];
             sprintf(path_data, "%s/%s" ,dir_name, de->d_name);
             sprintf(path_write1, "%s_position.dat" ,dir_name);
             sprintf(path_write2, "%s_momentum.dat" ,dir_name);
             //printf ("%s/%s\n",dir_name, de->d_name);
             data = fopen(path_data, "r"); 
             fseek(data, 0, SEEK_END);
             if (ftell(data) == 0) {return;}
             fclose(data);
             data = fopen(path_data, "r");
             res1 = fopen(path_write1, "a"); 
             res2 = fopen(path_write2, "a"); 
             if(data == NULL) {return;}
             calc_and_write(data, res1, res2);
     	            
         }
        
           
            
  }
    closedir(dr);  



}



int main (int argc, char *argv[]) {

    make_calculus("data_mu_plus1");
    make_calculus("data_mu_plus2");
    make_calculus("data_mu_minus1");
    make_calculus("data_mu_minus2");
    return 0; 

}