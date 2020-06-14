
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


void calc_and_write(FILE* data,FILE* res){




int n,nn,s, final;



double x[5000], y[5000], en[5000], px[5000], py[5000], pz[5000], xsoma, ysoma, ensoma,pxsoma,pysoma,pzsoma,xmed,ymed,enmed, pxmed, pymed, pzmed ,dx,dy,den, dpx, dpy, dpz;
double desx, vasx, desy, vasy, desen,vasen, despx, vaspx, despy, vaspy, despz, vaspz; 

        

	for(nn=0; !feof(data); nn++) { fscanf(data, "%i %lf %lf %lf %lf %lf %lf",&n, &x[nn], &y[nn], &en[nn], &px[nn], &py[nn], &pz[nn] ); final=nn;}
	fclose(data);


	xsoma = 0.0, ysoma = 0.0, ensoma = 0.0;

	for(s=1;s<= final;s++){
	xsoma = xsoma+x[s];
	ysoma = ysoma+y[s];
	ensoma = ensoma+en[s];
        pxsoma = pxsoma + px[s];
        pysoma = pysoma + py[s];
        pzsoma = pzsoma + pz[s];
	}

	final= final-1;
	xmed=xsoma/final;
	ymed=ysoma/final;
	enmed=ensoma/final;
        pxmed = pxsoma/final;
        pymed = pysoma/final;
        pzmed = pzsoma/final;

	vasx=0;
	vasy=0;
	vasen=0;
        
	
for(s=1;s<= final;s++){
	vasx = (xmed-x[s])*(xmed-x[s])+ vasx;
	vasy = (ymed-y[s])*(ymed-y[s])+ vasy;
	vasen = (enmed-en[s])*(enmed-en[s])+ vasen;
        vaspx = (pxmed-px[s])*(pxmed-px[s])+ vaspx;
        vaspy = (pymed-py[s])*(pymed-py[s])+ vaspy;
        vaspz = (pzmed-pz[s])*(pzmed-pz[s])+ vaspz;
	}

	desx=sqrt(vasx)/final;
	desy=sqrt(vasy)/final;
	desen=sqrt(vasen)/final;
        despx=sqrt(vaspx)/final;
        despy=sqrt(vaspy)/final;
        despz=sqrt(vaspz)/final;

dx=x[0]-xmed;
dy=y[0]-ymed;
den=en[0]-enmed;
dpx = px[0] - pxmed;
dpy = py[0] - pymed;
dpz = pz[0] - pzmed;

//printf("final, xmed, dx,ymed,dy,enmed,den\n");
printf("%f %f %f %f %f %f %f %f %f %f %f %f %i %f %f %f %f %f %f %f %f %f %f %f %f\n",en[0],x[0],y[0], xmed,ymed,enmed, dx,dy,den, desx, desy,desen,final, px[0], py[0], pz[0], pxmed, pymed, pzmed, dpx,dpy,dpz, despx,despy,despz);
printf("==================================================== \n");

fprintf(res, "%f %f %f %f %f %f %f %f %f %f %f %f %i %f %f %f %f %f %f %f %f %f %f %f %f\n",en[0],x[0],y[0], xmed,ymed,enmed, dx,dy,den, desx, desy,desen,final, px[0], py[0], pz[0], pxmed, pymed, pzmed, dpx,dpy,dpz, despx,despy,despz);
     
printf("==================================================== \n");


}


void make_calculus(char dir_name[]) {

struct dirent *de;  // Pointer for directory entry  
    DIR *dr = opendir(dir_name);
    FILE *data; FILE *res;
    
    if (dr == NULL)  // opendir returns NULL if couldn't open directory 
    { 
        printf("Could not open current directory" ); 
        return;
    } 
     printf("Inside %s directory\n", dir_name);
     printf(" \n", dir_name);
     while ((de = readdir(dr)) != NULL) {
          bool is_dat = has_txt_extension(de->d_name);
          
          if (de->d_type == DT_REG && is_dat == true) {
             char path_data[100];
             char path_write[100];
             sprintf(path_data, "%s/%s" ,dir_name, de->d_name);
             sprintf(path_write, "%s.dat" ,dir_name);
             printf ("%s/%s\n",dir_name, de->d_name);
             data = fopen(path_data, "r");
            res = fopen(path_write, "a"); 
           if(data == NULL || res == NULL) {return;}
           calc_and_write(data, res);
     	            
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
