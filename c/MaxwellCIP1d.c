#define X_m_max 101
#define Y_m_max 101

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define NZ 201

  float Ey[NZ];
  float Bx[NZ];
  float Eynew[NZ];
  float Bxnew[NZ];
  float gEy[NZ];
  float gBx[NZ];
  float gEynew[NZ];
  float gBxnew[NZ];
  float kai;
  float dt;

  float P1[NZ][NZ];
  float Vx[NZ][NZ];
  float Vy[NZ][NZ];
  float P1new[NZ][NZ];
  float Vxnew[NZ][NZ];
  float Vynew[NZ][NZ];

  float C1=1.0;
  float C2=1.0;
  float C3=1.0;

  void FDTD1d(){
    
    int i, j;

    for(i=0; i<NZ; i++){
      for(j=0; j<=NZ; j++){
        Vxnew[i+1][j]=Vx[i+1][j]-C1*(P1[i-1][j]-P1[i][j]);
      }
    }

    for(i=0; i<=NZ; i++){
      for(j=0; j<NZ; j++){
        Vynew[i][j+1]=Vy[i][j+1]-C1*(P1[i][j+1]-P1[i][j]);
      }
    }

    for(i=0; i<=NZ; i++){
      for(j=0; j<=NZ; j++){
        P1new[i][j]=P1[i][j]-C2*(Vxnew[i+1][j]-Vxnew[i][j])-C3*(Vynew[i][j+1]-Vynew[i][j]);
      }
    }
  }


  /*************** CIP algorithm ***************************/
  void CIP1d(float fi, float fiup, float gi, float giup,
                            float gzi, float D, float *temp )
  {
    float a,b;

    a=(gi+giup)/(D*D)+2.*(fi-fiup)/(D*D*D);
    b=3.*(fiup-fi)/(D*D)-(2.*gi+giup)/D;
    temp[1]=((a*gzi+b)*gzi+gi)*gzi+fi ;
    temp[2]=(3.*a*gzi+2.*b)*gzi+gi;
  }

  void keisan(){
    int i,j=0;
    int ischeme=4;
    float dx=1.0;
    float tempEp[3];
    float tempEm[3];
    float tempBp[3];
    float tempBm[3];

    for(i=1;i<NZ;i++){

    /* CIP1d(Ey[i], Ey[i-1], gEy[i], gEy[i-1], -dx*kai,-dx,tempEp) ;  */
    /* CIP1d(Ey[i], Ey[i+1], gEy[i], gEy[i+1],  dx*kai, dx,tempEm) ;  */
    /* CIP1d(Bx[i], Bx[i-1], gBx[i], gBx[i-1], -dx*kai,-dx,tempBp) ;  */
    /* CIP1d(Bx[i], Bx[i+1], gBx[i], gBx[i+1],  dx*kai, dx,tempBm) ;  */

    CIP1d(Ey[i], Ey[i-1], gEy[i], gEy[i-1], -dt*kai,-dx,tempEp) ; 
    CIP1d(Ey[i], Ey[i+1], gEy[i], gEy[i+1],  dt*kai, dx,tempEm) ; 
    CIP1d(Bx[i], Bx[i-1], gBx[i], gBx[i-1], -dt*kai,-dx,tempBp) ; 
    CIP1d(Bx[i], Bx[i+1], gBx[i], gBx[i+1],  dt*kai, dx,tempBm) ; 



    Bxnew[i]=(tempBm[1]+tempBp[1]+tempEm[1]-tempEp[1] )/2.0 ;
    Eynew[i]=(tempBm[1]-tempBp[1]+tempEm[1]+tempEp[1] )/2.0 ;
    gBxnew[i]=(tempBm[2]+tempBp[2]+tempEm[2]-tempEp[2] )/2.0 ;
    gEynew[i]=(tempBm[2]-tempBp[2]+tempEm[2]+tempEp[2] )/2.0 ;
    }


    /* for(i=1;i<=1000000;i++)  j=i+j*j; */

  }

  void shift_fdtd() {
    int i, j;
    for(i=0;i<NZ;i++) {
      for(j=0; j<NZ; j++){
        P1[i][j]=P1new[i][j];
        Vx[i][j]=Vxnew[i][j];
        Vy[i][j]=Vynew[i][j];
      }
    }
  }

  void shift_cip() {
    int i;
    for(i=0;i<NZ;i++) {
      Ey[i]=Eynew[i];
      gEy[i]=gEynew[i];
      Bx[i]=Bxnew[i];
      gBx[i]=gBxnew[i];
    }
  }
   

  void output_fdtd(){
    int i, j;
    char file[10]="field.dat";
    FILE *fp;
    
    if(NULL == (fp = fopen(file, "w"))){
      printf("\n\n Can not open file : %s\n", file);;
      exit(1);
    }

      for(i=0; i<=NZ; i++){
        for(j=0; j<=NZ; j++){
          fprintf(fp,"%d %f\t\n",i, (float)P1[i][j]);
        }
      fprintf(fp, "\n");
      }
    fclose(fp);
  }


  void output_cip(float *out){
    int i;
    char file[10]="field.dat";
    FILE *fp;
    
    if(NULL == (fp = fopen(file, "w"))){
      printf("\n\n Can not open file : %s\n", file);;
      exit(1);
    }

      for(i=0; i<=NZ; i++){
        fprintf(fp,"%d %f\t\n", i, (float)Ey[i]);
      }
    /* fprintf(fp, "\n"); */
    fclose(fp);
  }

  void init_cip1(){
    int i;
    float zz;
    float sigma=0.1;
    dt=0.5;
    kai=0.2;
    float z=-0.1;
    for(i=0; i<NZ; i++){
      zz=(float)(z-1)/sigma;
      /* Ey[i]=0.0-0.5*exp(-zz*zz);  */
      Ey[i]=0.5*exp(-zz*zz); 
      gEy[i]=0.0;
      Bx[i]=0.0;
      gBx[i]=0.0;
      z+=0.01;
    }
  }


  void init_cip2(){
    int i;
    float zz;
    float sigma=25.0;
    dt=0.0;
    kai=0.2;
    float z=-0.1;
    for(i=0; i<NZ; i++){
      zz=(float)(z-1)/sigma;
      /* Ey[i]=0.0-0.5*exp(-zz*zz);  */
      Ey[i]=0.5*exp(-zz*zz); 
      gEy[i]=0.0;
      Bx[i]=0.0;
      gBx[i]=0.0;
      z+=0.01;
    }
  }

void init_fdtd(){
  int i, j ;
  float z=0.0;

  for(i=0; i<=NZ; i++){
    for(j=0; j<=NZ; j++){
      if(z>=0.4 && z<=0.6)
        P1[i][j]=1.0; 
      else
        P1[i][j]=0.0;
      Vx[i][j]=0.0;
      Vy[i][j]=0.0;
     z+=0.02;
    }
  }

}

/* void init(){ */
/*   int i; */
/*   float zz; */
/*   int WN=4; */
/*   kai=0.2; */
/*   for(i=0; i<NZ; i++){ */
/*     #<{(| Ey[i]=0.0-0.5*exp(-zz*zz);  |)}># */
/*     if(z>=0.4 && z<=0.6) */
/*       Ey[i]=1.0;  */
/*     else */
/*       Ey[i]=0.0; */
/*     gEy[i]=0.0; */
/*     Bx[i]=0.0; */
/*     gBx[i]=0.0; */
/*     z+=0.02; */
/*   } */
/* } */


/************** main *****************************/
int main()
{
  int kk;
  init_cip1();
  /* init_fdtd(); */
  for(kk=0;kk<=160;kk++) {
 // for(kk=0;kk<=3;kk++) {
/* 
  FDTD1d();
  shift_fdtd();
  output_fdtd();
 */   
  keisan();
  shift_cip();
  output_cip(Ey);
	}

}

