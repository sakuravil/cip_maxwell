/* 計算領域の最大グリッド数 */
#define X_m_max 101 
#define Y_m_max 101

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

double EX[X_m_max][Y_m_max];
double EY[X_m_max][Y_m_max];
double EZ[X_m_max][Y_m_max];
double BX[X_m_max][Y_m_max];
double BY[X_m_max][Y_m_max];
double BZ[X_m_max][Y_m_max];
double EXN[X_m_max][Y_m_max];
double EYN[X_m_max][Y_m_max];
double EZN[X_m_max][Y_m_max];
double BXN[X_m_max][Y_m_max];
double BYN[X_m_max][Y_m_max];
double BZN[X_m_max][Y_m_max];
double gxEX[X_m_max][Y_m_max];
double gxEY[X_m_max][Y_m_max];
double gxEZ[X_m_max][Y_m_max];
double gxBX[X_m_max][Y_m_max];
double gxBY[X_m_max][Y_m_max];
double gxBZ[X_m_max][Y_m_max];
double gxEXN[X_m_max][Y_m_max];
double gxEYN[X_m_max][Y_m_max];
double gxEZN[X_m_max][Y_m_max];
double gxBXN[X_m_max][Y_m_max];
double gxBYN[X_m_max][Y_m_max];
double gxBZN[X_m_max][Y_m_max];
double gyEX[X_m_max][Y_m_max];
double gyEY[X_m_max][Y_m_max];
double gyEZ[X_m_max][Y_m_max];
double gyBX[X_m_max][Y_m_max];
double gyBY[X_m_max][Y_m_max];
double gyBZ[X_m_max][Y_m_max];
double gyEXN[X_m_max][Y_m_max];
double gyEYN[X_m_max][Y_m_max];
double gyEZN[X_m_max][Y_m_max];
double gyBXN[X_m_max][Y_m_max];
double gyBYN[X_m_max][Y_m_max];
double gyBZN[X_m_max][Y_m_max];
double out[X_m_max][Y_m_max];

double c, dt, dx, dy;
int imax=X_m_max-1, jmax=Y_m_max-1;

void init(){
	int i, j;
	double radius;
	double sigma = 25.0;

	c = 1.0;
	dx = 1.0;
	dy = 1.0;
	dt = 0.1;

	for(j = 0; j <= jmax; j++){
		for(i = 0; i <= imax; i++){
			
			radius=(j-jmax/2)*(j-jmax/2)+(i-imax/2)*(i-imax/2) ;
			
			/* if(radius<4.0) BZ[i][j]=1.0 ; */
			BZ[i][j] = exp( (-(radius-4.0)/sigma) );
			EX[i][j]=0.0;  EY[i][j]=0.0;  EZ[i][j]=0.0; 
			gxEX[i][j]=0.0; gxEY[i][j]=0.0; gxEZ[i][j]=0.0;
			gyEX[i][j]=0.0; gyEY[i][j]=0.0; gyEZ[i][j]=0.0;
			BY[i][j]=0.0;  BX[i][j]=0.0; 
			gxBX[i][j]=0.0; gxBY[i][j]=0.0; gxBZ[i][j]=0.0;
			gyBX[i][j]=0.0; gyBY[i][j]=0.0; gyBZ[i][j]=0.0;
		}
	}
}


void CIP1dm(double fi, double fiup, double gi, double giup, double hi, double hiup, double gzi, double D, double *val, int s){
	
	double a, b;



			a=(gi+giup)/(D*D)+2.*(fi-fiup)/(D*D*D);
			b=3.*(fiup-fi)/(D*D)-(2.*gi+giup)/D;
			val[1]=((a*gzi+b)*gzi+gi)*gzi+fi ;
		if(s==1) {
			val[2]=(3.*a*gzi+2.*b)*gzi+gi;
			val[3]=hi-(hi-hiup)*gzi/D ;
		}
		else {
			val[3]=(3.*a*gzi+2.*b)*gzi+gi;
			val[2]=hi-(hi-hiup)*gzi/D ;
		}
}

void shift(double f[][Y_m_max], double gxf[][Y_m_max], double gyf[][Y_m_max], double fn[][Y_m_max], double gxfn[][Y_m_max], double  gyfn[][Y_m_max]){
	int i, j;

		for(j = 1; j < jmax; j++){
			for(i = 1; i < imax; i++){
				f[i][j] = fn[i][j];
				gxf[i][j]= gxfn[i][j];
				gyf[i][j] = gyfn[i][j];
			} 
		}
}


void Calc(){
	int i, j;
	double apl, ami;
	double EXp[4];
	double EXm[4];
	double EYp[4];
	double EYm[4];
	double EZp[4];
	double EZm[4];
	double BXp[4];
	double BXm[4];
	double BYp[4];
	double BYm[4];
	double BZp[4];
	double BZm[4];

	for(j = 1; j < jmax; j++){
		for(i = 1; i < imax; i++){
			CIP1dm(EY[i][j], EY[i-1][j], gxEY[i][j], gxEY[i-1][j], gyEY[i][j], gyEY[i-1][j],-dt*c,-dx,EYp, 1); 
			CIP1dm(EY[i][j], EY[i+1][j], gxEY[i][j], gxEY[i+1][j], gyEY[i][j], gyEY[i+1][j], dt*c, dx,EYm, 1); 
			CIP1dm(EZ[i][j], EZ[i-1][j], gxEZ[i][j], gxEZ[i-1][j], gyEZ[i][j], gyEZ[i-1][j],-dt*c,-dx,EZp, 1); 
			CIP1dm(EZ[i][j], EZ[i+1][j], gxEZ[i][j], gxEZ[i+1][j], gyEZ[i][j], gyEZ[i+1][j], dt*c, dx,EZm, 1); 
			CIP1dm(BY[i][j], BY[i-1][j], gxBY[i][j], gxBY[i-1][j], gyBY[i][j], gyBY[i-1][j],-dt*c,-dx,BYp, 1); 
			CIP1dm(BY[i][j], BY[i+1][j], gxBY[i][j], gxBY[i+1][j], gyBY[i][j], gyBY[i+1][j],  dt*c, dx,BYm, 1); 
			CIP1dm(BZ[i][j], BZ[i-1][j], gxBZ[i][j], gxBZ[i-1][j], gyBZ[i][j], gyBZ[i-1][j],-dt*c,-dx,BZp, 1); 
			CIP1dm(BZ[i][j], BZ[i+1][j], gxBZ[i][j], gxBZ[i+1][j], gyBZ[i][j], gyBZ[i+1][j], dt*c, dx,BZm, 1); 

			apl = 1.0;
			ami = 1,0;
			
      EYN[i][j]=(EYp[1]/apl+EYm[1]/ami+BZp[1]-BZm[1])/(1.0/apl+1.0/ami) ;
      gxEYN[i][j]=(EYp[2]/apl+EYm[2]/ami+BZp[2]-BZm[2])/(1.0/apl+1.0/ami) ;
      gyEYN[i][j]=(EYp[3]/apl+EYm[3]/ami+BZp[3]-BZm[3])/(1.0/apl+1.0/ami) ;
		
			BZN[i][j]=(BZp[1]*apl+BZm[1]*ami+EYp[1]-EYm[1])/(apl+ami);
      gxBZN[i][j]=(BZp[2]*apl+BZm[2]*ami+EYp[2]-EYm[2])/(apl+ami);
      gyBZN[i][j]=(BZp[3]*apl+BZm[3]*ami+EYp[3]-EYm[3])/(apl+ami);

      EZN[i][j]=(EZp[1]/apl+EZm[1]/ami-BYp[1]+BYm[1])/(1.0/apl+1.0/ami);
			gxEZN[i][j]=(EZp[2]/apl+EZm[2]/ami-BYp[2]+BYm[2])/(1.0/apl+1.0/ami);
	    gyEZN[i][j]=(EZp[3]/apl+EZm[3]/ami-BYp[3]+BYm[3])/(1.0/apl+1.0/ami);

      BYN[i][j]=(BYp[1]*apl+BYm[1]*ami-EZp[1]+EZm[1])/(apl+ami);
      gxBYN[i][j]=(BYp[2]*apl+BYm[2]*ami-EZp[2]+EZm[2])/(apl+ami);
      gyBYN[i][j]=(BYp[3]*apl+BYm[3]*ami-EZp[3]+EZm[3])/(apl+ami);
		}
	}

	for(j = 1; j < jmax; j++){
		for(i = 1; i < imax; i++){
			
      CIP1dm(EX[i][j], EX[i][j-1], gyEX[i][j], gyEX[i][j-1], gxEX[i][j], gxEX[i][j-1],-dt*c,-dy,EXp , 2) ; 
      CIP1dm(EX[i][j], EX[i][j+1], gyEX[i][j], gyEX[i][j+1], gxEX[i][j], gxEX[i][j+1], dt*c, dy,EXm , 2) ; 
      CIP1dm(EZN[i][j], EZN[i][j-1], gyEZN[i][j], gyEZN[i][j-1], gxEZN[i][j], gxEZN[i][j-1], -dt*c,-dy,EZp, 2) ; 
      CIP1dm(EZN[i][j], EZN[i][j+1], gyEZN[i][j], gyEZN[i][j+1], gxEZN[i][j], gxEZN[i][j+1], dt*c, dy,EZm, 2) ; 
      CIP1dm(BX[i][j], BX[i][j-1], gyBX[i][j], gyBX[i][j-1], gxBX[i][j], gxBX[i][j-1],-dt*c,-dy,BXp, 2) ; 
      CIP1dm(BX[i][j], BX[i][j+1], gyBX[i][j], gyBX[i][j+1], gxBX[i][j], gxBX[i][j+1], dt*c, dy,BXm, 2) ; 
      CIP1dm(BZN[i][j], BZN[i][j-1], gyBZN[i][j], gyBZN[i][j-1], gxBZN[i][j], gxBZN[i][j-1], -dt*c,-dy,BZp, 2) ; 
      CIP1dm(BZN[i][j], BZN[i][j+1], gyBZN[i][j], gyBZN[i][j+1], gxBZN[i][j], gxBZN[i][j+1], dt*c, dy,BZm, 2) ; 

      apl=1.0;
			ami =1.0;

	    EZ[i][j]=(EZp[1]/apl+EZm[1]/ami+BXp[1]-BXm[1])/(1.0/apl+1.0/ami);
	    gxEZ[i][j]=(EZp[2]/apl+EZm[2]/ami+BXp[2]-BXm[2])/(1.0/apl+1.0/ami);
	    gyEZ[i][j]=(EZp[3]/apl+EZm[3]/ami+BXp[3]-BXm[3])/(1.0/apl+1.0/ami) ;

	    BXN[i][j]=(BXp[1]*apl+BXm[1]*ami+EZp[1]-EZm[1])/(apl+ami)  ;
	    gxBXN[i][j]=(BXp[2]*apl+BXm[2]*ami+EZp[2]-EZm[2])/(apl+ami)  ;
	    gyBXN[i][j]=(BXp[3]*apl+BXm[3]*ami+EZp[3]-EZm[3])/(apl+ami)  ;
 
      BZ[i][j]=(BZp[1]*apl+BZm[1]*ami-EXp[1]+EXm[1])/(apl+ami)  ;
      gxBZ[i][j]=(BZp[2]*apl+BZm[2]*ami-EXp[2]+EXm[2])/(apl+ami)  ;
      gyBZ[i][j]=(BZp[3]*apl+BZm[3]*ami-EXp[3]+EXm[3])/(apl+ami)  ;

      EXN[i][j]=(EXp[1]/apl+EXm[1]/ami-BZp[1]+BZm[1])/(1.0/apl+1.0/ami)  ;
      gxEXN[i][j]=(EXp[2]/apl+EXm[2]/ami-BZp[2]+BZm[2])/(1.0/apl+1.0/ami)  ;
      gyEXN[i][j]=(EXp[3]/apl+EXm[3]/ami-BZp[3]+BZm[3])/(1.0/apl+1.0/ami)  ;
		}
	}
			shift(EX,gxEX,gyEX,EXN,gxEXN,gyEXN);
      shift(BX,gxBX,gyBX,BXN,gxBXN,gyBXN);
      shift(EY,gxEY,gyEY,EYN,gxEYN,gyEYN);
			shift(BY,gxBY,gyBY,BYN,gxBYN,gyBYN);

		//	for(i = 1; i < 1000000; i++)
			//	j = i+j*j;
}


void dataset(){
	int i,j;
	double fmax, fmin ;

	fmax=0.0;
	fmin=1.e+30;
	for(j=0;j<=jmax;j++) {
		for(i=0;i<=imax;i++) {
			if(BZ[i][j]>=fmax) fmax=BZ[i][j];
			if(BZ[i][j]<=fmin) fmin=BZ[i][j];
		}
	}
	for(j=0;j<=jmax;j++) {
		for(i=0;i<=imax;i++) {
		out[i][j]=(BZ[i][j]-fmin)/(fmax-fmin);
		}
	}
}

void output(double out[][Y_m_max]){
	int i, j;
	char file[10]="field.dat";
	FILE *fp;
	
	if(NULL == (fp = fopen(file, "w"))){
		printf("\n\n Can not open file : %s\n", file);;
		exit(1);
	}

	for(j = 0; j <= jmax; j++){
		for(i = 0; i <= imax; i++){
			fprintf(fp,"%d %d %f\t\n", i, j, (float)out[i][j]);
		}
	fprintf(fp, "\n");
	}
	fclose(fp);
}

int main(){
	int k;
	init();
	for(int k =0; k <= 160; k++){
		Calc();
		dataset();
	}

		output(out);
}



