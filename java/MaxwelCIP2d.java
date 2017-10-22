import java.awt.*;
import java.awt.event.*;


public class MaxwelCIP2d extends Frame{
	private double [][] EX,EY,EZ,BX,BY,BZ;
	private double [][] EXN,EYN,EZN,BXN,BYN,BZN;
	private double [][] gxEX,gxEY,gxEZ,gxBX,gxBY,gxBZ;
	private double [][] gyEX,gyEY,gyEZ,gyBX,gyBY,gyBZ;
	private double [][] gxEXN,gxEYN,gxEZN,gxBXN,gxBYN,gxBZN;
	private double [][] gyEXN,gyEYN,gyEZN,gyBXN,gyBYN,gyBZN;
	private double [][] out;
        private double c, dt, dx, dy;
        private int imax, jmax ;
         
	public MaxwelCIP2d()	{
                imax=100 ; jmax=100 ;
		EX = new double [101][101];
		EY = new double [101][101];
		EZ = new double [101][101];
		BX = new double [101][101];
		BY = new double [101][101];
		BZ = new double [101][101];
		EXN= new double [101][101];
		EYN= new double [101][101];
		EZN= new double [101][101];
		BXN= new double [101][101];
		BYN= new double [101][101];
		BZN= new double [101][101];
		gxEX = new double [101][101];
		gxEY = new double [101][101];
		gxEZ = new double [101][101];
		gxBX = new double [101][101];
		gxBY = new double [101][101];
		gxBZ = new double [101][101];
		gxEXN= new double [101][101];
		gxEYN= new double [101][101];
		gxEZN= new double [101][101];
		gxBXN= new double [101][101];
		gxBYN= new double [101][101];
		gxBZN= new double [101][101];
		gyEX = new double [101][101];
		gyEY = new double [101][101];
		gyEZ = new double [101][101];
		gyBX = new double [101][101];
		gyBY = new double [101][101];
		gyBZ = new double [101][101];
		gyEXN= new double [101][101];
		gyEYN= new double [101][101];
		gyEZN= new double [101][101];
		gyBXN= new double [101][101];
		gyBYN= new double [101][101];
		gyBZN= new double [101][101];
		out= new double [101][101];

		setSize(400,300);
		addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){
				System.exit(0);
			}
		});
	}

        public void keisan( ) {
           int i,j=0;
           double apl, ami;
           double[] EXp=new double[4];
           double[] EXm=new double[4];
           double[] EYp=new double[4];
           double[] EYm=new double[4];
           double[] EZp=new double[4];
           double[] EZm=new double[4];
           double[] BXp=new double[4];
           double[] BXm=new double[4];
           double[] BYp=new double[4];
           double[] BYm=new double[4];
           double[] BZp=new double[4];
           double[] BZm=new double[4];

           for(j=1;j<jmax;j++) {
           for(i=1;i<imax;i++) {

             CIP1dm(EY[i][j], EY[i-1][j], gxEY[i][j], gxEY[i-1][j], gyEY[i][j],
                         gyEY[i-1][j],-dt*c,-dx,EYp, 1) ; 
             CIP1dm(EY[i][j], EY[i+1][j], gxEY[i][j], gxEY[i+1][j], gyEY[i][j],
                         gyEY[i+1][j], dt*c, dx,EYm, 1) ; 
             CIP1dm(EZ[i][j], EZ[i-1][j], gxEZ[i][j], gxEZ[i-1][j], gyEZ[i][j],
                         gyEZ[i-1][j],-dt*c,-dx,EZp, 1) ; 
             CIP1dm(EZ[i][j], EZ[i+1][j], gxEZ[i][j], gxEZ[i+1][j], gyEZ[i][j],
                         gyEZ[i+1][j], dt*c, dx,EZm, 1) ; 
             CIP1dm(BY[i][j], BY[i-1][j], gxBY[i][j], gxBY[i-1][j], gyBY[i][j],
                         gyBY[i-1][j],-dt*c,-dx,BYp, 1) ; 
             CIP1dm(BY[i][j], BY[i+1][j], gxBY[i][j], gxBY[i+1][j], gyBY[i][j],
                         gyBY[i+1][j],  dt*c, dx,BYm, 1) ; 
             CIP1dm(BZ[i][j], BZ[i-1][j], gxBZ[i][j], gxBZ[i-1][j], gyBZ[i][j],
                         gyBZ[i-1][j],-dt*c,-dx,BZp, 1) ; 
             CIP1dm(BZ[i][j], BZ[i+1][j], gxBZ[i][j], gxBZ[i+1][j], gyBZ[i][j],
                         gyBZ[i+1][j], dt*c, dx,BZm, 1) ; 

             apl=1.0;   ami =1.0;
             EYN[i][j]=(EYp[1]/apl+EYm[1]/ami+BZp[1]-BZm[1])/
                          (1.0/apl+1.0/ami) ;
             gxEYN[i][j]=(EYp[2]/apl+EYm[2]/ami+BZp[2]-BZm[2])/
                          (1.0/apl+1.0/ami) ;
             gyEYN[i][j]=(EYp[3]/apl+EYm[3]/ami+BZp[3]-BZm[3])/
                          (1.0/apl+1.0/ami) ;

             BZN[i][j]=(BZp[1]*apl+BZm[1]*ami+EYp[1]-EYm[1])/(apl+ami)   ;
             gxBZN[i][j]=(BZp[2]*apl+BZm[2]*ami+EYp[2]-EYm[2])/(apl+ami)   ;
             gyBZN[i][j]=(BZp[3]*apl+BZm[3]*ami+EYp[3]-EYm[3])/(apl+ami)   ;

             EZN[i][j]=(EZp[1]/apl+EZm[1]/ami-BYp[1]+BYm[1])/
                          (1.0/apl+1.0/ami)  ;
             gxEZN[i][j]=(EZp[2]/apl+EZm[2]/ami-BYp[2]+BYm[2])/
                          (1.0/apl+1.0/ami)  ;
             gyEZN[i][j]=(EZp[3]/apl+EZm[3]/ami-BYp[3]+BYm[3])/
                          (1.0/apl+1.0/ami)  ;

             BYN[i][j]=(BYp[1]*apl+BYm[1]*ami-EZp[1]+EZm[1])/(apl+ami)   ;
             gxBYN[i][j]=(BYp[2]*apl+BYm[2]*ami-EZp[2]+EZm[2])/(apl+ami)   ;
             gyBYN[i][j]=(BYp[3]*apl+BYm[3]*ami-EZp[3]+EZm[3])/(apl+ami)   ;
             
            }}
 
           for(j=1;j<jmax;j++) {
           for(i=1;i<imax;i++) {

             CIP1dm(EX[i][j], EX[i][j-1], gyEX[i][j], gyEX[i][j-1], gxEX[i][j],
                    gxEX[i][j-1],-dt*c,-dy,EXp , 2) ; 
             CIP1dm(EX[i][j], EX[i][j+1], gyEX[i][j], gyEX[i][j+1], gxEX[i][j],
                    gxEX[i][j+1], dt*c, dy,EXm , 2) ; 
             CIP1dm(EZN[i][j], EZN[i][j-1], gyEZN[i][j], gyEZN[i][j-1], 
                    gxEZN[i][j], gxEZN[i][j-1], -dt*c,-dy,EZp, 2) ; 
             CIP1dm(EZN[i][j], EZN[i][j+1], gyEZN[i][j], gyEZN[i][j+1], 
                    gxEZN[i][j], gxEZN[i][j+1], dt*c, dy,EZm, 2) ; 
             CIP1dm(BX[i][j], BX[i][j-1], gyBX[i][j], gyBX[i][j-1], 
                    gxBX[i][j], gxBX[i][j-1],-dt*c,-dy,BXp, 2) ; 
             CIP1dm(BX[i][j], BX[i][j+1], gyBX[i][j], gyBX[i][j+1], 
                    gxBX[i][j], gxBX[i][j+1], dt*c, dy,BXm, 2) ; 
             CIP1dm(BZN[i][j], BZN[i][j-1], gyBZN[i][j], gyBZN[i][j-1], 
                    gxBZN[i][j], gxBZN[i][j-1], -dt*c,-dy,BZp, 2) ; 
             CIP1dm(BZN[i][j], BZN[i][j+1], gyBZN[i][j], gyBZN[i][j+1], 
                    gxBZN[i][j], gxBZN[i][j+1], dt*c, dy,BZm, 2) ; 

             apl=1.0;   ami =1.0;

	     EZ[i][j]=(EZp[1]/apl+EZm[1]/ami+BXp[1]-BXm[1])/(1.0/apl+1.0/ami) ;
	     gxEZ[i][j]=(EZp[2]/apl+EZm[2]/ami+BXp[2]-BXm[2])/
                        (1.0/apl+1.0/ami) ;
	     gyEZ[i][j]=(EZp[3]/apl+EZm[3]/ami+BXp[3]-BXm[3])/
                        (1.0/apl+1.0/ami) ;

	     BXN[i][j]=(BXp[1]*apl+BXm[1]*ami+EZp[1]-EZm[1])/(apl+ami)  ;
	     gxBXN[i][j]=(BXp[2]*apl+BXm[2]*ami+EZp[2]-EZm[2])/(apl+ami)  ;
	     gyBXN[i][j]=(BXp[3]*apl+BXm[3]*ami+EZp[3]-EZm[3])/(apl+ami)  ;
 
             BZ[i][j]=(BZp[1]*apl+BZm[1]*ami-EXp[1]+EXm[1])/(apl+ami)  ;
             gxBZ[i][j]=(BZp[2]*apl+BZm[2]*ami-EXp[2]+EXm[2])/(apl+ami)  ;
             gyBZ[i][j]=(BZp[3]*apl+BZm[3]*ami-EXp[3]+EXm[3])/(apl+ami)  ;

             EXN[i][j]=(EXp[1]/apl+EXm[1]/ami-BZp[1]+BZm[1])/
                       (1.0/apl+1.0/ami)  ;
             gxEXN[i][j]=(EXp[2]/apl+EXm[2]/ami-BZp[2]+BZm[2])/
                       (1.0/apl+1.0/ami)  ;
             gyEXN[i][j]=(EXp[3]/apl+EXm[3]/ami-BZp[3]+BZm[3])/
                       (1.0/apl+1.0/ami)  ;

           }}

             shift(EX,gxEX,gyEX,EXN,gxEXN,gyEXN);
             shift(BX,gxBX,gyBX,BXN,gxBXN,gyBXN);
             shift(EY,gxEY,gyEY,EYN,gxEYN,gyEYN);
             shift(BY,gxBY,gyBY,BYN,gxBYN,gyBYN);
//           shift(EZ,gxEZ,gyEZ,EZN,gxEZN,gyEZN);
//           shift(BZ,gxBZ,gyBZ,BZN,gxBZN,gyBZN);

          for(i=1;i<=1000000;i++)  j=i+j*j;

 }

 public void shift(double [][] f,double [][] gxf,double [][]gyf,
                         double [][] fn, double [][] gxfn,double [][] gyfn) 
{
    int i,j ;

           for(j=1;j<jmax;j++) {
           for(i=1;i<imax;i++) {
              f[i][j]=fn[i][j];
              gxf[i][j]=gxfn[i][j];
              gyf[i][j]=gyfn[i][j];

           }}
}
  
  public void CIP1dm(double fi, double fiup, double gi, double giup,
                          double hi, double hiup, double gzi, double D, 
                          double [] val, int s)
  {
          double a,b;

             a=(gi+giup)/(D*D)+2.*(fi-fiup)/(D*D*D);
             b=3.*(fiup-fi)/(D*D)-(2.*gi+giup)/D;
             val[1]=((a*gzi+b)*gzi+gi)*gzi+fi ;
           if(s==1) {
             val[2]=(3.*a*gzi+2.*b)*gzi+gi;
             val[3]=hi+(hi-hiup)*gzi/Math.abs(D) ; }
           else {
             val[3]=(3.*a*gzi+2.*b)*gzi+gi;
             val[2]=hi+(hi-hiup)*gzi/Math.abs(D) ; }
           
  }

   public void init(){
           int i,j;
           double radius;
           double sigma=25.0;
            
           c=1.0 ; dx=1.0 ; dy=1.0 ; dt= 0.2 ;

             for(j=0;j<=jmax;j++) { 
             for(i=0;i<=imax;i++) { 
                radius=(j-jmax/2)*(j-jmax/2)+(i-imax/2)*(i-imax/2) ;
                if(radius<4.0) BZ[i][j]=1.0 ;
                else BZ[i][j]=Math.exp( -(radius-4.0)/sigma );
                EX[i][j]=0.0;  EY[i][j]=0.0;  EZ[i][j]=0.0; 
                gxEX[i][j]=0.0; gxEY[i][j]=0.0; gxEZ[i][j]=0.0;
                gyEX[i][j]=0.0; gyEY[i][j]=0.0; gyEZ[i][j]=0.0;
                BY[i][j]=0.0;  BX[i][j]=0.0; 
                gxBX[i][j]=0.0; gxBY[i][j]=0.0; gxBZ[i][j]=0.0;
                gyBX[i][j]=0.0; gyBY[i][j]=0.0; gyBZ[i][j]=0.0;
             }}
   }

	public void dataset( )
	{
          int i,j;
          double fmax, fmin ;

            fmax=0.0; fmin=1.e+30;
          for(j=0;j<=jmax;j++) {
          for(i=0;i<=imax;i++) {
            if(BZ[i][j]>=fmax) fmax=BZ[i][j];
            if(BZ[i][j]<=fmin) fmin=BZ[i][j];
          }}
          for(j=0;j<=jmax;j++) {
          for(i=0;i<=imax;i++) {
            out[i][j]=(BZ[i][j]-fmin)/(fmax-fmin);
           
          }}
 
        }

    private double[] x = new double[7];
    private double[] y = new double[7];
    private double[] d = new double[7];

public void paint(Graphics g)
{
    int i,j,jj;
    int k,icount,iab,iae,ibe,iysize=300,incont=10;
    double dab,xab=0.0,yab=0.0,dae,xae=0.0,yae=0.0,dbe,xbe=0.0,ybe=0.0;
    double xcent,ycent,f0;
    double x1,y1,x2,y2;
    double deltax=2.0, deltay=2.0 ;
 
    dataset();

    for(jj=1;jj<=incont;jj=jj+1) {
       f0=1./(double)incont*(double)jj;

      for(j=0;j<jmax;j++) {
      for(i=0;i<imax;i++) {

	x1=i*deltax; x2=(i+1)*deltax;
	y1=j*deltay; y2=(j+1)*deltay;   
	xcent=(x1+x2)/2.0;  ycent=(y1+y2)/2.0;

	x[1]=x1; x[2]=x2; x[3]=x2; x[4]=x1; x[5]=x[1]; x[6]=xcent;
	y[1]=y1; y[2]=y1; y[3]=y2; y[4]=y2; y[5]=y[1]; y[6]=ycent;
	
	d[1]=out[i][j]; d[2]=out[i+1][j]; d[3]=out[i+1][j+1];
	d[4]=out[i][j+1]; d[5]=d[1];
	d[6]=(d[1]+d[2]+d[3]+d[4])/4.0;

        if(d[1]!=d[2] || d[2]!=d[3] || d[3]!=d[4]) {

          for(k=1;k<=4;k++) {
                   icount=0; iae=0; iab=0; ibe=0;

	        if(d[k+1]!=d[k]) {        
		  dab=(f0-d[k+1])/(d[k]-d[k+1]);
		  if(0.0<=dab && dab<=1.0) {
	               xab=(1.0-dab)*x[k+1]+dab*x[k];
	               yab=(1.0-dab)*y[k+1]+dab*y[k];
	               icount=icount+1; iab=1;
		  }
	        } 
	   if(d[k]!=d[6]) {              
		  dae=(f0-d[k])/(d[6]-d[k]);  
		  if(0.0<=dae && dae<=1.0) {
	               xae=(1.0-dae)*x[k]+dae*x[6];
	               yae=(1.0-dae)*y[k]+dae*y[6];
	               icount=icount+1; iae=1;
		  }
		}

       if(icount!=0) {   
         if(iab==1 && iae==1)  
          g.drawLine((int)xae,iysize-(int)yae,(int)xab,iysize-(int)yab);

		if(d[6]!=d[k+1]) { 
                   dbe=(f0-d[6])/(d[k+1]-d[6]); 
		   if(0.0<=dbe && dbe<=1.0) {
	               xbe=(1.0-dbe)*x[6]+dbe*x[k+1];
	               ybe=(1.0-dbe)*y[6]+dbe*y[k+1];
		       ibe=1;
		   }
	        }	          
		    
	if(ibe==1 && iae==1)
         g.drawLine((int)xbe,iysize-(int)ybe,(int)xae,iysize-(int)yae);
	if(ibe==1 && iab==1) 
         g.drawLine((int)xab,iysize-(int)yab,(int)xbe,iysize-(int)ybe);

        if(d[k]==d[k+1] && d[k]==f0)
         g.drawLine((int)x[k],iysize-(int)y[k],(int)x[k+1],iysize-(int)y[k+1]);

        if(d[k+1]==d[6] && d[6]==f0)
         g.drawLine((int)x[6],iysize-(int)y[6],(int)x[k+1],iysize-(int)y[k+1]);
		    
        if(d[6]==d[k] && d[6]==f0)
         g.drawLine((int)x[k],iysize-(int)y[k],(int)x[6],iysize-(int)y[6]);

            }    
          }  
        }

      } 
      }
    }
}
	
  public static void main(String[] args)
	{
          int kk;
	     MaxwelCIP2d w = new MaxwelCIP2d();
		w.show();
                w.init();
                for(kk=0;kk<=160;kk++) {
                  w.keisan();
                  w.repaint();

	        }

}
}
