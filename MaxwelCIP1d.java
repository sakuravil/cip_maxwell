import java.awt.*;
import java.awt.event.*;

public class MaxwelCIP1d extends Frame{
	private double [] Ey,Bx;
	private double [] Eynew,Bxnew;
	private double [] gEy,gBx;
	private double [] gEynew,gBxnew;
        private double kai;
         
	public MaxwelCIP1d()	{
		Ey = new double [101];
		Eynew = new double [101];
		gEy = new double [101];
		gEynew = new double [101];
		Bx = new double [101];
		Bxnew = new double [101];
		gBx = new double [101];
		gBxnew = new double [101];
		setSize(400,300);
		addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){
				System.exit(0);
			}
		});
	}

        public void keisan( ) {
           int i,j=0;
           int ischeme=4;
           double dx=1.0;
           double[] tempEp=new double[3];
           double[] tempEm=new double[3];
           double[] tempBp=new double[3];
           double[] tempBm=new double[3];

           for(i=1;i<100;i++) {

             CIP1d(Ey[i], Ey[i-1], gEy[i], gEy[i-1], -dx*kai,-dx,tempEp) ; 
             CIP1d(Ey[i], Ey[i+1], gEy[i], gEy[i+1],  dx*kai, dx,tempEm) ; 
             CIP1d(Bx[i], Bx[i-1], gBx[i], gBx[i-1], -dx*kai,-dx,tempBp) ; 
             CIP1d(Bx[i], Bx[i+1], gBx[i], gBx[i+1],  dx*kai, dx,tempBm) ; 
             Bxnew[i]=( tempBm[1]+tempBp[1]+tempEm[1]-tempEp[1] )/2.0 ;
             Eynew[i]=( tempBm[1]-tempBp[1]+tempEm[1]+tempEp[1] )/2.0 ;
             gBxnew[i]=( tempBm[2]+tempBp[2]+tempEm[2]-tempEp[2] )/2.0 ;
             gEynew[i]=( tempBm[2]-tempBp[2]+tempEm[2]+tempEp[2] )/2.0 ;
           }


          for(i=1;i<=1000000;i++)  j=i+j*j;

           }
        
        public void CIP1d(double fi, double fiup, double gi, double giup,
                          double gzi, double D, double[] temp )
       {
          double a,b;

             a=(gi+giup)/(D*D)+2.*(fi-fiup)/(D*D*D);
             b=3.*(fiup-fi)/(D*D)-(2.*gi+giup)/D;
             temp[1]=((a*gzi+b)*gzi+gi)*gzi+fi ;
             temp[2]=(3.*a*gzi+2.*b)*gzi+gi;
        }

        public void shift( ) {
            int i;
           for(i=0;i<=100;i++) {
              Ey[i]=Eynew[i];
              gEy[i]=gEynew[i];
              Bx[i]=Bxnew[i];
              gBx[i]=gBxnew[i];
            }
        }
 
	public void init(){
           int i;
           double zz;
           double sigma=5.0;
            kai=0.1;
               
             for(i=0;i<=100;i++) { 
                zz=(double)(i-50)/sigma ;
                Ey[i]=0.0-0.5*Math.exp(-zz*zz); 
                gEy[i]=0.0;
                Bx[i]=0.0+0.2*Math.exp(-zz*zz); gBx[i]=0.0;
             }
	}

	public void paint(Graphics g)
	{
          int i,x1,y1,x2,y2,kk;
          int totalw=400,totalh=300,w=6,h=6;

       /*    g.setColor(Color.white);
           g.fillRect(0,0,totalw,totalh);
           g.setColor(Color.black);*/
           g.drawLine(40,totalh-50,totalw-40,totalh-50);
           g.drawLine(40,totalh-20,40,50);
           g.drawString("x",totalw-40,totalh-30);
           g.drawString("f",50,50);

           for(i=0;i<=100;i++){
             x1=50+i*3-w/2;  y1=totalh-50-(int)(Bx[i]*100.)-h/2;
  //           x1=50+i*3-w/2;  y1=totalh-50-(int)(fnew[i]*100.)-h/2;
             g.fillOval(x1,y1, w, h);            
         }
	}

	public static void main(String[] args)
	{
          int kk;
	     MaxwelCIP1d w = new MaxwelCIP1d();
		w.show();
                w.init();
                for(kk=0;kk<=350;kk++) {
       //        for(kk=0;kk<=3;kk++) {
                  w.keisan();
                  w.shift();
                  w.repaint();

	        }

}
}
