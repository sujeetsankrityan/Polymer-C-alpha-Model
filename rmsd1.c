/****************************************************************************/
/* rmsd calculation                                                         */
/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NMAX 500 

double center_of_mass(double x[],double y[],double z[],int n);
double correlation(double x1[],double y1[],double z1[],
		   double x2[],double y2[],double z2[],int n);

double rmsd2(int iflag,double x[],double y[], double z[],int n){
  int i;
  double xt[NMAX],yt[NMAX],zt[NMAX],gyr2;
  static double xnat[NMAX],ynat[NMAX],znat[NMAX],gyr2nat;

  if(iflag<0){
    for(i=0;i<n;i++){
      xnat[i]=x[i];
      ynat[i]=y[i];
      znat[i]=z[i];
    }
    gyr2nat=center_of_mass(xnat,ynat,znat,n);
    return 0;
  }

  for(i=0;i<n;i++){
    xt[i]=x[i];
    yt[i]=y[i];
    zt[i]=z[i];
  }
  gyr2=center_of_mass(xt,yt,zt,n);
  return gyr2nat+gyr2-correlation(xnat,ynat,znat,xt,yt,zt,n);
}
/****************************************************************************/
double center_of_mass(double x[],double y[],double z[], int n)
{
  double xcm=0,ycm=0,zcm=0,gyr2=0;
  int j;
  for (j=0;j<n;j++) {
    xcm+=x[j]; 
    ycm+=y[j]; 
    zcm+=z[j];
  }
  xcm*=1.0/n; ycm*=1.0/n; zcm*=1.0/n;
  
  for (j=0;j<n;j++){
    x[j]-=xcm; 
    y[j]-=ycm; 
    z[j]-=zcm;
    gyr2+=x[j]*x[j]+y[j]*y[j]+z[j]*z[j];
  }
  return gyr2/n;
}
/****************************************************************************/
double correlation(double x1[],double y1[],double z1[],
		   double x2[],double y2[],double z2[],int n)
{
  double R[3][3],RtR[3][3];
  double a,b,c,q,r,theta,detR,w[3];
  double pi2=2*acos(-1.);
  int i,j,im;
 
  for (i=0;i<3;i++) for (j=0;j<3;j++) R[i][j]=0;
  for (i=0;i<n;i++) {
    R[0][0]+=x1[i]*x2[i]; R[0][1]+=x1[i]*y2[i]; R[0][2]+=x1[i]*z2[i];
    R[1][0]+=y1[i]*x2[i]; R[1][1]+=y1[i]*y2[i]; R[1][2]+=y1[i]*z2[i];
    R[2][0]+=z1[i]*x2[i]; R[2][1]+=z1[i]*y2[i]; R[2][2]+=z1[i]*z2[i];
  }
  detR=+R[0][0]*R[1][1]*R[2][2]
       +R[0][1]*R[1][2]*R[2][0]
       +R[0][2]*R[1][0]*R[2][1]
       -R[2][0]*R[1][1]*R[0][2]
       -R[2][1]*R[1][2]*R[0][0]
       -R[2][2]*R[1][0]*R[0][1];
  if (detR==0) {printf("****** WARNING! *******\ndetR=0\n");}

  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      RtR[i][j]=R[0][i]*R[0][j]+R[1][i]*R[1][j]+R[2][i]*R[2][j];
    }
  }

  a=-(RtR[0][0]+RtR[1][1]+RtR[2][2]);
  b=RtR[0][0]*RtR[1][1]+RtR[0][0]*RtR[2][2]+RtR[1][1]*RtR[2][2]-
    RtR[0][1]*RtR[0][1]-RtR[0][2]*RtR[0][2]-RtR[1][2]*RtR[1][2];
  c=-detR*detR;

  q=(a*a-3*b)/9;
  r=(2*a*a*a-9*a*b+27*c)/54;
  if (r*r/(q*q*q)>=1) {
    printf("error r^2=%e q^3=%e %e\n",r*r,q*q*q,r*r/(q*q*q));
    exit(-1);
  }
  q=sqrt(q);
  theta=acos(r/(q*q*q));
  w[0]=sqrt(fabs(-2*q*cos(theta/3)-a/3));
  w[1]=sqrt(fabs(-2*q*cos((theta+pi2)/3)-a/3));
  w[2]=sqrt(fabs(-2*q*cos((theta-pi2)/3)-a/3));
  
  if (detR<0) {
    im = w[1] > w[0] ? 1 : 0;
    return 2.*(w[im]+fabs(w[(im+1)%3]-w[(im+2)%3]))/n;
  }else
    return 2.*(w[0]+w[1]+w[2])/n;
}

#undef NMAX
