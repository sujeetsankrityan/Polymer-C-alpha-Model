# include <stdio.h>
# include <stdlib.h>
# include <math.h>
#define MC_STEP 100

int main() {
  int NTA,i,j;       //Number of Total Atoms
  NTA=MC_STEP;
  
  double theta,theta2,pi,radius;
  double X[NTA],Y[NTA],Z[NTA];
  // double bead_dist=32,bead_rad=16.0;
  double bead_dist=32,bead_rad=16.0;


  FILE *fp,*fp1,*fp2,*fp3;

  fp = fopen("initial_circular_coordinate","w");
  fp1 = fopen("Initial_circular_polymer_pdb_structure.pdb","w");

  pi=acos(-1);
  theta=(2*pi)/MC_STEP;
  theta2=pi/6;
  radius = (NTA)*bead_dist/(2*pi);
  printf("%d %f %f %f\n",NTA,pi,theta,radius);

  for (i=0;i<NTA;i++) {
    X[i]=radius*cos(i*theta);
    Y[i]=radius*sin(i*theta);//*cos(theta2);
    Z[i]=0;
    //    Z[i]=radius*sin(i*theta)*sin(theta2);
    //printf("%f %f %f %f\n",i*theta,X[i],Y[i],Z[i]);
 
 }


  for (j=0;j<NTA;j++) {//fprintf(fp1,"%s  %4d %s %s %s %4d    %8.1f %8.1f %8.1f   0.00   1.00\n","ATOM",j,"CA","AAA","A",j,X[j],Y[j],Z[j]);
fprintf(fp1,"ATOM   %4u  CA AAA A %4u    %8.1f %8.1f %8.1f  1.00\n",j,j,X[j],Y[j],Z[j]);
    fprintf(fp,"%f %f %f\n",X[j],Y[j],Z[j]);}
  fprintf(fp1,"END");

  }
