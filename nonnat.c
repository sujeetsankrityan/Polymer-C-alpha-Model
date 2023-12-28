# include <time.h>
# include <stdio.h>
# include <ctype.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <defs.h>
/****************************************************************************/
/* Program nonnat.c                                                         */
/* Ca chain representation, Go-type potential                               */
/* Langevin dynamics                                                        */ 
/************* geometry *****************************************************/
double x[N+NP],y[N+NP],z[N+NP];             /* atom coordinates                      */
double xb[N+NP],yb[N+NP],zb[N+NP];             /* atom coordinates                      */
double b[N+NP];                       /* pseudo bond lenghts (0...N-1)         */
double th[N];                      /* pseudo bond angles (1...N-2)          */
double ph[N];                      /* pseudo torsion angles (1...N-3)       */
double bx[N+NP],by[N+NP],bz[N+NP];          /* pseudo bond vectors                   */
double fbx[N+NP],fby[N+NP],fbz[N+NP],fb[N+NP];          /* pseudo bond vectors                   */
double vbx[N+NP],vby[N+NP],vbz[N+NP],vb[N+NP];          /* pseudo bond vectors                   */
double sx[N+NP],sy[N+NP],sz[N+NP];          /* auxilary vector                       */
double bd_r[N+NP],Cq[N+NP];                 /*Bead radius & Charge                   */  
/************* energies and forces ******************************************/
double Ekin,Epot,Eben,Ebon,Erep;   /* energy terms                          */ 
double Etor,Eel,Ehp;              /* energy terms                          */ 
double fx[N+NP],fy[N+NP],fz[N+NP];          /* conformational force                  */
double fxo[N+NP],fyo[N+NP],fzo[N+NP];       /* conformational force old              */
double frdx[N+NP],frdy[N+NP],frdz[N+NP];    /* random force                          */
double frdxo[N+NP],frdyo[N+NP],frdzo[N+NP]; /* random force old                      */
/************* MD parameters ************************************************/
double vx[N+NP],vy[N+NP],vz[N+NP];          /* velocities                            */
double s[N];
const double tau=4.0;              /* time scale sqrt(m*l*l/eps0)           */
double gam;                        /* friction coefficient                  */
double dt;                         /* time step                             */
/************* interactions *************************************************/
const double eps=1.0;              /* energy unit                           */
double kbon=0.117;                 /* N.B. in unit of eps                   */
double kth=5.0;                   /* N.B. in unit of eps                   */ 
double kph1=0.5;                   /* N.B. in unit of eps                   */ 
double kph2=0.5;                   /* N.B. in unit of eps                   */
double kcon=1.0;                   /* N.B. in unit of eps                   */
double krep=0.01;                   /* N.B. in unit of eps                   */
double sigsa=4.0;                  /* atom radius g(in AA)                   */
double cut=40.0;                    /* cutoff self avoidance                 */
double sigev=28.37;                  /* hydrophobicity                        */
double cut_col=40.0;                  /* cutoff                                */
const double SFC=2529.11892;            /* Columb Screening Factor Constant    */
const double Er=80.0;              /* Permitivity                         */
const double Cs=0.02;               /* Ionic strength                      */
const double Ir=1.4;                 /* Ionic Radius                        */
const double sol_density=1.0;        /* Solvent density                     */

/************* sequence effects *********************************************/
int seq[N+NP];                        /* sequence                              */
int seqhp[N];                      /* sequence hp: >0 hp, 0 polar           */
double kap[N];                     /* hydrophobicity strengths              */
/************* native structure *********************************************/
double bn[N+NP],thn[N],phn[N];        /* bond lengths, bond & torsion angles   */
int npair;                         /* # native contacts                     */
int qpair;                         /* # native contacts                     */
int zpair;                         /* # native contacts                     */
int ip1[MAXP],ip2[MAXP];           /* list of contacts                      */
int iq1[MAXP],iq2[MAXP];           /* list of contacts                      */
int iz1[MAXP],iz2[MAXP];           /* list of contacts                      */
double iqkap[MAXP];
double distp[MAXP];                /* distances                             */
double distp2[MAXP];               /* distances                             */
double distz[MAXP];                /* distances                             */
short cc[N][N];                    /* 1 native contact, 0 otherwise         */
/************* miscellaneous ************************************************/
double c1,c2,c3;                   /* auxilary mdstep parameters            */
double pi,pi2,pid2;                /* pi,2*pi,pi/2                          */
double tconst;                     /* random force magnitude                */
double kcon60,krep12,kbon2,sig2;
double cthmin,sthmin;
int carterr=0,therr=0;
long seed=-13;                     /* random number seed                    */
//const double Cs=0.02;              /* Ionic strength (M)                    */
double dielectric;
double kcoulumb=322.0;             /*coulumb strength                       */
/*--------------------------------------------------------------------------*/
/* external functions */
double rmsd2(int iflag,double x[],double y[],double z[],int n);

/* internal functions */
void movedna(double dx,double dy,double dz);
void movepro(double dx,double dy,double dz);
double cmdna(double *xcm,double *ycm,double *zcm);
double cmpro(double *xcm,double *ycm,double *zcm);
void in2box();
void ch2box();
double gasdev2();
double delta2();
void printhead();
void ramachan(char *fn,double b,double th,double ph);
void runtime(char *fn,long it,double o[],int n);
void write_conf(int iflag,char *fn,double r);
void dumppdb(void);
void movie_make(int imd);
void dumppdb2(char *fn);
void init(void); 
double bond();
double bend();
double torsion();
double coulumb(int iflag);
double cont(int iflag);
double hp(int iflag);
double sac(int iflag,double emax);
double gyr2();
double writhe();
double twist();
int no_cont();
void mdstep();
int cart2dof();
void dof2cart();



int main (int argc,char *argv[])
{
  int i,j,k,iblock=0;
  double o[NOBS],so[NOBS],po[NOBS];
  double nn1,rmsd,rg;
  long imd;
  char fn[100];
  FILE *fp;

  //* printf("Exec: %s\n\n",argv[0]);

  for (i=0;i<NOBS;i++) o[i]=po[i]=so[i]=0;

  init();

  //printhead();
  
    double dx1,dy1,dz1,dx2,dy2,dz2,ds;
    int jst1=0,jst2=50,k1,k2,k3,k2_end,cnt_jxt[N][5];
    for (k3=0;k3<N;k3++) {cnt_jxt[k3][0]=0;cnt_jxt[k3][1]=0;cnt_jxt[k3][2]=0;cnt_jxt[k3][3]=0;cnt_jxt[k3][4]=0;}


        double wrt_o,twst_o,lk_no;
        wrt_o= writhe();
	twst_o = (twist());
	lk_no = wrt_o+twst_o;
    //    printf ("\n%f  %f  \n\n\n",wrt_o,twst_o); exit(-1);
        //printf("%f %f %f\n",wrt_o,twst_o,lk_no);

  for (imd=0;imd<MDCYC;imd++) {
    for (j=0;j<MDCON;j++) {
      mdstep();



    }

    /***************protein Distance********/
    double d_x1,d_y1,d_z1,d_x2,d_y2,d_z2,d1,d2,min1,min2;
    int ind1,ind2;
    min1=10000.0;min2=10000.0;
    for(i=0;i<N;i++) {
      d_x1=x[N]-x[i];d_x2=x[N+1]-x[i];
      d_y1=y[N]-y[i];d_y2=y[N+1]-y[i];
      d_z1=z[N]-z[i];d_z2=z[N+1]-z[i];
      d1=sqrt(d_x1*d_x1+d_y1*d_y1+d_z1*d_z1);
      d2=sqrt(d_x2*d_x2+d_y2*d_y2+d_z2*d_z2);
      if (d1<min1) {min1=d1;ind1=i;}
      if (d2<min2) {min2=d2;ind2=i;}
}

    /***************Juxtaposition Count********/
    double d_jx1,d_jy1,d_jz1,d_jx2,d_jy2,d_jz2,dj1,dj2,minj1,minj2,cut_j;
    int indj1[20],indj2[20],cntj;
    minj1=10000.0;minj2=10000.0;cut_j=40.0;cntj=0;
    for(i=0;i<N;i++) {
      for (j=i;j<N;j++) {
	if((j-i)<4) continue;
	d_jx1=x[j]-x[i];
	d_jy1=y[j]-y[i];
	d_jz1=z[j]-z[i];
	dj1=sqrt(d_jx1*d_jx1+d_jy1*d_jy1+d_jz1*d_jz1);
	//printf("%f %f\n",dj1);
	if (dj1<cut_j) {indj1[cntj]=j;indj2[cntj]=i;cntj+=1;}
      }
    }

    //if ((imd+1)%IRT==0) {
    // for (j=0;j<cntj;j++) {printf("%d %d ",indj1[j],indj2[j]);}
    // printf("\n");}
    // printf("%d %f %d %f\n",ind1,min1,ind2,min2);
        double wrt,twst,lk_no;
	wrt = writhe();
	twst = twist();
	lk_no=wrt+twst;
	//printf("%f %f %f\n",wrt,twst,lk_no);
      //printf("%f %f %f\n",Etor,Erep,Eel);
    o[0]=Ekin; o[1]=Epot; o[2]=Ebon; o[3]=Eben;o[4]=Etor; o[5]=Erep; o[6]=Eel;
    o[9]=ind1;o[10]=min1;o[11]=ind2;o[12]=min2;o[13]=wrt;o[14]=twst;//o[15]=indj1[0];o[16]=indj2[0];o[17]=indj1[1];o[18]=indj2[1];o[19]=indj1[2];o[20]=indj2[2];o[21]=indj1[3];o[22]=indj2[3];o[23]=indj1[4];o[24]=indj2[4];


    if ((imd+1)>NTHERM) {
      for (i=0;i<NOBS;i++) po[i]+=o[i];

      iblock++;
      if (iblock==IBLOCK) {
	for (i=0;i<NOBS;i++) so[i]+=po[i];
	fp=fopen(DATA,"a");
	fwrite(po,sizeof(double),NOBS,fp);
	for (i=0;i<NOBS;i++) po[i]=0;
	fclose(fp);
	iblock=0;
      }

    }

    if ((imd+1)%IRT==0) runtime(RT,imd+1,o,NOBS);
    if ((imd+1)%ICONF==0) write_conf(1,CONF,rg);
    if ((imd+1)%IPDB==0) dumppdb();
    if ((imd+1)%ICONF==0) movie_make(imd);
  }


  printf("\nRun over\n\n");
  dumppdb2("_stop.pdb");
  if (iblock!=0) printf("iblock!=0, %i\n",iblock);
  printf("carterr %i therr %i\n\n",carterr,therr);

  printf("Averages\n");
  for (i=0;i<NOBS;i++) printf("obs %i %.8f\n",i,so[i]/(MDCYC-NTHERM));

  printf("\nWriting histograms\n");
  return 0;
}
/****************************************************************************/
/***** ENERGIES & FORCES ****************************************************/
/****************************************************************************/
double bond() {
  int i,j;
  double fbx,fby,fbz,fb;
  double e=0,db;
  double f_bx,f_by,f_bz;
  f_bx=0;f_by=0;f_bz=0;

  db=b[0]-bn[0];
  e+=db*db;
  //x[0]+=kbon2*(b[0]-bn[0])*bx[0];
  //[0]+=kbon2*(b[0]-bn[0])*by[0];
  //[0]+=kbon2*(b[0]-bn[0])*bz[0];
  fx[0]+=kbon2*((b[0]-bn[0])*bx[0]-(b[N-1]-bn[N-1])*bx[N-1]);
  fy[0]+=kbon2*((b[0]-bn[0])*by[0]-(b[N-1]-bn[N-1])*by[N-1]);
  fz[0]+=kbon2*((b[0]-bn[0])*bz[0]-(b[N-1]-bn[N-1])*bz[N-1]);
  
  //  f_bx+=fx[0];
  //  f_by+=fy[0];
  //  f_bz+=fz[0];


  for (i=1;i<N-1;i++) {
    
    //b[-1]=b[N-1]; bn[-1]=bn[N-1];bx[-1]=bx[N-1];by[-1]=by[N-1];bz[-1]=bz[N-1];
    db=b[i]-bn[i];
    e+=db*db;
   fx[i]-=kbon2*((b[i-1]-bn[i-1])*bx[i-1]-(b[i]-bn[i])*bx[i]);
    fy[i]-=kbon2*((b[i-1]-bn[i-1])*by[i-1]-(b[i]-bn[i])*by[i]);
    fz[i]-=kbon2*((b[i-1]-bn[i-1])*bz[i-1]-(b[i]-bn[i])*bz[i]);
    // printf("***** %d %f %f %f %f %f %f\n",i, fx[i], fy[i], fz[i], b[i], bn[i], db);
    //printf("Sujeet");
    //    f_bx+=fx[i];
    //    f_by+=fy[i];
    //    f_bz+=fz[i];
    
    
  }

  fx[N-1]+=kbon2*(-(b[N-2]-bn[N-2])*bx[N-2]+(b[N-1]-bn[N-1])*bx[N-1]);
  fy[N-1]+=kbon2*(-(b[N-2]-bn[N-2])*by[N-2]+(b[N-1]-bn[N-1])*by[N-1]);
  fz[N-1]+=kbon2*(-(b[N-2]-bn[N-2])*bz[N-2]+(b[N-1]-bn[N-1])*bz[N-1]);
  //  f_bx+=fx[N-1];
  //  f_by+=fy[N-1];
  //  f_bz+=fz[N-1];




  // printf("#bond##%f %f %f\n",f_bx,f_by,f_bz);



  e+=(b[N]-bn[N])*(b[N]-bn[N]);
  fx[N]+=kbon2*(b[N]-bn[N])*bx[N];
  fy[N]+=kbon2*(b[N]-bn[N])*by[N];
  fz[N]+=kbon2*(b[N]-bn[N])*bz[N];
  fx[N+1]-=kbon2*(b[N]-bn[N])*bx[N];
  fy[N+1]-=kbon2*(b[N]-bn[N])*by[N];
  fz[N+1]-=kbon2*(b[N]-bn[N])*bz[N];
  
  return kbon*e;
}
/****************************************************************************/
double bend() {
  int i,j,k;
  double b1x,b1y,b1z,b1;
  double b2x,b2y,b2z,b2;
  double dix,diy,diz;
  double dkx,dky,dkz;
  double cth,sth,dth;
  double cph,sph,dph;
  double aa,bb;
  double e=0,fben;
  double f_bex,f_bey,f_bez;
  f_bex=0;f_bey=0;f_bez=0;

  for (i=0;i<N;i++) {
    if (i<N-2){
      j=i+1;      k=i+2;
    }else if (i==N-2) {
	j=i+1;	k=0;
    } else {j=0; k=1;}
    dth=(th[i]-pi);
    cth=cos(th[i]);
    sth=sin(th[i]);

    e+=dth*dth;
    aa=-2*kth*dth/(sth*b[i]); bb=-2*kth*dth/(sth*b[j]);
    
    dix=aa*(-bx[i]*cth - bx[j]); diy=aa*(-by[i]*cth - by[j]); 
    diz=aa*(-bz[i]*cth - bz[j]);
    dkx=bb*(bx[j]*cth + bx[i]); dky=bb*(by[j]*cth + by[i]); 
    dkz=bb*(bz[j]*cth + bz[i]);
    

    fx[i]+=dix;
    fy[i]+=diy;
    fz[i]+=diz;
    fx[j]+=(-dix-dkx);
    fy[j]+=(-diy-dky);
    fz[j]+=(-diz-dkz);
    fx[k]+=dkx;
    fy[k]+=dky;
    fz[k]+=dkz;

  }
 for (i=0;i<N;i++) {
    f_bex+=fx[i];
    f_bey+=fy[i];
    f_bez+=fz[i];
 }
 //printf("#bend##%f %f %f\n",f_bex,f_bey,f_bez);

  return kth*e;
  }
/****************************************************************************/
/****************************************************************************/
/*
double torsion() {
  int i,j,k;
  double e=0,fben,dph,fact;

  for (i=0;i<(N-1);i++) {
    dph=(ph[i]-pi);
    e+=dph*dph;
    
    fact=-(dph/sin(dph));
    if (sin(dph)<0.01) fact=1.0;

    fx[i]=2*kph1*(fact)*(1/pow((1+cos(th[i])),2))*vbx[i+1]*(-by[i]*fbz[i]+bz[i]*fby[i]);
    fy[i]=2*kph1*(fact)*(1/pow((1+cos(th[i])),2))*vby[i+1]*(-bz[i]*fbx[i]+bx[i]*fbz[i]);
    fz[i]=2*kph1*(fact)*(1/pow((1+cos(th[i])),2))*vbz[i+1]*(-bx[i]*fby[i]+by[i]*fbx[i]);

    fx[i+2]=2*kph1*(fact)*(1/pow((1+cos(th[i])),2))*vbx[i]*(-by[i+1]*fbz[i+1]+bz[i+1]*fby[i+1]);
    fy[i+2]=2*kph1*(fact)*(1/pow((1+cos(th[i])),2))*vby[i]*(-bz[i+1]*fbx[i+1]+bx[i+1]*fbz[i+1]);
    fz[i+2]=2*kph1*(fact)*(1/pow((1+cos(th[i])),2))*vbz[i]*(-bx[i+1]*fby[i+1]+by[i+1]*fbx[i+1]);

    fx[i+1]=2*kph1*(fact)*(1/pow((1+cos(th[i])),2))*((vbx[i+1]*(by[i]*fbz[i]-bz[i]*fby[i]))+vbx[i]*(-by[i+1]*fbz[i+1]+bz[i+1]*fby[i+1]));
    fy[i+1]=2*kph1*(fact)*(1/pow((1+cos(th[i])),2))*((vby[i+1]*(bz[i]*fbx[i]-bx[i]*fbz[i]))+vby[i]*(-bz[i+1]*fbx[i+1]+bx[i+1]*fbz[i+1]));
    fz[i+1]=2*kph1*(fact)*(1/pow((1+cos(th[i])),2))*((vbz[i+1]*(bx[i]*fby[i]-by[i]*fbx[i]))+vbz[i]*(-bx[i+1]*fby[i+1]+by[i+1]*fbx[i+1]));
    */
    
    /*
    fact=-(dph/tan(dph));
    if (tan(dph)<0.0001) fact=1.0;
    fx[i]=2*kph1*(fact)*(1/(1+cos(th[i])))*(bx[i+1]);
    fy[i]=2*kph1*(fact)*(1/(1+cos(th[i])))*(by[i+1]);
    fz[i]=2*kph1*(fact)*(1/(1+cos(th[i])))*(bz[i+1]);
    fx[i+2]=-2*kph1*(fact)*(1/(1+cos(th[i])))*(bx[i]);
    fy[i+2]=-2*kph1*(fact)*(1/(1+cos(th[i])))*(by[i]);
    fz[i+2]=-2*kph1*(fact)*(1/(1+cos(th[i])))*(bz[i]);
    fx[i+1]=2*kph1*(fact)*(1/(1+cos(th[i])))*(bx[i]-bx[i+1]);
    fy[i+1]=2*kph1*(fact)*(1/(1+cos(th[i])))*(by[i]-by[i+1]);
    fz[i+1]=2*kph1*(fact)*(1/(1+cos(th[i])))*(bz[i]-bz[i+1]); 
    
    // printf ("%d %f %f %f %f\n",i,dph,tan(dph),e,th[i]);
    // printf("%d %f %f %f %f %f %f\n",i,fact,fx[i],fx[i+1],fx[i+2],bx[i],bx[i+1]);
  }
  return kph1*e;
  }*/
  
double torsion() {
  int i,j,k;
  double e=0,fben,dph;
  double f_tx,f_ty,f_tz;
  f_tx=0.0;f_ty=0.0;f_tz=0.0;
    dph=(ph[0]);
    e+=dph*dph;
    fx[0]+=2*kph1*(bx[0]*ph[0]-bx[N-1]*ph[N-1]);;
    fy[0]+=2*kph1*(by[0]*ph[0]-by[N-1]*ph[N-1]);;
    fz[0]+=2*kph1*(bz[0]*ph[0]-bz[N-1]*ph[N-1]);;
    f_tx+=fx[0];
    f_ty+=fy[0];
    f_tz+=fz[0];

  for (i=1;i<(N-1);i++) {
    dph=(ph[i]);
    e+=dph*dph;

    // fx[i]+=2*kph1*bx[i]*(ph[i] - ph[i-1]);
    //      fy[i]+=2*kph1*by[i]*(ph[i] - ph[i-1]);
    //      fz[i]+=2*kph1*bz[i]*(ph[i] - ph[i-1]);


    fx[i]+=2*kph1*(bx[i]*ph[i] - bx[i-1]*ph[i-1]);
    fy[i]+=2*kph1*(by[i]*ph[i] - by[i-1]*ph[i-1]);
    fz[i]+=2*kph1*(bz[i]*ph[i] - bz[i-1]*ph[i-1]);

	    f_tx+=fx[i];
	    f_ty+=fy[i];
	    f_tz+=fz[i];
	    
      }
  dph=(ph[N-1]);
  e+=dph*dph;
  fx[N-1]-=2*kph1*(bx[N-2]*ph[N-2]-bx[N-1]*ph[N-1]);
  fy[N-1]-=2*kph1*(by[N-2]*ph[N-2]-by[N-1]*ph[N-1]);
  fz[N-1]-=2*kph1*(bz[N-2]*ph[N-2]-bz[N-1]*ph[N-1]);
  f_tx+=fx[N-1];
  f_ty+=fy[N-1];
  f_tz+=fz[N-1];
 
  //printf("#torsion#%f %f %f\n",f_tx,f_ty,f_tz);

   
  return kph1*e;
  }
/****************************************************************************/
/****** Writhe Calculation ********/
/** PAPER - Transcription driven twin supercoiling of DNA loop **/

double writhe () {
  int i,j;
  float b1x_w,b1y_w,b1z_w,b2x_w,b2y_w,b2z_w,bx_w,by_w,bz_w,b12x_w,b12y_w,b12z_w;
  float bw,mag_bw,b_w;
  float wr=0;
  
  for(i=0;i<N-1;i++) {
    for(j=0;j<N-1;j++) {
      if (i==j) continue;
      b1x_w = x[i+1]-x[i];b1y_w = y[i+1]-y[i];b1z_w = z[i+1]-z[i];
      b2x_w = x[j+1]-x[j];b2y_w = y[j+1]-y[j];b2z_w = z[j+1]-z[j];
      bx_w = x[j]-x[i];by_w = y[j]-y[i];bz_w = z[j]-z[i];
      bw = sqrt(bx_w*bx_w+by_w*by_w+bz_w*bz_w);
      // mag_bw = 1/bw;
      
    b12x_w = b1y_w*b2z_w-b1z_w*b2y_w;
    b12y_w = b1z_w*b2x_w-b1x_w*b2z_w;
    b12z_w = b1x_w*b2y_w-b1y_w*b2x_w;
    
    b_w = b12x_w*bx_w + b12y_w*by_w + b12z_w*bz_w;
    wr+=b_w/pow(bw,3);

    // printf("%d %d %f %f %f %f\n",i,j,wr,mag_bw,bw,b_w);

  }
 }

  return wr/(4*pi);
}

/******************************************************************/
double twist() {
  int i;
  double tw,dh;
  for(i=0;i<N;i++) {
    dh=(ph[i]);
    tw+=dh;
  }
  return tw/(2*pi);
}


 /****************************************************************************/
double  sac(int iflag,double emax) { 
  /* iflag<0 initialize                                    */
  /* iflag=0 calculate full energy                         */
  /* iflag>0 calculate energy but stop if e>emax           */

  extern double cellcalc(int inc,int nl,short list[]);
  static int nx,ny,nz,nc;
  static int h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13;
  int i,j,ix,iy,iz,ic,nec=0;
  int lx[N+NP],ly[(N+NP)],lz[(N+NP)],lc[(N+NP)];
  int in_cell,nl,cn;
  double etest,e=0;
  short a,list[(N+NP)],pnt[(N+NP)];              /* pointers to atoms */
  static short cell[MAXCELL];  /* division into cells */
  static double cutgx,cutgy,cutgz;

  if (iflag<0) {
    for (i=0;i<MAXCELL;i++) cell[i]=-1;
    cellcalc(-1,0,list);
    nx=XBOX/cut; ny=YBOX/cut; nz=ZBOX/cut;
    nc=nx*ny*nz;
    cutgx=XBOX/(double)nx;  cutgy=YBOX/(double)ny;  cutgz=ZBOX/(double)nz;
    if (nc>MAXCELL) {printf("nc=%i TOO BIG, nx=%i ny=%i nz=%i\n",nc,nx,ny,nz); exit(-1);}
    h3=nx*ny; h4=1+nx; h5=h3+1; h6=h3+nx; h7=1-nx; h8=1-h3; h9=nx-h3;
    h10=h4+h3; h11=h4-h3; h12=h3+h7; h13=h3-h7;
    return 0;
  }
  etest=emax/krep;
  if (etest<0 && iflag>0) return krep*e;
  // for(j=0;j<N+NP;j++){printf("%f %f %f\n",x[j],y[j],z[j]);}
  in2box();

  
  for (i=(N+NP)-1;i>=0;i--) {
    ix=xb[i]/cutgx; iy=yb[i]/cutgy; iz=zb[i]/cutgz;
    ic=ix+nx*(iy+ny*iz);
    pnt[i]=cell[ic];
    if (cell[ic]<0) {lx[nec]=ix; ly[nec]=iy; lz[nec]=iz; lc[nec++]=ic;}
    cell[ic]=i;
  }
  i=0;
  while (i<nec && (e<etest || iflag==0)) {
    nl=0;
    list[nl++]=a=cell[(j=lc[i])]; 
    while ((a=pnt[a])>=0) {list[nl++]=a;}
    in_cell=nl;
    cn=j+1; 
    if (lx[i]+1==nx) cn-=nx;
    if ((a=cell[cn])>=0) { 
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+nx;
    if (ly[i]+1==ny) cn-=h3;
    if ((a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+h3;
    if (lz[i]+1==nz) cn-=nc;
    if ( (a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+h4;
    if (lx[i]+1==nx) cn-=nx; if (ly[i]+1==ny) cn-=h3; 
    if ((a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+h5;
    if (lx[i]+1==nx) cn-=nx; if (lz[i]+1==nz) cn-=nc;
    if ((a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+h6;
    if (ly[i]+1==ny) cn-=h3; if (lz[i]+1==nz) cn-=nc; 
    if ((a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+h7;
    if (lx[i]+1==nx) cn-=nx; if (ly[i]==0) cn+=h3; 
    if ((a=cell[cn])>=0) { 
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+h8;
    if (lx[i]+1==nx) cn-=nx; if (lz[i]==0) cn+=nc;
    if ((a=cell[cn])>=0) {  
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+h9;
    if (ly[i]+1==ny) cn-=h3; if (lz[i]==0) cn+=nc; 
    if ((a=cell[cn])>=0) {  
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+h10;
    if (lx[i]+1==nx) cn-=nx; if (ly[i]+1==ny) cn-=h3; if(lz[i]+1==nz) cn-=nc; 
    if ((a=cell[cn])>=0) { 
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+h11;
    if (lx[i]+1==nx) cn-=nx; if (ly[i]+1==ny) cn-=h3; if (lz[i]==0) cn+=nc; 
    if ((a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+h12;
    if (lx[i]+1==nx) cn-=nx; if (ly[i]==0) cn+=h3; if (lz[i]+1==nz) cn-=nc; 
    if ((a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    cn=j+h13;
    if (lx[i]==0) cn+=nx; if (ly[i]+1==ny) cn-=h3; if (lz[i]+1==nz) cn-=nc; 
    if((a=cell[cn])>=0) {
      list[nl++]=a; while ((a=pnt[a])>=0) list[nl++]=a;
    }
    e+=cellcalc(in_cell,nl,list);
    i++;
  }

  for (i=0;i<nec;i++) cell[lc[i]]=-1;

  return krep*e;
  }
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
 double cellcalc(int inc,int nl,short list[])
{
  int i,j,n,m;
  double r2,r6,r12,rx,ry,rz,fr,e=0;
  static double cut2,sigsa2;
  static double asa,bsa,bsa2;
  double len,leni,slen,sleni;
  static double cons;
  static double realT;
  static double ldby,ildby;
  static double cons_pd,cons_dd,diel_inv;



  if (inc<0) {
     cut2=cut*cut;
    sig2=sigev*sigev;
    asa=-(krep/cut2)*6*(-2*pow(sig2/cut2,6.0) - pow(sig2/cut2,3.0));
    bsa=-krep*(13*pow(sig2/cut2,6.0) - 7*pow(sig2/cut2,3.0));
    bsa2=2*bsa;
    //    printf("sac: asa %e bsa %e\n",asa,bsa);
    /*realT=Temp;
    // Re-calculate the dielectric                                                                                                                          
    dielectric = 249.4 - 7.88E-01 * realT + 7.20E-04 * realT * realT;
    dielectric *= (1.000 - (0.2551 * Cs) + (0.05151 * Cs * Cs) - (0.006889 * Cs * Cs * Cs));
    // Calculate Deybe length
    ldby=sqrt(0.000396*dielectric*realT/Cs);
    //ildby=1/ldby;
    cons_pd=sqrt(0.000396*78*realT/Cs);
    cons_dd=sqrt(0.000396*dielectric*realT/Cs);
    diel_inv=1/dielectric;*/
    return 0.0;
  }


  for (n=0;n<inc;n++) {
    i=list[n];
    for (m=n+1;m<nl;m++) {
      j=list[m];
      // if (abs(j-i)<3) continue;
      rx=x[j]-x[i];
      ry=y[j]-y[i];
      rz=z[j]-z[i];
      r2=rx*rx+ry*ry+rz*rz;
     
      if (r2<1024) {
	  //	e-=krep*pow(r2,0.5);
	  //	fr=krep;
	  //printf("\n ***********************************");
	//	if (i<N && j<N) {e+=0; continue;}// r6=sig2/r2;}
	//	if (i>=N || j>=N) {r6=0.66*sig2/r2;}
	r6=sig2/r2;
	  r6=r6*r6*r6; r12=r6*r6;
	  //printf("\n %f %f %i %i",r2,r12-r6, i, j);
	  e+=(r12-r6) + asa*r2 +bsa ;
	  fr=krep*((12*r12-6*r6)/r2 - 2*asa);
	  // if (i>N || j>N) printf("%d %d %f %f\n",i,j,(r12-r6) + asa*r2 +bsa,fr);
	  fx[i]-=fr*rx;
	  fy[i]-=fr*ry;
	  fz[i]-=fr*rz;
	  fx[j]+=fr*rx;
	  fy[j]+=fr*ry;
	  fz[j]+=fr*rz;
      }
    }
  }
  return e;
}

 /***************************************************************************/
/****************************************************************************/
double coulumb(int iflag) {
  int i,j;
  double rx,ry,rz,fcoul,cons=0;
  double r2,e=0,len,leni,slen,sleni,fbx,fby,fbz;
  static double SF,realT,Bk;
  static double cutg2;

  if (iflag < 0) {
    cutg2=cut_col*cut_col;
    /* We consider T_cg=350K and for each 0.1 change in T_cg, real temp changes by 10 K */
    realT=Temp;//-((1-(1/beta[ind]))*100);
    SF=sqrt((SFC*Cs*sol_density)/(Er*realT));
    Bk=(exp(SF*Ir))/(1+SF*Ir);
    return 0.0;
  }

  cons=(kcoulumb*Bk/Er);
  // printf("%f %f %f",SF,Bk,cons);exit(-1);
  for (i=0;i<N;i++) {
    for (j=i+1;j<N;j++) {
      //    j=N;
      rx=x[j]-x[i]; //while (rx>(0.5*XBOX)) rx-=XBOX; while (rx<(-0.5*XBOX)) rx+=XBOX;
      ry=y[j]-y[i]; //while (ry>(0.5*YBOX)) ry-=YBOX; while (ry<(-0.5*YBOX)) ry+=YBOX;
      rz=z[j]-z[i]; //while (rz>(0.5*ZBOX)) rz-=ZBOX; while (rz<(-0.5*ZBOX)) rz+=ZBOX;
      
      if (((r2=rx*rx+ry*ry+rz*rz) > cutg2)) continue;
      //      if ((i<N && j<N) || (i>=N && j>= N)) continue;
      len=sqrt(r2); leni=1/len; slen=SF*len; sleni=SF*leni;
      e+=leni*Cq[i]*Cq[j]*exp(-slen);
      fcoul=(leni*cons*Cq[i]*Cq[j]*exp(-slen))*(sleni + (1/r2));
      
      fbx=fcoul*rx;
      fby=fcoul*ry;
      fbz=fcoul*rz;
      
      fx[i]-=fbx;
      fy[i]-=fby;
      fz[i]-=fbz;
      fx[j]+=fbx;
      fy[j]+=fby;
      fz[j]+=fbz;
     
    }
  }

  cons=(kcoulumb*Bk/Er);
  for (i=0;i<N;i++) {
    for (j=N;j<N+NP;j++) {
      //    j=N;
      rx=x[j]-x[i]; //while (rx>(0.5*XBOX)) rx-=XBOX; while (rx<(-0.5*XBOX)) rx+=XBOX;
      ry=y[j]-y[i]; //while (ry>(0.5*YBOX)) ry-=YBOX; while (ry<(-0.5*YBOX)) ry+=YBOX;
      rz=z[j]-z[i]; //while (rz>(0.5*ZBOX)) rz-=ZBOX; while (rz<(-0.5*ZBOX)) rz+=ZBOX;
      
      if (((r2=rx*rx+ry*ry+rz*rz) > cutg2)) continue;
      //      if ((i<N && j<N) || (i>=N && j>= N)) continue;
      len=sqrt(r2); leni=1/len; slen=SF*len; sleni=SF*leni;
      e+=leni*Cq[i]*Cq[j]*exp(-slen);
      fcoul=(leni*cons*Cq[i]*Cq[j]*exp(-slen))*(sleni + (1/r2));
      
      fbx=fcoul*rx;
      fby=fcoul*ry;
      fbz=fcoul*rz;
      
      fx[i]-=fbx;
      fy[i]-=fby;
      fz[i]-=fbz;
      fx[j]+=fbx;
      fy[j]+=fby;
      fz[j]+=fbz;
     
    }
  }



  return cons*e;
}
/****************************************************************************/

/****************************************************************************/

/*    Electrostatic       */

/*double coulumb(int iflag) {
  int i,j;
  double rx,ry,rz,fcoul,cons;
  double r2,e=0,len,leni,slen,sleni,fbx,fby,fbz;
  static double cutg2;
  static double realT;
  static double ldby,ildby;
  static double cons_pd,cons_dd,diel_inv;




  if (iflag < 0) {
    realT=300;   
    // Re-calculate the dielectric                                                                                                                          
    dielectric = 249.4 - 7.88E-01 * realT + 7.20E-04 * realT * realT;
    dielectric *= (1.000 - (0.2551 * Cs) + (0.05151 * Cs * Cs) - (0.006889 * Cs * Cs * Cs));

    // Calculate Deybe length
    ldby=sqrt(0.000396*dielectric*realT/Cs);
    ildby=1/ldby;
    cons_pd=sqrt(0.000396*78*realT/Cs);
    cons_dd=sqrt(0.000396*dielectric*realT/Cs);
    diel_inv=1/dielectric;
    return 0.0;
  }
  // printf ("%f %f %f %f %f %f\n",dielectric,ldby,cons_pd,cons_dd,diel_inv,kcoulumb);exit(-1);
  for (i=0; i<(N+NP-1); i++) {
    for (j=i+1; j<N+NP; j++) {
      if(i>=N && j>=N) continue;
      cons=(kcoulumb*diel_inv); 
      rx=x[j]-x[i];//  while (rx>(0.5*XBOX)) rx-=XBOX; while (rx<(0.5*XBOX)) rx+=XBOX;
      ry=y[j]-y[i];//  while (ry>(0.5*YBOX)) ry-=YBOX; while (ry<(0.5*YBOX)) ry+=YBOX;
      rz=z[j]-z[i];//  while (rz>(0.5*ZBOX)) rz-=ZBOX; while (rz<(0.5*ZBOX)) rz+=ZBOX;
  
      if ((r2=rx*rx+ry*ry+rz*rz) >= 2500) continue;
      
      if ((i>=N || j>=N)) {
	ldby=cons_pd;
	cons=cons*1.67;
      // printf("*****%f %f %f\n ",sqrt(r2),Cq[i],Cq[j]);
      } else if ((i<N && j<N)) {
	e+=0; continue;
	ldby=cons_dd; 
	cons=cons;
      }

      
     //printf ("%f  \n",sqrt(r2)); 
      ildby=1/ldby;
      len=sqrt(r2); leni=1/len; slen=ildby*len; sleni=ildby*leni;
      e+=cons*leni*Cq[i]*Cq[j]*exp(-slen);
      fcoul=(leni*cons*Cq[i]*Cq[j]*exp(-slen))*(sleni + (1/r2));
      if ((i>=N || j>=N)) { printf("%d %d %f %f %f %f %f\n",i,j,len,cons*leni*Cq[i]*Cq[j]*exp(-slen),fcoul,Cq[i],Cq[j]);}
      //  printf("%f %f\n",cons*leni*Cq[i]*Cq[j]*exp(-slen),e);
       
      fbx=fcoul*rx;
      fby=fcoul*ry;
      fbz=fcoul*rz;
      fx[i]+=fbx;
      fy[i]+=fby;
      fz[i]+=fbz;
      fx[j]-=fbx;
      fy[j]-=fby;
      fz[j]-=fbz;

      //      printf("%f %f %f\n",fx[i],fy[i],fz[i]);

      
      // rx=x[Iq[j]]-x[Iq[i]]; while (rx>xboxhf) rx-=XBOX; while (rx<-xboxhf) rx+=XBOX;
      //ry=y[Iq[j]]-y[Iq[i]]; while (ry>yboxhf) ry-=YBOX; while (ry<-yboxhf) ry+=YBOX;
      //rz=z[Iq[j]]-z[Iq[i]]; while (rz>zboxhf) rz-=ZBOX; while (rz<-zboxhf) rz+=ZBOX;
      
    
    }
  }
  
  return e;
}
*/


/****************************************************************************/
double cont(int iflag) {
  int i,j,m;
  double r2,r4,r6,rx,ry,rz,fr,e=0;
  double static asa,bsa,bsa2;

  if (iflag<0) {
    bsa=1905*T*(2*pow(0.893,14)-pow(0.893,8));
    asa=-2400*T*(13*pow(0.893,12)-7*pow(0.893,6));
    bsa2=2*bsa;
    //    printf("cont: asa %f bsa %f\n",asa,bsa);
    }

  for (i=0;i<N+NP-1;i++) {
    for (j=i+1;j<N+NP;j++) {
      //      if (j-i <3) continue;
      rx=x[j]-x[i]; 
      ry=y[j]-y[i]; 
      rz=z[j]-z[i];
      if ((r2=rx*rx+ry*ry+rz*rz)>=1.2599*sig2) continue; 
      r6=sig2/r2; r4=r6*r6; r6=r4*r6;
      e+=r6*(r6-1);//+asa+bsa*r2/sig2;
      fr=400*T*r6*(12*r6-6)/r2;//-bsa2/sig2;
      fx[i]-=fr*rx;
      fy[i]-=fr*ry;
      fz[i]-=fr*rz;
      fx[j]+=fr*rx;
      fy[j]+=fr*ry;
      fz[j]+=fr*rz;
    }
  }
  return 400*T*e;
  }

/****************************************************************************/
/*double coulumb(int iflag) {
  int i,j;
  double rx,ry,rz,fcoul;
  double r2,e=0,len,leni,slen,sleni,fbx,fby,fbz;
  static double SF,realT,Bk;
  static double cut2=4.0,cons;

  if (iflag < 0) {
    //    cut2=cut*cut; 
    //    realT=300-((1-(1/T))*100);
    //    SF=(sqrt((2529.11829*0.2)/(80.0*realT)));
    //    Bk=(exp(SF*1.4))/(1+SF*1.4);
    return 0.0;
  }

  //  cons=(332*Bk/80.0);
  cons=85.7;
  for (i=0;i<N-1;i++) {
    for (j=i+1;j<N;j++) {

      rx=x[i]-x[j]; 
      ry=y[i]-y[j]; 
      rz=z[i]-z[j]; 

      if (((r2=rx*rx+ry*ry+rz*rz) < cut2)) continue;
      len=sqrt(r2); leni=1/len; slen=1.4711*len; sleni=1.4711*leni;   
      e+=leni*exp(-slen);
      fcoul=(leni*cons*exp(-slen))*(sleni + (1/r2));

      fbx=fcoul*rx;
      fby=fcoul*ry;
      fbz=fcoul*rz;
      fx[i]-=fbx;
      fy[i]-=fby;
      fz[i]-=fbz;
      fx[j]+=fbx;
      fy[j]+=fby;
      fz[j]+=fbz;
    }
  }

  return cons*e;
  }*/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/***** OBSERVABLES **********************************************************/
/****************************************************************************/
/*int no_cont() {
  int i,j,m,nc=0;

  for (m=0;m<npair;m++) {
    i=ip1[m]; j=ip2[m];
    if ((x[i]-x[j])*(x[i]-x[j])+
	(y[i]-y[j])*(y[i]-y[j])+
	(z[i]-z[j])*(z[i]-z[j])>1.44*distp2[m]) continue;
    nc++;
  }
  return nc;
  }*/
/****************************************************************************/
/****************************************************************************/
/*double delta2() {
  return rmsd2(1,x,y,z,N);
  }*/


/****************************************************************************/
/***** GEOMETRY *************************************************************/
/****************************************************************************/
void dof2cart() {
  int i,j,k;
  double vx,vy,vz;
  double sth1,sth2,cth1,cth2,cph,sph;

  x[0]=y[0]=z[0]=0;
  bx[0]=by[0]=0; bz[0]=1.0;
  x[1]=y[1]=0; z[1]=b[0];
  
  bx[1]=0;
  by[1]=sin(th[1]);
  bz[1]=-cos(th[1]);

  x[2]=0;
  y[2]=b[1]*by[1];
  z[2]=z[1]+b[1]*bz[1];

  for (k=2;k<N-1;k++) {
    j=k-1;
    i=k-2;
    
    sth1=sin(th[j]);
    sth2=sin(th[k]);
    cth1=cos(th[j]);
    cth2=cos(th[k]);
    sph=sin(ph[j]);
    cph=cos(ph[j]);

    sx[j]=(by[i]*bz[j]-bz[i]*by[j])/sth1;
    sy[j]=(bz[i]*bx[j]-bx[i]*bz[j])/sth1;
    sz[j]=(bx[i]*by[j]-by[i]*bx[j])/sth1;

    vx=(bx[i]+bx[j]*cth1)/sth1;
    vy=(by[i]+by[j]*cth1)/sth1;
    vz=(bz[i]+bz[j]*cth1)/sth1;

    bx[k]=-cth2*bx[j]+sth2*(sph*sx[j]-cph*vx);
    by[k]=-cth2*by[j]+sth2*(sph*sy[j]-cph*vy);
    bz[k]=-cth2*bz[j]+sth2*(sph*sz[j]-cph*vz);
    
    x[k+1]=x[k]+b[k]*bx[k];
    y[k+1]=y[k]+b[k]*by[k];
    z[k+1]=z[k]+b[k]*bz[k];
  } 
}
/****************************************************************************/
int cart2dof() {
  int i,j,ok=1;
  double b1x,b1y,b1z,b1;
  double b2x,b2y,b2z;
  double ux,uy,uz,u;
  double tmp1;


  /* Charge Distribution */
  for (i=0;i<N+NP;i++) {
    if (i<N) Cq[i]=-10.0;
    else if (i==N) Cq[i]=6.0;
    else if (i==N+1) Cq[i]=5.0;}

  for (i=0;i<N+NP;i++){
    //    if (i<N) bd_r[i]=15.92;
    if (i<N) bd_r[i]=15.9;
    else if (i==N) bd_r[i]=10;
    else if (i==N+1) bd_r[i]=10;}



  /* bond vectors for DNA*/
  double f1, v1;
  for (i=0;i<N;i++) {
    if (i<(N-1)) j=i+1;
    if (i==(N-1)) j=0;
      
    b1x=x[j]-x[i];
    b1y=y[j]-y[i];
    b1z=z[j]-z[i];
    b1=sqrt(b1x*b1x+b1y*b1y+b1z*b1z);
    
    bx[i]=b1x/b1;
    by[i]=b1y/b1;
    bz[i]=b1z/b1;
    b[i]=b1;

    fbx[i]=1; fby[i]=-((bx[i]+bz[i])/by[i]); fbz[i]=1; 
    f1=sqrt(fbx[i]*fbx[i]+fby[i]*fby[i]+fbz[i]*fbz[i]);
    fbx[i]=fbx[i]/f1;
    fby[i]=fby[i]/f1;
    fbz[i]=fbz[i]/f1;
    fb[i]=f1;

    vbx[i]=by[i]*fbz[i]-bz[i]*fby[i];
    vby[i]=bz[i]*fbx[i]-bx[i]*fbz[i];
    vbz[i]=bx[i]*fby[i]-by[i]*fbx[i];
    v1=sqrt(vbx[i]*vbx[i]+vby[i]*vby[i]+vbz[i]*vbz[i]);
    vbx[i]=vbx[i]/v1; vby[i]=vby[i]/v1; vbz[i]=vbz[i]/v1;
  }


  /* bond vectors for protein*/
  b1x=x[N+1]-x[N]; while (b1x>(0.5*XBOX)) b1x-=XBOX; while (b1x<-(0.5*XBOX)) b1x+=XBOX;
  b1y=y[N+1]-y[N]; while (b1y>(0.5*YBOX)) b1y-=YBOX; while (b1y<-(0.5*YBOX)) b1y+=YBOX;
  b1z=z[N+1]-z[N]; while (b1z>(0.5*ZBOX)) b1z-=ZBOX; while (b1z<-(0.5*ZBOX)) b1z+=ZBOX;
  b1=sqrt(b1x*b1x+b1y*b1y+b1z*b1z);
  
  bx[N]=b1x/b1;
  by[N]=b1y/b1;
  bz[N]=b1z/b1;
  b[N]=b1;
  


  /* Bend/Torsion angles */
  int k;
  double tmp2,tmp3, tmp4, tmp5;
  double sum=0,ag;
  for (i=0;i<N;i++) {
    if(i==(N-1)) {
      j=0;
    } else j=i+1; 
    b1x=bx[i];
    b1y=by[i];
    b1z=bz[i];
    
    b2x=bx[j];
    b2y=by[j];
    b2z=bz[j];
  
    /* calculation of bend angles */
    tmp1=b1x*b2x+b1y*b2y+b1z*b2z;
    th[i]=pi-acos(tmp1);
    // printf("%d %f %f %f\n",i,tmp1,acos(tmp1),th[i]);
    /* calculation of torsion angles */    
    tmp2=(fbx[j]*fbx[i]+fby[j]*fby[i]+fbz[j]*fbz[i]);
    tmp3=(bx[i]*(fby[i]*fbz[j] - fby[j]*fbz[i])+by[i]*(fbz[i]*fbx[j] - fbz[j]*fbx[i])+bz[i]*(fbx[i]*fby[j] - fbx[j]*fby[i]));
    //    tmp3=(vbx[i+1]*vbx[i]+vby[i+1]*vby[i]+vbz[i+1]*vbz[i]);

    //    tmp4=(fbx[i+1]*vbx[i]+fby[i+1]*vby[i]+fbz[i+1]*vbz[i]);
    //    tmp5=(vbx[i+1]*fbx[i]+vby[i+1]*fby[i]+vbz[i+1]*fbz[i]);

    //ph[i]=asin(tmp3);
    if (tmp2 > 0 ){    ph[i]=asin(tmp3);}
    else if (tmp2 < 0 && tmp3>0 ){    ph[i]=pi - asin(tmp3);}
    else if (tmp2 < 0 && tmp3<0 ){    ph[i]=-pi - asin(tmp3);}
    if (ph[i] < 0){s[i]= -1;} else {s[i]=1;}

    //    ph[i]=acos((tmp2+tmp3)/(1+tmp1));
    //printf("%d %f %f %f\n",i, acos((tmp2+tmp3)/(1+tmp1)), asin((tmp4-tmp5)/(1+tmp1)), ph[i]);
    // printf("%d %f %f %f %f\n",i,ph[i],tmp2,tmp3,tmp2+tmp3);
    sum+=ph[i];
  }
  //printf(" %f \n", (sum/6.28));
   //exit(-1);

  return ok;
}
/****************************************************************************/
/***** INPUT/OUTPUT *********************************************************/
/****************************************************************************/
void printhead() {
  printf("\n***************************************************\n");
  printf("nonnat.c (2019-06-11)\n");
  printf("N %i MDCYC %i MDCON %i NTHERM %i\n",N,MDCYC,MDCON,NTHERM);
  printf("IRT %i ICONF %i IBLOCK %i ISTART %i\n",IRT,ICONF,IBLOCK,ISTART);
  printf("T %f\n",T);
  printf("Interaction parameters\n");
  printf("kbon %f kth %f krep %f \n",kbon,kth,krep);
  printf("kph1 %f kph2 %f kcon %f \n",kph1,kph2,kcon);
  printf("eps %f sigsa %f \n",eps,sigsa);
  printf("MD parameters\n");
  printf("tau %f dt %f gam %f\n",tau,dt,gam);
  printf("c1 %f c2 %f c3 %f\n",c1,c2,c3);
  printf("----\n\n");
  fflush(0);
}
/****************************************************************************/
/*void ramachan(char *fn,double b,double th,double ph) {
  FILE *fp;

  fp=fopen(fn,"a");
  fprintf(fp,"%.4f %.4f %.4f\n",b,th*180/pi,ph*180/pi);
  fclose(fp);
  }*/
/****************************************************************************/
void runtime(char *fn,long it,double o[],int n) {
  int i;
  FILE *fp;
  
  fp=fopen(fn,"a");
  fprintf(fp,"%li ",it);
  for (i=0;i<n;i++)
    fprintf(fp,"%.4f ",o[i]);
  fprintf(fp,"\n");
  fclose(fp);
} 
/****************************************************************************/
void write_conf(int iflag,char *fn,double r) {
  FILE *fp;
  
  if (iflag<0) {
    fp=fopen(fn,"w");
    fclose(fp);
    return;
  }

  fp=fopen(fn,"a");
  fwrite(b,sizeof(double),N-1,fp);
  fwrite(&th[1],sizeof(double),N-2,fp);
  fwrite(&ph[1],sizeof(double),N-3,fp);
  fwrite(&r,sizeof(double),1,fp);
  fclose(fp);
} 
/****************************************************************************/
void dumppdb(void) {
  int i;
  FILE *fp;
  double xcm=0,ycm=0,zcm=0;

  fp=fopen(PDB,"w");
  for (i=0;i<N+NP;i++) {
    fprintf(fp,"ATOM   %4u  CA  ALA  %4u    %8.3f%8.3f%8.3f  1.00\n",
	    i+1,i+1,x[i],y[i],z[i]); 
  }
  fclose(fp);
}
/****************************************************************************/
/****************************************************************************/
void movie_make(int imd) {
  int i;
  FILE *fp;
  double xcm=0,ycm=0,zcm=0;
  fp=fopen(MOVIE,"a");
 
  ch2box();
  for (i=0;i<N;i++) {
    fprintf(fp,"ATOM   %4u  CA  ALA  %4u    %8.1f%8.1f%8.1f  1.00\n",
	    i+1,i+1,x[i]-xcm,y[i]-ycm,z[i]-zcm); 
  }
  fprintf(fp,"TER\n");
  for (i=N;i<N+NP;i++) {
    fprintf(fp,"ATOM   %4u  CA  LYS  %4u    %8.1f%8.1f%8.1f  1.00\n",
	    i+1,i+1,x[i],y[i],z[i]); 
  }

  fprintf(fp,"ENDMDL\n");
  fclose(fp);
}
/****************************************************************************/

void dumppdb2(char *fn) {
  int i;
  FILE *fp;
  double xcm=0,ycm=0,zcm=0;

  fp=fopen(fn,"w");
  for (i=0;i<N+NP;i++) {
    fprintf(fp,"ATOM   %4u  CA  ALA  %4u    %8.1f%8.1f%8.1f  1.00\n",
	    i+1,i+1,x[i],y[i],z[i]); 
  }
  fclose(fp);
}
/****************************************************************************/
/***** MOLECULAR DYNAMICS ***************************************************/
/****************************************************************************/
void mdstep() {
  int i;

  for (i=0;i<N+NP;i++) {
    x[i]=x[i]+dt*vx[i]+c1*(fx[i]-gam*vx[i]+frdx[i]);  
    y[i]=y[i]+dt*vy[i]+c1*(fy[i]-gam*vy[i]+frdy[i]);
    z[i]=z[i]+dt*vz[i]+c1*(fz[i]-gam*vz[i]+frdz[i]);
    // printf("x[%i]=%f \n",i,x[i]);
  }

  ch2box();  

  if (1!=cart2dof()) carterr++;

  for (i=0;i<N+NP;i++) {
    frdxo[i]=frdx[i];
    frdyo[i]=frdy[i];
    frdzo[i]=frdz[i];
    fxo[i]=fx[i];
    fyo[i]=fy[i];
    fzo[i]=fz[i];
    frdx[i]=gasdev2()*tconst;
    frdy[i]=gasdev2()*tconst;
    frdz[i]=gasdev2()*tconst;
  }

  for (i=0;i<N+NP;i++) {fx[i]=fy[i]=fz[i]=0;}

  Epot=(Ebon=bond())+(Eben=bend())+(Etor=torsion())+(Eel=coulumb(0))+(Erep=sac(0,0));//
  // Epot=(Ebon=bond())+(Eben=bend())+(Eel=coulumb(0))+(Erep=sac(0,0));


  Ekin=0;
  for (i=0;i<N+NP;i++) {
    vx[i]=c2*vx[i]+c3*(fxo[i]+fx[i]+frdxo[i]+frdx[i]);
    vy[i]=c2*vy[i]+c3*(fyo[i]+fy[i]+frdyo[i]+frdy[i]);
    vz[i]=c2*vz[i]+c3*(fzo[i]+fz[i]+frdzo[i]+frdz[i]);
    Ekin+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
  }
  Ekin*=0.5;
  ch2box();

  //  printf("%.15e\n",Ehp-hp2(0));
}
/****************************************************************************/
/***** INITIALIZATION *******************************************************/
/****************************************************************************/
void init(void) {
  int i,j,m;
  double xn[N+NP],yn[N+NP],zn[N+NP];
  double vxsum,vysum,vzsum;
  double r2,r2min,r2max,c0;
  char c;
  FILE *fp1;

  /* read sequence */

  fp1=fopen(INPUT,"r"); 
  for (i=0;i<N+NP;i++) {                       
    while ((c=getc(fp1))==' ' || c==',' || c=='\n') {}
    seq[i]=c;
    //* putchar(c);
  }
  //* printf("\n");
  //* for (i=0;i<N+NP;i++) printf("%i",i%10); 
  //* printf("\n");

  pi=acos(-1.);
  pi2=2.*pi;
  pid2=pi/2.;
  cthmin=cos(pi/90.);
  sthmin=sin(pi/90.);

  /* energy parameters */
  kbon*=eps;
  kth*=eps;
  kph1*=eps;
  kph2*=eps;
  kcon*=eps;
  kcon60=kcon*60;
  krep*=eps;
  krep12=krep*12;
  kbon2=kbon*2;
  sig2=sigev*sigev;
  
  /* MD parameters */
  dt=0.001*tau;
  gam=5/tau;
  tconst=sqrt(2*gam*T/dt);
  c0=gam*dt/2.;
  c1=dt*dt/2.;
  c2=(1-c0)*(1-c0+c0*c0);
  c3=dt*(1-c0+c0*c0)/2.;

  //* printf("seed %li\n",seed);

  /* native structure */
  fp1=fopen(NATIVE,"r");
  for (i=0;i<N+NP;i++) fscanf(fp1,"%lf %lf %lf",&x[i],&y[i],&z[i]);
  fclose(fp1);

  //  for (i=0;i<N+NP;i++) for (j=0;j<N+NP;j++) cc[i][j]=0;
  /* native dof values */
  for (i=0;i<N+NP;i++) {
    xn[i]=x[i];
    yn[i]=y[i];
    zn[i]=z[i];
  }
  double dx,dy,dz;
  cmdna(&dx,&dy,&dz); // printf("%f %f %f\n",dx,dy,dz);
  movedna(-dx+0.5*XBOX,-dy+0.5*YBOX,-dz+0.5*ZBOX);
  //cmdna(&dx,&dy,&dz);  //printf("%f %f %f\n",dx,dy,dz);//exit(-1);
  cmpro(&dx,&dy,&dz); movepro(-dx+0.5*XBOX,-dy+0.5*YBOX+100,-dz+0.5*ZBOX);


  // cont(-1);
   coulumb(-1);
  sac(-1,0);
  write_conf(-1,CONF,0);
  //  rmsd2(-1,xn,yn,zn,N);


  if (1!=cart2dof()) printf("Error native configuration");
  for (i=0;i<=N; i++) bn[i]=b[i];
  //  for (i=0;i<N;i++) bn[i]=31.84;
  //    for (i=0;i<N;i++) bn[i]=4.0;
  //  for (i=1;i<N-1;i++) thn[i]=th[i];
  //  for (i=1;i<N-2;i++) phn[i]=ph[i];

  //  for (i=0;i<N-1;i++) printf("\n b[%i]=%f",i,b[i]);
  //  exit(-1);
  /* initialize velocities */
  vxsum=vysum=vzsum=0;
  for (i=0;i<N+NP;i++) {
    vx[i]=sqrt(T)*gasdev2();
    vy[i]=sqrt(T)*gasdev2();
    vz[i]=sqrt(T)*gasdev2();
    vxsum+=vx[i];
    vysum+=vy[i];
    vzsum+=vz[i];
  }

  Ekin=0;
  for (i=0;i<N+NP;i++) {
    vx[i]-=vxsum/N;
    vy[i]-=vysum/N;
    vz[i]-=vzsum/N;
    Ekin+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
  }
  Ekin*=0.5;


  if (1!=cart2dof()) printf("Error initial configuration");

  for (i=0;i<N+NP;i++) fx[i]=fy[i]=fz[i]=0;
  //Epot=(Etor=torsion())+(Ebon=bond());//+(Erep=sac(0,0));
  Epot=(Ebon=bond())+(Eben=bend())+(Etor=torsion())+(Erep=sac(0,0))+(Eel=coulumb(0));
  //Epot=(Ebon=bond())+(Eben=bend())+(Erep=sac(0,0));
  

  for (i=0;i<N+NP;i++) {
    frdx[i]=gasdev2()*tconst;
    frdy[i]=gasdev2()*tconst;
    frdz[i]=gasdev2()*tconst;
    dumppdb2("_start.pdb");
    //  exit(-1);
  }

}
/****************************************************************************/
double gasdev2() {
  double ran3n(long *seed);
  static int iset=0;
  static double gcos;
  double tmp1,tmp2;

  if (iset==0) { 
    tmp1=sqrt(-2*log(ran3n(&seed)));
    tmp2=pi2*ran3n(&seed);
    gcos=tmp1*cos(tmp2);
    iset=1;
    return tmp1*sin(tmp2);
  }else{
    iset=0;
    return gcos;
  }
}
/****************************************************************************/
void in2box() {
  int i,j,k=0;

  for (i=0; i<N+NP; i++) {
      xb[k]=x[i]; yb[k]=y[i]; zb[k]=z[i];
      while (xb[k]>=XBOX) xb[k]-=XBOX; while (xb[k]<0) xb[k]+=XBOX;
      while (yb[k]>=YBOX) yb[k]-=YBOX; while (yb[k]<0) yb[k]+=YBOX;
      while (zb[k]>=ZBOX) zb[k]-=ZBOX; while (zb[k]<0) zb[k]+=ZBOX;
      k++;
    }
  }
/****************************************************************************/
void ch2box()
{
  int i,j,k=0,nx,ny,nz;

  for (j=0;j<N;j++) {
    nx=ny=nz=0;
    k=(int)N/2;

    while (x[k]+nx*XBOX>XBOX) nx--; while (x[k]+nx*XBOX<0) nx++;
    while (y[k]+ny*YBOX>YBOX) ny--; while (y[k]+ny*YBOX<0) ny++;
    while (z[k]+nz*ZBOX>ZBOX) nz--; while (z[k]+nz*ZBOX<0) nz++;
    for (i=0;i<N;i++) {
      x[i]+=XBOX*nx;
      y[i]+=YBOX*ny;
      z[i]+=ZBOX*nz;
    }
  }

  /*  for (i=0;i<(N+NP);i++) {
    xb[i]=x[i]; 
    yb[i]=y[i]; 
    zb[i]=z[i];
    }*/

}
/****************************************************************************/
/****************************************************************************/
double cmdna(double *xcm,double *ycm,double *zcm) 
{
  int i,n=0;
  (*xcm)=(*ycm)=(*zcm)=0;
  for (i=0;i<N;i++) {
    (*xcm)+=x[i];
    (*ycm)+=y[i];
    (*zcm)+=z[i];
    n++;
  }
  (*xcm)/=n; (*ycm)/=n; (*zcm)/=n;

  return 0;
}
/****************************************************************************/
/****************************************************************************/
void movedna(double dx,double dy,double dz) 
{
  int i;
  for (i=0;i<N;i++) {
    x[i]+=dx;
    y[i]+=dy;
    z[i]+=dz;
  }
}

/****************************************************************************/
/****************************************************************************/
double cmpro(double *xcm,double *ycm,double *zcm) 
{
  int i,n=0;
  (*xcm)=(*ycm)=(*zcm)=0;
  for (i=N;i<N+NP;i++) {
    (*xcm)+=x[i];
    (*ycm)+=y[i];
    (*zcm)+=z[i];
    n++;
  }
  (*xcm)/=n; (*ycm)/=n; (*zcm)/=n;

  return 0;
}
/****************************************************************************/
/****************************************************************************/
void movepro(double dx,double dy,double dz) 
{
  int i;
  for (i=N;i<N+NP;i++) {
    x[i]+=dx;
    y[i]+=dy;
    z[i]+=dz;
  }
}

/****************************************************************************/
