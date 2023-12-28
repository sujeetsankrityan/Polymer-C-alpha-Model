/************* settings *****************************************************/
# define N 100                    /* # amino acids                         */
# define T 0.6                     /* temperature                           */
# define NP 2                      /* 2*Number of Protein residues          */
#define Temp 300                   /* Temperature                           */
/************* MD parameters ************************************************/
# define MDCYC 10000000             /* # cycles                              */
# define MDCON 1000                 /* # md steps per cycle                  */
# define NTHERM 10000              /* # discarded cycles                    */
# define IBLOCK 10000              /* size of data blocks                   */
# define IRT 100                  /* runtime                               */
# define IRAMA 100000              /* ramachandran                          */
# define IPDB 10000                /* pdb                                   */
# define ICONF 100               /* configurations                        */
# define ISTART 0                  /* 0 native, 1 read                      */
# define ISEED 1
/************* measurements *************************************************/
# define NBIN 200                  /* # bins                                */
# define NOBS 25                   /* # observables                         */
# define MAXCELL 2000000             /* max # cells                           */
# define MAXP 500                  /* max # contact pairs                   */
# define MAXNC 1500                 /* max # contacts                        */
# define XBOX 11000     /* # Box length along X-axis */
# define YBOX 11000     /* # Box length along X-axis */
# define ZBOX 1000     /* # Box length along X-axis */
/************* files input **************************************************/
# define NATIVE "native1"
# define CONTMAP "cmap4.5"
# define START "start"
# define INPUT "input_WT" 
# define MOVIE "movie.pdb" 
/************* files output *************************************************/
# define RT "_rt" 
# define RAMA "_rama" 
# define PDB "_pdb"
# define DATA "_data"
# define CONF "_conf"
/************* functions ****************************************************/
# define max(A,B) ((A)>(B)?(A):(B)) 
# define min(A,B) ((A)<(B)?(A):(B))
/****************************************************************************/
