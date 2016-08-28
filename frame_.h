/**************************************************************************************************************/
/*** EXTERNAL LIBRARIES ***/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/**************************************************************************************************************/
/*** FILE FORMATS AND ANALYSIS TYPES ***/

#define FORMAT		DAT
#define ANALYSIS	CORR

#define DAT 0		// Single Variable
#define XYZ 1		// {x, y, z}  Coordinates
#define DIM 2		// G87 Format Dimer Coordinates
#define ATM 3		// G87 Format Monomer Positions
#define G96 4		// G96 Format Monomer Positions and Velocities

#define CORR 0		// Variable: Correlation Function
#define RCOR 1		// Rotation: Correlation Functions
#define RMSD 2		// Rotation: Mean Squared Displacement
#define TMSD 3		// Translation: Mean Squared Displacement
#define TFQT 4		// Translation: Intermediate Scattering Function
#define TCHI 5		// Translation: Four Point Correlation
#define TCLS 6		// Translation: Highly Mobile Clusters
#define FAST 7          // Rot and Trans: MSD Coupling
#define TSHL 8          // Translation: Shell Dynamics
#define VISC 9          // Stress Tensor

/**************************************************************************************************************/
/*** FUNCTIONS: FRAME_SETM / FRAME_FREE / FRAME_READ / FRAME_NORM / FRAME_ANAL ***/

double	OFF;
double	WAV;
int 	MOB;
double	CUT;
int     PPM;
double	SHL;

int     file_mols;
FILE    *qfile;
#if ANALYSIS==FAST || ANALYSIS==TSHL
FILE    *shfile;
#endif

#if ANALYSIS==CORR
  #define NQ 1
  #include "frame_corr_.h"
#elif ANALYSIS==RCOR
  #define NQ 3
  #include "frame_rcor_.h"
#elif ANALYSIS==RMSD
  #define NQ 3
  #include "frame_rmsd_.h"
#elif ANALYSIS==TMSD
  #define NQ 2
  #include "frame_tmsd_.h"
#elif ANALYSIS==TFQT
  #define NQ 3
  #include "frame_tmsd_.h"
#elif ANALYSIS==TCHI
  #define NQ 5
  #include "frame_tmsd_.h"
#elif ANALYSIS==TCLS
  #define NQ 7
  #include "frame_tmsd_.h"
#elif ANALYSIS==FAST
  #define NQ 5
  #include "frame_fast_.h"
#elif ANALYSIS==TSHL
  #define NQ 4
  #include "frame_tshl_.h"
#elif ANALYSIS==VISC
  #define NQ 3
  #include "frame_visc_.h"
#endif

/**************************************************************************************************************/
/*** FUNCTION: FRAME_HEAD --- Reads through the input file and sets the following variables: ***/

long int file_end;     // file  size in bytes
long int file_head;    // head  size in bytes
long int file_bytes;   // frame size in bytes
long int file_frames;  // file  size in frames

#if ANALYSIS==FAST || ANALYSIS==TSHL
long int shfile_head;  // shell head  size in bytes
long int shfile_bytes; // shell frame size in bytes
#endif

void frame_head()
{

  fseek(qfile,0,SEEK_END);
  file_end=ftell(qfile);
  rewind(qfile);

  #if FORMAT<=XYZ
    file_mols=1;
    file_head=0;
    fscanf(qfile, "%*[^\n]%*c");
    file_bytes=ftell(qfile);
  #elif FORMAT<=ATM
    fseek(qfile,29,SEEK_SET);
    fscanf(qfile,"%d",&file_mols);
    file_head=1;
    while(file_mols/(int)pow(10,file_head))
      file_head++;
    file_head+=29+31+1;
    file_bytes=(3*8+1);                                                                  //box
    file_bytes+=((3*file_mols)/10)*(10*8+1);                                             //full
    file_bytes+=((3*file_mols)%10)*8;                                                    //extra mols
    file_bytes+=(long int)ceil((double)((3*file_mols)%10)/(double)((3*file_mols)%10+1)); //extra newline
  #elif FORMAT==G96
    fscanf(qfile,"%*[^\n]%*c");
    fscanf(qfile,"%*[^\n]%*c");
    fscanf(qfile,"%*[^\n]%*c");
    file_head=ftell(qfile);
    fscanf(qfile,"%*[^D]%*c%*c");
    fscanf(qfile,"%*[^D]%*c%*c");
    fscanf(qfile,"%*[^D]%*c%*c");
    fscanf(qfile,"%*[^D]%*c%*c");
    fscanf(qfile,"%*[^D]%*c%*c");
    fscanf(qfile,"%*[^D]%*c%*c");
    file_bytes=ftell(qfile)-file_head;
    file_mols=file_bytes-(9+31+4+12+4+12+4+4+4);
    file_mols/=(2*(3*15+1));
    fseek(qfile,file_head,SEEK_SET);
  #else
    exit(0);
  #endif

  #if FORMAT==DIM
    file_mols/=2;
  #endif

  file_frames=(file_end-file_head)/file_bytes;

  #if ANALYSIS==TSHL
    shfile_head=61+1;
    shfile_bytes=(6*8+1)+(3*8+1);
    fseek(shfile,62,SEEK_SET);
  #endif

}

/**************************************************************************************************************/
/*** FUNCTION: CHECK_FLAGS --- Scans therough the command line flags and sets the following: ***/

int 	INFILE;
int 	SHFILE;
int 	INIT;
int 	CHOP;
int 	NUMBER;
int 	TCYCLE;
int 	NCYCLE;
int 	DECIMALS;
int 	NOHEADER;
int 	QUIET;

void check_flags(int argc, char **argv){

  INFILE	= 0;		// WHATEVER
  SHFILE        = 0;		// WHATEVER

  INIT		= 0;		// WHATEVER
  CHOP		= 0;		// WHATEVER

  NUMBER	= 200;		// > 1

  TCYCLE	= 0;		// 0 <= T <= 31
  NCYCLE	= 20;		// 1 <= N <= INFTY

  OFF		= 0.00;		// FOR CORR:     OFFSET
  WAV		= 7.25;		// FOR FQT:      WAVE NUMBER
  MOB		= 192;		// FOR CLUSTERS: NUMBER MOBILE
  CUT		= 1.40;		// FOR CLUSTERS: CUTOFF DISTANCE SQ
  PPM           = 20;           // FOR SHELL:    FRAMERATE RATIO PROBE / A's
  SHL		= 5.76;		// FOR SHELL:    CUTOFF DISTANCE SQ

  DECIMALS	= 6;		// 0 <= D <= INF
  NOHEADER	= 0;		// 0 <= P <= 1
  QUIET		= 0;		// 0 <= Q <= 1
  
  for(int n=1;n<argc;n++)
  {
    if(strcmp(argv[n],"-f")==0){
      INFILE=n+1;
	  n++;
    }
    else if(strcmp(argv[n],"-s")==0){
      SHFILE=n+1;
	  n++;
    }
    else if(strcmp(argv[n],"-init")==0){
      INIT=atoi(argv[n+1]);
	  n++;
    }
    else if(strcmp(argv[n],"-chop")==0){
      CHOP=atoi(argv[n+1]);
      n++;
    }
    else if(strcmp(argv[n],"-n")==0){
      NUMBER=atoi(argv[n+1]);
	  n++;
    }
    else if(strcmp(argv[n],"-tc")==0){
      TCYCLE=atoi(argv[n+1]);
	  if(TCYCLE>31)
	    TCYCLE=31;
	  n++;
    }
    else if(strcmp(argv[n],"-nc")==0){
      NCYCLE=atoi(argv[n+1]);
      n++;
    }
    else if(strcmp(argv[n],"-off")==0){
      OFF=atof(argv[n+1]);
	  n++;
    }
    else if(strcmp(argv[n],"-wav")==0){
      WAV=atof(argv[n+1]);
	  n++;
    }
    else if(strcmp(argv[n],"-mob")==0){
      MOB=atoi(argv[n+1]);
	  n++;
    }
    else if(strcmp(argv[n],"-cut")==0){
      CUT=atof(argv[n+1]);
	  n++;
    }
    else if(strcmp(argv[n],"-ppm")==0){
      PPM=atoi(argv[n+1]);
	  n++;
    }
    else if(strcmp(argv[n],"-shl")==0){
      SHL=atof(argv[n+1]);
	  n++;
    }
    else if(strcmp(argv[n],"-d")==0){
      DECIMALS=1;
	  n++;
    }
    else if(strcmp(argv[n],"-h")==0){
      NOHEADER=1;
    }
    else if(strcmp(argv[n],"-q")==0){
	  QUIET=1;
    }
    else{
      printf("\n  Qu'est-ce que c'est: %s\n\n",argv[n]);
      exit(0);
    }

  }

 if(INFILE==0){
    printf("\n  Ah, merde!\n\n");
    exit(0);
  }
  #if ANALYSIS==TSHL
  if(SHFILE==0){
    printf("\n  Ah, merde!\n\n");
    exit(0);
  } 
  #endif

  return;
}






