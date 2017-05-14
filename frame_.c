
#include "frame_.h"

int main(int argc, char *argv[]){


  /***************************************************************************/
  /*** SET ANALYSIS OPTIONS ***/

    check_flags(argc,argv);

  /***************************************************************************/
  /*** OPEN INPUT FILE(S) ***/

    qfile=fopen(argv[INFILE],"r");
    if(qfile==NULL){
      printf("\n  Ah, merde!\n\n");
      exit(0);
    }

    #if ANALYSIS==TSHL
    shfile=fopen(argv[SHFILE],"r");
    if(shfile==NULL){
      printf("\n  Ah, merde!\n\n");
      exit(0);
    }
    #endif

    #if ANALYSIS==FAST
    int tmp;
    double (*ref)[3];
    if(SHFILE>0){
      ref=(double*)malloc(NQ*(TCYCLE+NCYCLE)*sizeof(double));
      shfile=fopen(argv[SHFILE],"r");
      if(shfile==NULL){
        printf("\n  Ah, merde!\n\n");
        exit(0);
      }
      for(int i=0;i<7;i++)
        fscanf(qfile,"%*[^\n]%*c");
      for(int delta=0;delta<TCYCLE+NCYCLE;delta++){
        fscanf(shfile,"%d",&tmp);
        for(int q=0;q<NQ;q++)
          fscanf(shfile,"%lf",&ref[delta][q]);
        fscanf(shfile,"%d",&tmp);
      }
      fclose(shfile);
    }
    else{
      for(int q=0;q<NQ;q++)
        AVG[q]=0.0;
    }
    #endif

  /***************************************************************************/
  /*** READ INPUT FILE INFO ***/

    frame_head();
    frame_setm();

  /***************************************************************************/
  /*** SETUP ANALYSIS ***/

    int i_init = INIT;
    int i_skip;

    if((file_frames-1)>(long int)(2*(long int)pow(2,TCYCLE)*(long int)NCYCLE))
      i_skip=(file_frames-1-((long int)pow(2,TCYCLE)*(long int)NCYCLE))/(long int)(NUMBER);
    else
      i_skip=(file_frames-1)/(2*(NUMBER));

    if(i_skip==0)
      i_skip=1;

    if(CHOP>0){
      file_frames = CHOP;
      file_end = file_head+file_bytes*(file_frames);
      i_skip=1;
      i_init=0;
    }

  /***************************************************************************/
  /*** SETUP SCREEN ***/

    int g=1;
    while((int)pow(2,TCYCLE)/(int)pow(10,g))
      g++;
    int h=1;
    while((int)pow(2,TCYCLE)*NCYCLE/(int)pow(10,h))
      h++;
    g+=h+4;

    if(QUIET==0){
      printf("\n");
      printf("    %d %d | %d : %d : %d\n",FORMAT,ANALYSIS,1,(int)pow(2,TCYCLE),(int)pow(2,TCYCLE)*NCYCLE);
      printf("\n");
      printf("    cadres:  %*ld\n\n",g,file_frames);
      printf("    debut:   %*d\n",g,i_init);
      printf("    ecart:   %*d\n",g,i_skip);
      printf("    fin:     %*d\n\n",g,i_init+(NUMBER-1)*i_skip);
      printf("    nombre:  %*d\n\n",g,NUMBER);
    }

  /***************************************************************************/
  /*** INITIALIZE STATS VARIABLES ***/

    double value_nrm[NQ];
    double value_avg[NQ];
    double (*array_nrm)[NQ]  = (double*)calloc(NQ*(TCYCLE+NCYCLE),sizeof(double));
    double (*array_avg)[NQ]  = (double*)calloc(NQ*(TCYCLE+NCYCLE),sizeof(double));
    int    *array_count      = (int*)calloc((TCYCLE+NCYCLE),sizeof(int));

  /***************************************************************************/
  /*** MAIN LOOP ***/

    /*** INITIAL ***************************************************/
    for(int initial=0; initial<NUMBER; initial++){

      // update screen
      if(QUIET==0){
        printf("\r    courant: %*d",g,initial+1);
        fflush(stdout);
      }

      // seek to initial frame
      if(file_head+file_bytes*(long int)(i_init+initial*i_skip+2)>file_end)
        break;
      fseek(qfile,file_head+file_bytes*(long int)(i_init+initial*i_skip),SEEK_SET);
      #if ANALYSIS==TSHL
      //assumes twenty P frames per A frame
      fseek(shfile,shfile_head+(PPM)*shfile_bytes*(long int)(i_init+initial*i_skip),SEEK_SET);
      #endif

      // calculate norms
      frame_read();
      frame_norm(value_nrm);

      /*** FINAL ***************************************************/
      int e=0;
      for(int delta=0; delta<TCYCLE+NCYCLE; delta++){

        #if ANALYSIS==FAST
        if(SHFILE>0){
          for(int q=0; q<NQ; q++)
            AVG[q]=ref[delta][q];
        }
        #endif

        // seek to final frame
        int d = (delta<TCYCLE) ? (int)pow(2,delta):(delta-TCYCLE+1)*(int)pow(2,TCYCLE);
        if(file_head+file_bytes*(long int)(i_init+initial*i_skip+d+1) > file_end)
          break;
        fseek(qfile,file_bytes*(long int)(d-e-1),SEEK_CUR);
        e=d;
  
        // calculate values
        frame_read();
        frame_anal(value_avg);

        // update running totals
        array_count[delta]++;
        for(int q=0; q<NQ; q++){
          array_nrm[delta][q]+=value_nrm[q];
          array_avg[delta][q]+=value_avg[q];
        }

      }

    }

  /***************************************************************************/
  /*** FINISH CALCULATING AVERAGES ***/

  for(int delta=0; delta<TCYCLE+NCYCLE; delta++)
    for(int q=0; q<NQ; q++)
      array_avg[delta][q]/=array_nrm[delta][q];

  /***************************************************************************/
  /*** CLOSE INPUT FILE ***/

  frame_free();

  fclose(qfile);
  #if ANALYSIS==TSHL
  fclose(shfile);
  #endif

  #if ANALYSIS==FAST
  if(SHFILE>0)
    free(ref);
  #endif

  /***************************************************************************/
  /*** OUTPUT ***/

  char fname[100];
  strcpy(fname,argv[INFILE]);
  #if ANALYSIS==CORR
    strcat(fname,".corr");
  #elif ANALYSIS==RCOR
    strcat(fname,".rcor");
  #elif ANALYSIS==RMSD
    strcat(fname,".rmsd");
  #elif ANALYSIS==TMSD
    strcat(fname,".tmsd");
  #elif ANALYSIS==TFQT
    strcat(fname,".tfqt");
  #elif ANALYSIS==TCHI
    strcat(fname,".tchi");
  #elif ANALYSIS==TCLS
    strcat(fname,".tcls");
  #elif ANALYSIS==FAST
    strcat(fname,".fnew");
  #elif ANALYSIS==TSHL
    strcat(fname,".tshl");
  #elif ANALYSIS==VISC
    strcat(fname,".visc");
  #endif
  strcat(fname,".avg");

  qfile=fopen(fname,"w");

  if(NOHEADER==0){
    fprintf(qfile,"#\n");
    fprintf(qfile,"# INIT\t\t%d\n",INIT);
    fprintf(qfile,"# CHOP\t\t%d\n",CHOP);
    fprintf(qfile,"# NUMBER\t%d\n",NUMBER);
    fprintf(qfile,"# TCYCLE\t%d\n",TCYCLE);
    fprintf(qfile,"# NCYCLE\t%d\n",NCYCLE);
    #if ANALYSIS==CORR
    fprintf(qfile,"# OFFSET\t%lf\n",OFF);
    #elif ANALYSIS>=TFQT && ANALYSIS<=TCLS
    fprintf(qfile,"# WAV\t\t%lf\n",WAV);
    #if ANALYSIS==TCLS
    fprintf(qfile,"# MOB\t\t%d\n", MOB);
    fprintf(qfile,"# CUT\t\t%lf\n",CUT);
    #endif
    #elif ANALYSIS==TSHL
    fprintf(qfile,"# PPM\t\t%d\n", PPM);
    fprintf(qfile,"# SHL\t\t%lf\n",SHL);
    #endif
    fprintf(qfile,"#\n");
  }

  for(int delta=0; delta<TCYCLE+NCYCLE; delta++){
    if(array_count[delta]>0){
      fprintf(qfile,"%10d",(delta<TCYCLE)?(int)pow(2,delta):(delta-TCYCLE+1)*(int)pow(2,TCYCLE));
      for(int q=0; q<NQ; q++)
        fprintf(qfile," %.*lf",DECIMALS,array_avg[delta][q]);
      fprintf(qfile," %10d\n",array_count[delta]);
    }
  }

  fclose(qfile);

  /***************************************************************************/
  /*** FREE STATS VARIABLES ***/

  free(array_nrm);
  free(array_avg);

  /***************************************************************************/
  /*** CLEAR SCREEN ***/

  if(QUIET==0)
    printf("\n\n");

  /***************************************************************************/
  /*** DONE ***/

  return 0;
}



























