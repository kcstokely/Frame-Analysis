/**************************************************************************************************************/
/*** RCOR: XYZ, DIM ***/

double (*P)[3];  // current  coordinates
double (*Q)[3];  // previous coordinates
double *D;       // current  dichroism
double *E;       // previous dichroism

#if FORMAT==DIM
  double (*pos)[2][3];
#endif

void frame_setm()
{
  P=(double*)malloc(3*file_mols*sizeof(double));
  Q=(double*)malloc(3*file_mols*sizeof(double));
  D=(double*)malloc(file_mols*sizeof(double));
  E=(double*)malloc(file_mols*sizeof(double));
  #if FORMAT==DIM
    pos=(double*)malloc(6*file_mols*sizeof(double));
  #endif
}

void frame_free()
{
  free(P);
  free(Q);
  free(D);
  free(E);
  #if FORMAT==DIM
    free(pos);
  #endif
}

void frame_read()
{
  #if FORMAT==XYZ

  for(int i=0;i<file_mols;i++){
    for(int j=0;j<3;j++)
      fscanf(qfile,"%lf",&P[i][j]);
    D[i]=(P[i][0]*P[i][0]-P[i][1]*P[i][1])/(P[i][0]*P[i][0]+P[i][1]*P[i][1]);
  }

  #elif FORMAT==DIM

  double box[3];
  for(int i=0;i<file_mols;i++)
    for(int j=0;j<6;j++)
      fscanf(qfile,"%lf",&pos[i][j/3][j%3]);
  fscanf(qfile,"%lf %lf %lf",&box[0],&box[1],&box[2]);
  for(int i=0;i<file_mols;i++){
    P[i][0]=pos[i][0][0]-pos[i][1][0]-box[0]*floor((pos[i][0][0]-pos[i][1][0])/box[0]+0.5);
    P[i][1]=pos[i][0][1]-pos[i][1][1]-box[1]*floor((pos[i][0][1]-pos[i][1][1])/box[1]+0.5);     
    P[i][2]=pos[i][0][2]-pos[i][1][2]-box[2]*floor((pos[i][0][2]-pos[i][1][2])/box[2]+0.5);
    D[i]=(P[i][0]*P[i][0]-P[i][1]*P[i][1])/(P[i][0]*P[i][0]+P[i][1]*P[i][1]);
  }

  #endif
}

void frame_norm(double *result)
{
  result[0]=(double)file_mols;
  result[1]=(double)file_mols;
  result[2]=0.0;
  for(int i=0;i<file_mols;i++){
    result[2]+=D[i]*D[i];
    Q[i][0]=P[i][0];
    Q[i][1]=P[i][1];
    Q[i][2]=P[i][2];
    E[i]=D[i];
  }
}

void frame_anal(double *result)
{
  result[0]=0.0;
  result[1]=0.0;
  result[2]=0.0;
  for(int i=0;i<file_mols;i++){
    double f=P[i][0]*Q[i][0]+P[i][1]*Q[i][1]+P[i][2]*Q[i][2];
    f/=sqrt(P[i][0]*P[i][0]+P[i][1]*P[i][1]+P[i][2]*P[i][2]);
    f/=sqrt(Q[i][0]*Q[i][0]+Q[i][1]*Q[i][1]+Q[i][2]*Q[i][2]);
    result[0]+=f;
    result[1]+=0.5*(3.0*f*f-1.0);
    result[2]+=D[i]*E[i];
  }
}




















