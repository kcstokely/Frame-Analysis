/**************************************************************************************************************/
/*** RMSD: XYZ, DIM  ***/

double (*A)[3];  // current  angle
double (*B)[3];  // previous angle
double (*C)[3];  // cumulative angular displacement

#if FORMAT==DIM
  double (*pos)[2][3];
#endif

void frame_setm()
{
  A=(double*)malloc(3*file_mols*sizeof(double));
  B=(double*)malloc(3*file_mols*sizeof(double));
  C=(double*)malloc(3*file_mols*sizeof(double));
  #if FORMAT==DIM
    pos=(double*)malloc(6*file_mols*sizeof(double));
  #endif
}

void frame_free()
{
  free(A);
  free(B);
  free(C);
  #if FORMAT==DIM
    free(pos);
  #endif
}


void frame_read()
{

  double sep[3];

  #if FORMAT==XYZ

  for(int i=0;i<file_mols;i++){
    for(int j=0;j<3;j++)
      fscanf(qfile,"%lf",&sep[j]);
    A[i][0]=atan2(sep[1],sep[0]);
    A[i][1]=atan2(sep[2],sep[0]);
    A[i][2]=atan2(sep[2],sep[1]);
  }

  #elif FORMAT==DIM

  double box[3];
  for(int i=0;i<file_mols;i++)
    for(int j=0;j<6;j++)
      fscanf(qfile,"%lf",&pos[i][j/3][j%3]);
  fscanf(qfile,"%lf %lf %lf",&box[0],&box[1],&box[2]);

  for(int i=0;i<file_mols;i++){

    for(int j=0;j<2;j++){
      for(int k=0;k<3;k++){
        if(pos[i][j][k]<0.0)
          pos[i][j][k]+=box[k];
        if(pos[i][j][k]>=box[k])
          pos[i][j][k]-=box[k];
      }
    }

    sep[0]=pos[i][0][0]-pos[i][1][0]-box[0]*floor((pos[i][0][0]-pos[i][1][0])/box[0]+0.5);
    sep[1]=pos[i][0][1]-pos[i][1][1]-box[1]*floor((pos[i][0][1]-pos[i][1][1])/box[1]+0.5);     
    sep[2]=pos[i][0][2]-pos[i][1][2]-box[2]*floor((pos[i][0][2]-pos[i][1][2])/box[2]+0.5);

    A[i][0]=atan2(sep[1],sep[0]);
    A[i][1]=atan2(sep[2],sep[0]);
    A[i][2]=atan2(sep[2],sep[1]);

  }

  #endif

}

void frame_norm(double *result)
{
  for(int i=0;i<file_mols;i++){
    C[i][0]=0.0;
    C[i][1]=0.0;
    C[i][2]=0.0;
    B[i][0]=A[i][0];
    B[i][1]=A[i][1];
    B[i][2]=A[i][2];
  }
  result[0]=(double)file_mols;
  result[1]=(double)file_mols;
  result[2]=(double)file_mols;
}

void frame_anal(double *result)
{
  const double PI = 3.14159265;
  result[0]=0.0;
  result[1]=0.0;
  result[2]=0.0;
  for(int i=0;i<file_mols;i++){

    C[i][0]+=(A[i][0]-B[i][0])-2.0*PI*floor((A[i][0]-B[i][0])/(2.0*PI)+0.5);
    C[i][1]+=(A[i][1]-B[i][1])-2.0*PI*floor((A[i][1]-B[i][1])/(2.0*PI)+0.5);
    C[i][2]+=(A[i][2]-B[i][2])-2.0*PI*floor((A[i][2]-B[i][2])/(2.0*PI)+0.5);

    //double xy = (A[i][0]-B[i][0])-2.0*PI*floor((A[i][0]-B[i][0])/(2.0*PI)+0.5);
    //double xz = (A[i][1]-B[i][1])-2.0*PI*floor((A[i][1]-B[i][1])/(2.0*PI)+0.5);
    //double yz = (A[i][2]-B[i][2])-2.0*PI*floor((A[i][2]-B[i][2])/(2.0*PI)+0.5);

    //double magsq = sep[0]*sep[0]+sep[1]*sep[1]+sep[2]*sep[2];

    //xy *= 0.5*sqrt(magsq-sep[2]*sep[2]);
    //xz *= 0.5*sqrt(magsq-sep[1]*sep[1]);
    //yz *= 0.5*sqrt(magsq-sep[0]*sep[0]);

    //C[i][0]+=xy;
    //C[i][1]+=xz;
    //C[i][2]+=yz;

    result[0]+=C[i][0]*C[i][0];
    result[1]+=C[i][1]*C[i][1];
    result[2]+=C[i][2]*C[i][2];
    B[i][0]=A[i][0];
    B[i][1]=A[i][1];
    B[i][2]=A[i][2];
  }
}













