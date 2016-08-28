/**************************************************************************************************************/
/*** FAST: DIM ***/

double (*A)[3]; // current angular coordinates
double (*B)[3]; // previous angular coordinates
double (*C)[3]; // cumulative arc length displacement
double (*R)[3]; // current atomic coordinates
double (*S)[3]; // previous atomic coordinates
double (*T)[3]; // cumulative atomic displacement

double (*pos)[2][3];
double (*sep)[3];

double AVG[NQ];

void frame_setm()
{
  A=(double*)malloc(3*file_mols*sizeof(double));
  B=(double*)malloc(3*file_mols*sizeof(double));
  C=(double*)malloc(3*file_mols*sizeof(double));
  R=(double*)malloc(3*file_mols*sizeof(double));
  S=(double*)malloc(3*file_mols*sizeof(double));
  T=(double*)malloc(3*file_mols*sizeof(double));
  pos=(double*)malloc(6*file_mols*sizeof(double));
  sep=(double*)malloc(3*file_mols*sizeof(double));
}

void frame_free()
{
  free(A);
  free(B);
  free(C);
  free(R);
  free(S);
  free(T);
  free(pos);
  free(sep);
}

double box[3];

void frame_read()
{

  for(int i=0;i<file_mols;i++)
    for(int j=0;j<6;j++)
      fscanf(qfile,"%lf",&pos[i][j/3][j%3]);
  fscanf(qfile,"%lf %lf %lf",&box[0],&box[1],&box[2]);
 
  for(int i=0;i<file_mols;i++){

    sep[i][0]=pos[i][0][0]-pos[i][1][0]-box[0]*floor((pos[i][0][0]-pos[i][1][0])/box[0]+0.5);
    sep[i][1]=pos[i][0][1]-pos[i][1][1]-box[1]*floor((pos[i][0][1]-pos[i][1][1])/box[1]+0.5);     
    sep[i][2]=pos[i][0][2]-pos[i][1][2]-box[2]*floor((pos[i][0][2]-pos[i][1][2])/box[2]+0.5);

    A[i][0]=atan2(sep[i][1],sep[i][0]);
    A[i][1]=atan2(sep[i][2],sep[i][0]);
    A[i][2]=atan2(sep[i][2],sep[i][1]);

    R[i][0]=pos[i][1][0]+sep[i][0]/2.0;
    R[i][1]=pos[i][1][1]+sep[i][1]/2.0;
    R[i][2]=pos[i][1][2]+sep[i][2]/2.0;

    if(R[i][0]>box[0])
      R[i][0]-=box[0];
    if(R[i][0]<0.0)
      R[i][0]+=box[0];

    if(R[i][1]>box[1])
      R[i][1]-=box[1];
    if(R[i][1]<0.0)
      R[i][1]+=box[1];

    if(R[i][2]>box[2])
      R[i][2]-=box[2];
    if(R[i][2]<0.0)
      R[i][2]+=box[2];

  }

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
    T[i][0]=0.0;
    T[i][1]=0.0;
    T[i][2]=0.0;
    S[i][0]=R[i][0];
    S[i][1]=R[i][1];
    S[i][2]=R[i][2];
  }

  for(int q=0;q<NQ;q++)
    result[q]=(double)file_mols;

}

void frame_anal(double *result)
{

  const double PI = 3.14159265;

  for(int q=0;q<NQ;q++)
    result[q]=0.0;

  for(int i=0;i<file_mols;i++){

    double arc[3];

    arc[0]=(A[i][0]-B[i][0])-2.0*PI*floor((A[i][0]-B[i][0])/(2.0*PI)+0.5);
    arc[1]=(A[i][1]-B[i][1])-2.0*PI*floor((A[i][1]-B[i][1])/(2.0*PI)+0.5);
    arc[2]=(A[i][2]-B[i][2])-2.0*PI*floor((A[i][2]-B[i][2])/(2.0*PI)+0.5);

    arc[0]*=sqrt(sep[i][0]*sep[i][0]+sep[i][1]*sep[i][1]);
    arc[1]*=sqrt(sep[i][0]*sep[i][0]+sep[i][2]*sep[i][2]);
    arc[2]*=sqrt(sep[i][1]*sep[i][1]+sep[i][2]*sep[i][2]);

    C[i][0]+=arc[0];
    C[i][1]+=arc[1];
    C[i][2]+=arc[2];

    B[i][0]=A[i][0];
    B[i][1]=A[i][1];
    B[i][2]=A[i][2];

    T[i][0]+=(R[i][0]-S[i][0])-box[0]*floor((R[i][0]-S[i][0])/box[0]+0.5);
    T[i][1]+=(R[i][1]-S[i][1])-box[1]*floor((R[i][1]-S[i][1])/box[1]+0.5);
    T[i][2]+=(R[i][2]-S[i][2])-box[2]*floor((R[i][2]-S[i][2])/box[2]+0.5);

    S[i][0]=R[i][0];
    S[i][1]=R[i][1];
    S[i][2]=R[i][2];

    double dr = sqrt(T[i][0]*T[i][0]+T[i][1]*T[i][1]+T[i][2]*T[i][2]);
    double da = sqrt(C[i][0]*C[i][0]+C[i][1]*C[i][1]+C[i][2]*C[i][2]);

    result[0] += fabs(T[i][0]) - AVG[0];
    result[1] += fabs(T[i][1]) - AVG[1];
    result[2] += fabs(T[i][2]) - AVG[2];
    result[3] += dr            - AVG[3];

    result[4] += fabs(C[i][0]) - AVG[4];
    result[5] += fabs(C[i][1]) - AVG[5];
    result[6] += fabs(C[i][2]) - AVG[6];
    result[7] += da            - AVG[7];

    result[8]  += pow((fabs(T[i][0])-AVG[0]),2);
    result[9]  += pow((fabs(T[i][1])-AVG[1]),2);
    result[10] += pow((fabs(T[i][2])-AVG[2]),2);
    result[11] += pow(dr-AVG[3],2);

    result[12] += pow((fabs(C[i][0])-AVG[4]),2);
    result[13] += pow((fabs(C[i][1])-AVG[5]),2);
    result[14] += pow((fabs(C[i][2])-AVG[6]),2);
    result[15] += pow(da-AVG[7],2);

    result[16] += (fabs(T[i][0])-AVG[0])*(fabs(C[i][2])-AVG[6]);
    result[17] += (fabs(T[i][1])-AVG[1])*(fabs(C[i][1])-AVG[5]);
    result[18] += (fabs(T[i][2])-AVG[2])*(fabs(C[i][0])-AVG[4]);
    result[19] += (dr-AVG[3])*(da-AVG[7]);

  }

}
























