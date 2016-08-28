/**************************************************************************************************************/
/*** FAST: DIM ***/

double (*A)[3]; // current  angular coordinates
double (*B)[3]; // original angular coordinates
double (*C);
double (*R)[3]; // current  atomic  coordinates
double (*S)[3]; // previous atomic  coordinates
double (*T)[3]; // cumulative atomic displacement
double (*D);
double (*flag);

double (*pos)[2][3];

double AVG[NQ];

void frame_setm()
{
  A=(double*)malloc(3*file_mols*sizeof(double));
  B=(double*)malloc(3*file_mols*sizeof(double));
  C=(double*)malloc(1*file_mols*sizeof(double));
  R=(double*)malloc(3*file_mols*sizeof(double));
  S=(double*)malloc(3*file_mols*sizeof(double));
  T=(double*)malloc(3*file_mols*sizeof(double));
  D=(double*)malloc(1*file_mols*sizeof(double));
  flag=(double*)malloc(1*file_mols*sizeof(double));
  pos=(double*)malloc(6*file_mols*sizeof(double));
}

void frame_free()
{
  free(A);
  free(B);
  free(C);
  free(R);
  free(S);
  free(T);
  free(D);
  free(flag);
  free(pos);
}

double box[3];

void frame_read()
{

  for(int i=0;i<file_mols;i++)
    for(int j=0;j<6;j++)
      fscanf(qfile,"%lf",&pos[i][j/3][j%3]);
  fscanf(qfile,"%lf %lf %lf",&box[0],&box[1],&box[2]);
 
  for(int i=0;i<file_mols;i++){

    double sep[3];

    sep[0]=pos[i][0][0]-pos[i][1][0]-box[0]*floor((pos[i][0][0]-pos[i][1][0])/box[0]+0.5);
    sep[1]=pos[i][0][1]-pos[i][1][1]-box[1]*floor((pos[i][0][1]-pos[i][1][1])/box[1]+0.5);     
    sep[2]=pos[i][0][2]-pos[i][1][2]-box[2]*floor((pos[i][0][2]-pos[i][1][2])/box[2]+0.5);

    A[i][0]=sep[0];
    A[i][1]=sep[1];
    A[i][2]=sep[2];

    R[i][0]=pos[i][1][0]+sep[0]/2.0;
    R[i][1]=pos[i][1][1]+sep[1]/2.0;
    R[i][2]=pos[i][1][2]+sep[2]/2.0;

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

  for(int q=0;q<5;q++)
    result[q]=(double)file_mols;

  for(int q=5;q<10;q++)
    result[q]=(double)MOB;

}

void frame_anal(double *result)
{

  const double PI = 3.14159265;

  for(int q=0;q<NQ;q++)
    result[q]=0.0;

  int    cut_ind;
  double cut_val = 1000000.0;

  for(int i=0;i<file_mols;i++){

    C[i]  = A[i][0]*B[i][0]+A[i][1]*B[i][1]+A[i][2]*B[i][2];
    C[i] /= sqrt(A[i][0]*A[i][0]+A[i][1]*A[i][1]+A[i][2]*A[i][2]);
    C[i] /= sqrt(B[i][0]*B[i][0]+B[i][1]*B[i][1]+B[i][2]*B[i][2]);

    if( 1.000-C < 0.0001 )
      C[i] = 0.0;
    else
      C[i]  = acos(C[i]);

    T[i][0]+=(R[i][0]-S[i][0])-box[0]*floor((R[i][0]-S[i][0])/box[0]+0.5);
    T[i][1]+=(R[i][1]-S[i][1])-box[1]*floor((R[i][1]-S[i][1])/box[1]+0.5);
    T[i][2]+=(R[i][2]-S[i][2])-box[2]*floor((R[i][2]-S[i][2])/box[2]+0.5);

    S[i][0]=R[i][0];
    S[i][1]=R[i][1];
    S[i][2]=R[i][2];

    D[i] = sqrt(T[i][0]*T[i][0]+T[i][1]*T[i][1]+T[i][2]*T[i][2]);

    if(i<MOB){
      flag[i]=1;
      if(D[i]<cut_val){
        cut_ind=i;
        cut_val=D[i];
      }
    }
    else{
      flag[i]=0;
      if(D[i]>cutoff){
        flag[i]=1;
        flag[cut_ind]=0;
        cut_val=D[i];
        for(j=0;j<i;j++){
          if(flag[j]==1 && D[j]<cut_val){
            cut_ind=j;
            cut_val=D[j];
          }
        }
      }
    }

    result[0] += D[i]-AVG[0];
    result[1] += C-AVG[1];
    result[2] += pow(D[i]-AVG[0],2);
    result[3] += pow(C-AVG[1],2);
    result[4] += (D[i]-AVG[0])*(C-AVG[1]);

  }

  for(i=0;i<file_mols;i++)
    if(flag[i]==1){
      result[5] += D[i]-AVG[0];
      result[6] += C-AVG[1];
      result[7] += pow(D[i]-AVG[0],2);
      result[8] += pow(C-AVG[1],2);
      result[9] += (D[i]-AVG[0])*(C-AVG[1]);
    }

}
























