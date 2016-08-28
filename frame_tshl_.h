/**************************************************************************************************************/
/*** Wat  ***/

double (*R)[3]; // current positions
double (*S)[3]; // previous positions
double *dx;     // cumulative displacement
double *dy;     // cumulative displacement
double *dz;     // cumulative displacement

int    *flag;

void frame_setm()
{
  R= (double*)malloc(3*file_mols*sizeof(double));
  S= (double*)malloc(3*file_mols*sizeof(double));
  dx=(double*)malloc(file_mols*sizeof(double));
  dy=(double*)malloc(file_mols*sizeof(double));
  dz=(double*)malloc(file_mols*sizeof(double));
  flag= (int*)malloc(file_mols*sizeof(int));
}

void frame_free()
{
  free(R);
  free(S);
  free(dx);
  free(dy);
  free(dz);
  free(flag);
}

double box[3];

void frame_read()
{
  for(int i=0;i<3*file_mols;i++)
    fscanf(qfile,"%lf",&R[i/3][i%3]);
  fscanf(qfile,"%lf%lf%lf",&box[0],&box[1],&box[2]);
}

void frame_norm(double *result)
{

  // read in probe coords
  double pos[2][3];
  for(int i=0;i<6;i++)
    fscanf(shfile,"%lf",&pos[i/3][i%3]);
  fscanf(shfile,"%lf%lf%lf",&box[0],&box[1],&box[2]);

  // get orientation
  double sep[3];
  sep[0]=pos[0][0]-pos[1][0]-box[0]*floor((pos[0][0]-pos[1][0])/box[0]+0.5);
  sep[1]=pos[0][1]-pos[1][1]-box[1]*floor((pos[0][1]-pos[1][1])/box[1]+0.5);     
  sep[2]=pos[0][2]-pos[1][2]-box[2]*floor((pos[0][2]-pos[1][2])/box[2]+0.5);

  // get center of mass
  double com[3];
  com[0]=pos[1][0]+sep[0]/2.0;
  com[1]=pos[1][1]+sep[1]/2.0;
  com[2]=pos[1][2]+sep[2]/2.0;
  if(com[0]>=box[0])
    com[0]-=box[0];
  if(com[0]<0.0)
    com[0]+=box[0];
  if(com[1]>=box[1])
    com[1]-=box[1];
  if(com[1]<0.0)
    com[1]+=box[1];
  if(com[2]>=box[2])
    com[2]-=box[2];
  if(com[2]<0.0)
    com[2]+=box[2];

  // flag shell molecules
  int shl_mols = 0;
  for(int i=0;i<file_mols;i++){
    dx[i]=0.0;
    dy[i]=0.0;
    dz[i]=0.0;
    double x=(R[i][0]-com[0])-box[0]*floor((R[i][0]-com[0])/box[0]+0.5);
    double y=(R[i][1]-com[1])-box[1]*floor((R[i][1]-com[1])/box[1]+0.5);
    double z=(R[i][2]-com[2])-box[2]*floor((R[i][2]-com[2])/box[2]+0.5);
    if(x*x+y*y+z*z<SHL){
      shl_mols++;
      flag[i]=1;
    }
    else
      flag[i]=0;
  }

  for(int i=0;i<file_mols;i++){
    dx[i]=0.0;
    dy[i]=0.0;
    dz[i]=0.0;
    S[i][0]=R[i][0];
    S[i][1]=R[i][1];
    S[i][2]=R[i][2];
  }

  result[0]=(double)shl_mols;
  result[1]=(double)shl_mols;
  result[2]=(double)(file_mols-shl_mols);
  result[3]=(double)(file_mols-shl_mols);
}

void frame_anal(double *result)
{

  double shl_dr_square = 0.0;
  double shl_dr_fourth = 0.0;
  double non_dr_square = 0.0;
  double non_dr_fourth = 0.0;

  for(int i=0;i<file_mols;i++){

    dx[i]+=(R[i][0]-S[i][0])-box[0]*floor((R[i][0]-S[i][0])/box[0]+0.5);
    dy[i]+=(R[i][1]-S[i][1])-box[1]*floor((R[i][1]-S[i][1])/box[1]+0.5);
    dz[i]+=(R[i][2]-S[i][2])-box[2]*floor((R[i][2]-S[i][2])/box[2]+0.5);
    S[i][0]=R[i][0];
    S[i][1]=R[i][1];
    S[i][2]=R[i][2];

    double dr=dx[i]*dx[i]+dy[i]*dy[i]+dz[i]*dz[i];

    if(flag[i]==1){
      shl_dr_square+=dr;
      shl_dr_fourth+=dr*dr;
    }
    else{
      non_dr_square+=dr;
      non_dr_fourth+=dr*dr;
    }

  }

  result[0]=shl_dr_square;
  result[1]=shl_dr_fourth;
  result[2]=non_dr_square;
  result[3]=non_dr_fourth;

}














