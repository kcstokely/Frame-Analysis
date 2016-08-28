/**************************************************************************************************************/
/*** VISC: G96 ***/

double (*R)[3]; // atomic positions
double (*V)[3]; // atomic velocities

void frame_setm()
{
  R=(double*)malloc(3*file_mols*sizeof(double));
  V=(double*)malloc(3*file_mols*sizeof(double));
}

void frame_free()
{
  free(R);
  free(V);
}

double box[3];

void frame_read()
{
  fseek(qfile,9+31+4+12,SEEK_CUR);
  for(int i=0;i<3*file_mols;i++)
    fscanf(qfile,"%lf",&R[i/3][i%3]);
  fseek(qfile,4+12,SEEK_CUR);
  for(int i=0;i<3*file_mols;i++)
    fscanf(qfile,"%lf",&V[i/3][i%3]);
  fseek(qfile,4+4,SEEK_CUR);
  fscanf(qfile,"%lf %lf %lf",&box[0],&box[1],&box[2]);
  fseek(qfile,4,SEEK_CUR);

  frame_calc();

}

const int type[3002] = {[0 ... 1] = 0, [2 ... 2401] = 1, [2402 ... 3001] = 2};

const double sigma[3][3]   = {{0.00, 1.00, 0.80},
                              {1.00, 1.00, 0.80},
                              {0.80, 0.80, 0.88}};

const double epsilon[3][3] = {{0.00, 1.00, 1.50},
                              {1.00, 1.00, 1.50},
                              {1.50, 1.50, 0.50}};

double sxy[2];
double sxz[2];
double syz[2];

void frame_calc(){
  sxy[1]=0.0;
  sxz[1]=0.0;
  syz[1]=0.0;
  for(int i=0;i<file_mols;i++){
    sxy[1]+=V[i][0]*V[i][1];
    sxy[1]+=V[i][0]*V[i][2];
    sxy[1]+=V[i][1]*V[i][2];
    for(int j=i+1;j<file_mols;j++){
      double dx=(R[i][0]-R[j][0])-box[0]*floor((R[i][0]-R[j][0])/box[0]+0.5);
      double dy=(R[i][1]-R[j][1])-box[1]*floor((R[i][1]-R[j][1])/box[1]+0.5);
      double dz=(R[i][2]-R[j][2])-box[2]*floor((R[i][2]-R[j][2])/box[2]+0.5);
      double dr=sqrt(dx*dx+dy*dy+dz*dz);
      double ff=12.0*pow(sigma[type[i]][type[j]],12)/pow(dr,13);
      ff-=6.0*pow(sigma[type[i]][type[j]],6)/pow(dr,7);
      ff*=4.0*epsilon[type[i]][type[j]];
      sxy[1]+=ff*dx*dy/dr;
      sxy[1]+=ff*dx*dz/dr;
      sxy[1]+=ff*dy*dz/dr;
    }
  }
}

void frame_norm(double *result){
  result[0]=1.0;
  result[1]=1.0;
  result[2]=1.0;
  sxy[0]=sxy[1];
  sxz[0]=sxz[1];
  syz[0]=syz[1];
}

void frame_anal(double *result){
  result[0]=sxy[0]*sxy[1];
  result[1]=sxz[0]*sxz[1];
  result[2]=syz[0]*syz[1];
}


































