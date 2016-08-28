/**************************************************************************************************************/
/*** TMSD, TFQT, TCHI, TCLS: DIM, ATM  ***/

double box[3];  // box dimensions
double (*R)[3]; // current  atomic coordinates
double (*S)[3]; // previous atomic coordinates
double (*T)[3]; // cumulative atomic displacement

#if FORMAT==DIM
  double (*pos)[2][3];
#endif

void frame_setm()
{
  R=(double*)malloc(3*file_mols*sizeof(double));
  S=(double*)malloc(3*file_mols*sizeof(double));
  T=(double*)malloc(3*file_mols*sizeof(double));
  #if FORMAT==DIM
    pos=(double*)malloc(6*file_mols*sizeof(double));
  #endif
}

void frame_free()
{
  free(R);
  free(S);
  free(T);
  #if FORMAT==DIM
    free(pos);
  #endif
}

void frame_read()
{
  #if FORMAT==DIM

  for(int i=0;i<file_mols;i++){
    for(int j=0;j<6;j++)
      fscanf(qfile,"%lf",&pos[i][j/3][j%3]);
  }
  fscanf(qfile,"%lf %lf %lf",&box[0],&box[1],&box[2]);

  double sep[3];

  for(int i=0;i<file_mols;i++){

    sep[0]=pos[i][0][0]-pos[i][1][0]-box[0]*floor((pos[i][0][0]-pos[i][1][0])/box[0]+0.5);
    sep[1]=pos[i][0][1]-pos[i][1][1]-box[1]*floor((pos[i][0][1]-pos[i][1][1])/box[1]+0.5);     
    sep[2]=pos[i][0][2]-pos[i][1][2]-box[2]*floor((pos[i][0][2]-pos[i][1][2])/box[2]+0.5);

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

  #elif FORMAT==ATM

  for(int i=0;i<3*file_mols;i++)
    fscanf(qfile,"%lf",&R[i/3][i%3]);
  fscanf(qfile,"%lf %lf %lf",&box[0],&box[1],&box[2]);

  #endif

}

void frame_norm(double *result)
{

  for(int i=0;i<file_mols;i++){
    T[i][0]=0.0;
    T[i][1]=0.0;
    T[i][2]=0.0;
    S[i][0]=R[i][0];
    S[i][1]=R[i][1];
    S[i][2]=R[i][2];
  }

  result[NQ%2]=(double)file_mols;
  result[NQ%2+1]=(double)file_mols;
  #if ANALYSIS>=TFQT
  result[0]=(double)(3*file_mols);
  #endif
  #if ANALYSIS>=TCHI
  result[3]=1.0;
  result[4]=1.0;
  #endif
  #if ANALYSIS>=TCLS
  result[5]=1.0;
  result[6]=1.0;
  #endif

}

void frame_anal(double *result)
{

  double dr_square=0.0;
  double dr_fourth=0.0;
  #if ANALYSIS>=TFQT
  double dr_fqtfqt=0.0;
  #endif
  #if ANALYSIS>=TCHI
  double dr_chione=0.0;
  double dr_chitwo=0.0;
  double *ds;
  ds=malloc(file_mols*sizeof(double));
  #endif
  #if ANALYSIS>=TCLS
  int current;
  int cluster_ind;
  int limit_ind;
  int mobil_ind[MOB];
  int molecule_ind[MOB];
  int cluster[MOB][MOB];
  double mobil_val[MOB];
  double limit_val=box[0]*box[0];
  #endif

  for(int i=0;i<file_mols;i++){

    T[i][0]+=(R[i][0]-S[i][0])-box[0]*floor((R[i][0]-S[i][0])/box[0]+0.5);
    T[i][1]+=(R[i][1]-S[i][1])-box[1]*floor((R[i][1]-S[i][1])/box[1]+0.5);
    T[i][2]+=(R[i][2]-S[i][2])-box[2]*floor((R[i][2]-S[i][2])/box[2]+0.5);
    S[i][0]=R[i][0];
    S[i][1]=R[i][1];
    S[i][2]=R[i][2];

    #if ANALYSIS>=TFQT
    dr_fqtfqt+=cos(WAV*T[i][0])+cos(WAV*T[i][1])+cos(WAV*T[i][2]);
    #endif

    double dr=T[i][0]*T[i][0]+T[i][1]*T[i][1]+T[i][2]*T[i][2];

    dr_square+=dr;
    dr_fourth+=dr*dr;

    #if ANALYSIS>=TCHI
    ds[i]=sqrt(dr);
    dr_chione+=ds[i];
    for(int j=0;j<=i;j++)
      dr_chitwo+=ds[i]*ds[j];
    #endif

    #if ANALYSIS>=TCLS
    if(i<MOB){
      mobil_val[i]=dr;
      mobil_ind[i]=i;
      if(dr<limit_val){
        limit_val=dr;
        limit_ind=i;
      }
    }
    else{
      if(dr>limit_val){
        mobil_val[limit_ind]=dr;
        mobil_ind[limit_ind]=i;
        limit_val=dr;
        for(int j=0;j<MOB;j++){
          if(mobil_val[j]<limit_val){
            limit_val=mobil_val[j];
            limit_ind=j;
          }
        }
      }
    }
    #endif

  }

  #if ANALYSIS>=TCLS
  for(int i=0;i<MOB;i++){
    molecule_ind[i]=-1;  //the clusters have no number of molecules yet
    for(int j=0;j<MOB;j++)
      cluster[i][j]=-1;  //the clusters have no molecules yet
  }
  cluster_ind=0;  //we are on the first cluster
  for(int i=0;i<MOB;i++){  //go through each mobile molecule
    if(mobil_ind[i]>=0){  //if it hasn't been added to a cluster yet:
      molecule_ind[cluster_ind]=0;  //it is the first molecule in this cluster
      cluster[cluster_ind][molecule_ind[cluster_ind]]=mobil_ind[i];  //add the mobile molecule
      mobil_ind[i]=-1;  //mark that it has been added to a cluster
      molecule_ind[cluster_ind]++;  //the next molecule added to this cluster will be this index
      current=0;  //we're going to start by searching the first molecule for mobile neighbors
      while(current<molecule_ind[cluster_ind]){  //until we've searched everything we've added for neighbors
        for(int j=0;j<MOB;j++){  //go through each mobile molecule
          if(mobil_ind[j]>=0){  //if it hasn't been added to a cluster yet
            double dx=S[cluster[cluster_ind][current]][0]-S[mobil_ind[j]][0];
            dx=dx-box[0]*floor(dx/box[0]+0.5);
            double dy=S[cluster[cluster_ind][current]][1]-S[mobil_ind[j]][1];
            dy=dy-box[1]*floor(dy/box[1]+0.5);
            double dz=S[cluster[cluster_ind][current]][2]-S[mobil_ind[j]][2];
            dz=dz-box[2]*floor(dz/box[2]+0.5);
            double dr=dx*dx+dy*dy+dz*dz;  //calculate it's (original) distance from the current search molecule
            if(dr<CUT*CUT){  //if it is a neighbor
              cluster[cluster_ind][molecule_ind[cluster_ind]]=mobil_ind[j];  //add it to the cluster
              mobil_ind[j]=-1;  //mark that it has been added to a cluster
              molecule_ind[cluster_ind]++;  //mark that we have added a molecule to the cluster
            }
          }
        }  //we've now searched all the unadded mobile molecules to see if they are a n.n. of the current
        current++;  //move on to the next
      }  //we've now searched the whole list
      cluster_ind++;  //move on to the next cluster
    }
  }  //we've now started as many clusters as possible
  #endif

  result[NQ%2]=dr_square;
  result[NQ%2+1]=dr_fourth;
  #if ANALYSIS>=TFQT
  result[0]=dr_fqtfqt;
  #endif
  #if ANALYSIS>=TCHI
  free(ds);
  result[3]=dr_chione;
  result[4]=dr_chitwo;
  #endif
  #if ANALYSIS>=TCLS
  for(int i=0;i<cluster_ind;i++){
    result[5]+=molecule_ind[i];
    result[6]+=molecule_ind[i]*molecule_ind[i];
  }
  #endif

}









