/**************************************************************************************************************/
/*** CORR: DAT  ***/

  double A;  // current  value
  double B;  // previous value

  void frame_setm(){}

  void frame_free(){}

  void frame_read()
  {
    fscanf(qfile,"%lf",&A);
  }

  void frame_norm(double *result)
  {
    // depending on normalization condition:
    result[0]=1.0;
    //result[0]=A*A;
    B=A;
  }

  void frame_anal(double *result)
  {
    result[0]=(A-OFF)*(B-OFF);
  }





















