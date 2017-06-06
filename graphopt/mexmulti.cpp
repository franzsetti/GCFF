#include"multi.cpp"
#include "mex.h"
#include<cstdio>
void mexFunction (int nlhs, mxArray *out[],
		  int nrhs, const mxArray *in[]) {
  //        clog<<"Code begins"<<endl;
  // Check to see if we have the correct number of input and output
  // arguments.
  if (nrhs != 8)
    mexErrMsgTxt("Incorrect number of input arguments.\n Correct form is (cost, overlapping neigbourhood, pairwise neighbourhood, pairwise costs, interior labels, weight of exterior points, cost of models,threshold).");
  if (nlhs != 2)
    mexErrMsgTxt("Incorrect number of output arguments should be of the form [Binary matrix,interior labels].");
  printf("check!\n");fflush(0);
  int points = mxGetM(in[0]);
  int hyp   = mxGetN(in[0]);
  int p2 = mxGetM(in[1]);
  int maxn   = mxGetN(in[1]);
  int p2p = mxGetM(in[2]);
  int maxnp   = mxGetN(in[2]);
  int cp2p = mxGetM(in[2]);
  int cmaxnp   = mxGetN(in[2]);
  double lambda=mxGetPr(in[5])[0];
  double* hweight=mxGetPr(in[6]);
  double thresh=mxGetPr(in[7])[0];
  double* pweight=mxGetPr(in[3]);
  if((lambda>1)||(lambda<0)){
    //clog<<"Lambda: "<<lambda<<endl;
    mexErrMsgTxt("Lambda must be in [0,1]");
  }
  ////clog<<"hyp "<<hyp<<endl;
  ////clog<<"points "<<points<<endl;
  ////clog<<"maxn "<<maxn<<endl;
  if (points!=p2){
    //cout<<"points "<<points<<endl;
    //cout<<"p2 "<<p2<<endl;
    mexErrMsgTxt("Matrixs 1 and 2 must be of same length.");
  }
  if (points!=p2p){
    //cout<<"points "<<points<<endl;
    //cout<<"p2 "<<p2<<endl;
    mexErrMsgTxt("Matrixs 1 and 3 must be of same length.");
  }
  if (cp2p!=p2p){
    mexErrMsgTxt("Matrixs 3 and 4 must be the same size.");
  }
  if (cmaxnp!=cmaxnp){
    mexErrMsgTxt("Matrixs 3 and 4 must be the same size.");
  }

  if(maxn>=points)
    mexErrMsgTxt("Matrix 2 has too many neighbours (no neigh must always be larger than points)");

  if(maxnp>=points)
    mexErrMsgTxt("Matrix 3 has too many neighbours (no neigh must always be larger than points)");

  if (points != static_cast<int> (mxGetM(in[4]))){
    //cout<<"points "<<points<<endl;
    //cout<<"mxGetM(in[2]) "<<mxGetM(in[2])<<endl;
    mexErrMsgTxt("Matrix 1, and vector 4 must be of same length.");
  }
  if (mxGetN(in[4])!=1)
    mexErrMsgTxt("Third element must be a vector.");

  if ((mxGetN(in[5])!=1)||(mxGetM(in[5])!=1))
    mexErrMsgTxt("Sixth element must be a scalar.");

  if (mxGetN(in[6])!=1)
    mexErrMsgTxt("7th element must be a vector.");
  if (static_cast<int>(mxGetM(in[6]))!=hyp)
    mexErrMsgTxt("7th element must be a vector of length hyp.");

  if ((mxGetN(in[7])!=1)||(mxGetM(in[7])!=1)||(thresh<0))
    mexErrMsgTxt("8th element must be a positive scalar.");

  for (int i = 0; i !=hyp; ++i)
    if (hweight[i]<0){
      mexErrMsgTxt("matrix 7 must be elementwise >= 0");
    }
  for (int i = 0; i !=hyp*maxnp; ++i)
    if (pweight[i]<0){
      mexErrMsgTxt("Matrix 3 must be elementwise >= 0");
    }
 
  hypothesis H(points,maxn,hyp,lambda,hweight,maxnp);

  H.un= mxGetPr(in[0]);
  for (int i = 0; i !=points*hyp; ++i)
    if(H.un[i]<0){
      mexErrMsgTxt("matrix 1 must be elementwise >= 0");
    }
    else if(H.un[i]>thresh)
      H.un[i]=thresh;
 
  for (int i = 0; i !=points;++i){
    H.label[i]=mxGetPr(in[4])[i];
    if (H.label[i]>=hyp){
      mexErrMsgTxt("All labels must be less than number of hypothesis");
    }
    if (H.label[i]<0)
      mexErrMsgTxt("All labels must be greater than or equal to number 0");
  }

  for (int i = 0; i !=points;++i)
    for (int j  = 0; j !=maxn; ++j){
      H.neigh_m[i*maxn+j]=mxGetPr(in[1])[i*maxn+j]-1;
      if (H.neigh_m[i*maxn+j]>=points)
	mexErrMsgTxt("overlap neighbours (matrix 2) must be less  or equal to than number of points");
      if (H.neigh_m[i*maxn+j]<-1)
	mexErrMsgTxt("overlap neighbours (matrix 2) must be greater than 0");
    }

  for (int i = 0; i !=points;++i)
    for (int j  = 0; j !=maxnp; ++j){
      H.neigh_p[i*maxnp+j]=mxGetPr(in[2])[i*maxnp+j]-1;
      if (H.neigh_p[i*maxnp+j]>=points)
	mexErrMsgTxt("pairwise neighbours (matrix 3) must be less  or equal to than number of points");
      if (H.neigh_p[i*maxnp+j]<-1){
	char a[80];
	sprintf(a,"points: %d H.neigh_p[%d]=%d i:%d maxnp:%d j:%d",points,i*maxnp+j,H.neigh_p[i*maxnp+j],i,maxnp,j);
	 mexErrMsgTxt(a);
	mexErrMsgTxt("pairwise neighbours (matrix 3) must be greater than or equal to 0");
      }
    }

  H._ncost=pweight;
 
  H.solve();
  
  out[0]=mxCreateDoubleMatrix(points,hyp,mxREAL);   
  H.annotate((double*)mxGetPr(out[0]));
  for (int i = 0; i !=hyp*points; ++i)
    if (H.un[i]==thresh)
      mxGetPr(out[0])[i]=0;

    
  out[1]=mxCreateDoubleMatrix(points,1,mxREAL);
  for (int i = 0; i !=points; ++i)
    mxGetPr(out[1])[i]=H.label[i];
 

  return ;
}
 
