#include"class.cpp"
#include "mex.h"
#include<cstdio>
void mexFunction (int nlhs, mxArray *out[],
		  int nrhs, const mxArray *in[]) {
// Check to see if we have the correct number of input and output
  // arguments.
  if (nrhs != 6)
    mexErrMsgTxt("Incorrect number of input arguments.\n Correct form is (cost, neigbourhood, interior labels, weight of exterior points, cost of models,threshold).");
  if (nlhs != 2)
    mexErrMsgTxt("Incorrect number of output arguments should be of the form [Binary matrix,interior labels].");


  int points =  mxGetM(in[0]);
  int hyp   =mxGetN(in[0]);
  int p2 = mxGetN(in[1]);
  int maxn   = mxGetM(in[1]);
  double lambda=mxGetPr(in[3])[0];
  double* hweight=mxGetPr(in[4]);
  double thresh=mxGetPr(in[5])[0];;
mexPrintf("thresh %f\n",thresh);
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
  if(maxn>=points)
    mexErrMsgTxt("Matrix 2 has too many neighbours (no neigh must always be less than points)");

  if (points != static_cast<int> (mxGetM(in[2]))){
    //cout<<"points "<<points<<endl;
    //cout<<"mxGetM(in[2]) "<<mxGetM(in[2])<<endl;
    mexErrMsgTxt("Matrix 1, and vector 3 must be of same length.");
  }
  if (mxGetN(in[2])!=1)
    mexErrMsgTxt("Third element must be a vector.");

  if (mxGetN(in[3])!=1)
    mexErrMsgTxt("Fourth element must be a scalar.");

  if (mxGetM(in[3])!=1)
    mexErrMsgTxt("Fourth element must be a scalar.");

  if (mxGetN(in[4])!=1)
    mexErrMsgTxt("Fifth element must be a vector.");
  if (static_cast<int>(mxGetM(in[4]))!=hyp)
    mexErrMsgTxt("Fifth element must be a vector of length hyp.");

  if (mxGetN(in[5])!=1)
    mexErrMsgTxt("Sixth element must be a scalar.");

  if (mxGetM(in[5])!=1)
    mexErrMsgTxt("Sixth element must be a scalar.");

  for (int i = 0; i !=hyp; ++i)
    if (hweight[i]<0)
      mexErrMsgTxt("hweight[i] must be >= 0");

  hypothesis H(points,maxn,hyp,lambda,hweight,0);

  vector<double>un(points*hyp);
 double* tmp= mxGetPr(in[0]);
 H.un=&un[0];
 for(unsigned int i=0;i!=points*hyp;++i)
      if(tmp<0){
      mexErrMsgTxt("hweight[i] must be >= 0");
    }
    else 
      un[i]=min(thresh,tmp[i]);
 
  for (int i = 0; i !=points;++i){
    H.label[i]=mxGetPr(in[2])[i];
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
	mexErrMsgTxt("neighbours must be less  or equal to than number of points");
      if (H.neigh_m[i*maxn+j]<-1)
	mexErrMsgTxt("neighbours must be greater than 0");
      if(H.neigh_m[i*maxn+j]==i){
	mexPrintf("Matrix 2 %d has a self neighbour in position %d\n",i+1,j+1);
	mexErrMsgTxt("No self neighbours");
      }
    }
 
  
  H.solve();
  
  out[0]=mxCreateDoubleMatrix(points,hyp,mxREAL);   
  H.annotate((double*)mxGetPr(out[0]),thresh);

    
  out[1]=mxCreateDoubleMatrix(points,1,mxREAL);
  for (int i = 0; i !=points; ++i)
    mxGetPr(out[1])[i]=H.label[i];
 

  return ;
}
 
