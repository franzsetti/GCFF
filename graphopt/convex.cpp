#include <vector>
#include "mex.h"
#include "math.h"
#include <iostream>
using namespace std;
#include<cstdio>
#define xstr(s) str(s)
#define str(s) #s

#define DEBUG __
#ifdef DEBUG
  #define mexassert(x) \
    {if(! (x))\
	mexErrMsgTxt("ERROR!! Mexassert: " #x  "failed\n  on line" xstr(__LINE__) "in file " __FILE__ "\n");}
#else
#define mexassert(x)
#endif
#include "graph.h"
#include "graph.cpp"
#include "maxflow.cpp"

using namespace std;
template class Graph<float,float,float>;

class conv{
public:
  unsigned int max_neigh_p;
  unsigned int nodes, hyp;
  double *un;

  unsigned int maxnodes, maxedges;

  vector<int> neigh_p;
  double* _ncost;
  Graph<float,float,float>* g;

  inline  void add_tweights(unsigned int i,float w1,float w2, char id){
    //    cout<<"id "<<id <<endl;

    if ((w1==0)&&(w2==0))
      return;
    mexassert(i>=0);
    //    if (i>=maxnodes)
    //  cout<<"i"<<i<<endl;

    mexassert(i<maxnodes);
     g-> add_tweights(i,w1,w2);
  }
  inline void add_edge(unsigned int i,unsigned int j,float w1,float w2,char id){
    // cout<<"edge id "<<i<<' '<<j <<endl;
    // cout<<"w1 "<<w1<<endl;
    // cout<<"w2 "<<w2<<endl;


    mexassert(i!=j);
    mexassert(i>=0);
    mexassert(j>=0);
    mexassert(i<maxnodes);
    mexassert(j<maxnodes);
    mexassert(w1>=0.0);
    mexassert(w2>=0.0);
    if ((w1==0)&&(w2==0))
      return;
    g->add_edge(i,j,w1,w2);
  }
  inline float unary(unsigned int node,unsigned int lab){
    mexassert(node<nodes);
    mexassert(lab<hyp);
    //    cout<<"un["<<node<<"+"<<lab<<"*nodes] "<<un[node+lab*nodes]<<endl;

    return un[node+lab*nodes];
  }
  int pair_neighbour(unsigned int n,unsigned int n2){
    //    cout<<"n"<<n<<endl;
    // cout<<"n2"<<n2<<endl;
    //    cout<<"max_neigh_p"<<max_neigh_p<<endl;

    mexassert(n<nodes);
    mexassert(n2<max_neigh_p);
    //return neigh[n*max_neigh+n2];
    //    cout<<"neigh_p[n+nodes*n2]"<<neigh_p[n+nodes*n2]<<endl;

    return neigh_p[n+nodes*n2];
  }
  float ncost(unsigned int n,unsigned int n2){
    mexassert(_ncost);
    mexassert(n<nodes);
    mexassert(n2<max_neigh_p);
    mexassert(!isnan(_ncost[n+nodes*n2]));    
    mexassert(_ncost[n+nodes*n2]>=0);
    mexassert(_ncost[n+nodes*n2]!=INFINITY);
    return _ncost[n+nodes*n2];
  }
  conv(int nodes,int hyp,int maxn_p){
    max_neigh_p=maxn_p;
    this->nodes=nodes;
    this->hyp=hyp;
    //    cout<<"hyp"<<hyp<<endl;

    neigh_p.resize(nodes*max_neigh_p,-1);
    //    max_neigh_p=0;
    maxnodes=(hyp-1)*nodes;
    //    cout<<"maxnodes "<<maxnodes <<endl;
    maxedges=((hyp-1)*nodes*max_neigh_p+1)/2+maxnodes;
    // cout<<"maxedges "<<maxedges <<endl;

    g=new Graph<float,float,float> (maxnodes,maxedges);
    g->add_node(maxnodes);
    cout<<"allocated memory"<<endl;

  }

  ~conv(){
    delete g;
  }
  void construct(){
    for (unsigned int i = 0; i !=nodes; ++i){
   
      cout<<"i"<<i<<endl;

      add_tweights(i*(hyp-1),0,unary(i,0),'a');
      add_tweights(i*(hyp-1)+hyp-2,unary(i,hyp-1),0,'a');
      for (unsigned int j = 1; j !=hyp-1; ++j)
	add_edge(i*(hyp-1)+j-1,i*(hyp-1)+j,INFINITY,unary(i,j),'a');
          cout<<"pair"<<endl; 
   
       for (unsigned int j = 0; (j !=max_neigh_p)&&(pair_neighbour(i,j)!=-1);++j)
	 if (((unsigned int) pair_neighbour(i,j))<i)
       	  for (unsigned int k = 0; k !=hyp-1; ++k){
	      	    add_edge(i*(hyp-1)+k,pair_neighbour(i,j)*(hyp-1)+k,ncost(i,j),ncost(i,j),'a');//Must be symetric
	  }

      
    }
  }
  void solve(){
        cout<<"calling construct"<<endl;
    construct();
    cout<<"constructed"<<endl;
    g->maxflow();
    cout<<"flow complete"<<endl;
  }

  void annotate (double* out){
    mexassert(out);
    for (unsigned int i = 0; i !=nodes; ++i){
      out[i]=0;
      for (unsigned int j = 0; j !=hyp-1; ++j){
	//	cout<<"i*(hyp-1)+j"<<i*(hyp-1)+j<<endl;
	bool ass=g->what_segment(i*(hyp-1)+j,Graph<float,float,float>::SINK)==Graph<float,float,float>::SINK;
	//	cout<<"ass "<<ass <<endl;
	if(ass)
	  out[i]=j+1;
	 else
	   break;
      }
    }
    return ;
  }

};

void mexFunction (int nlhs, mxArray *out[],
		  int nrhs, const mxArray *in[]) {
 if (nrhs != 4)
    mexErrMsgTxt("Incorrect number of input arguments.\n Correct form is (1:unary cost, 2:neigbourhood, 3:neighbour costs, 4:outlier rejection).");
 if (nlhs != 1)
    mexErrMsgTxt("Incorrect number of output arguments should be of the form [labels]");
 unsigned int points = mxGetM(in[0]);
 unsigned int hyp   = mxGetN(in[0]);
 unsigned int p2 = mxGetM(in[1]);
 unsigned int maxn   = mxGetN(in[1]);
 double thresh=mxGetPr(in[3])[0];
 double* pweight=mxGetPr(in[2]);
 ////clog<<"hyp "<<hyp<<endl;
 ////clog<<"points "<<points<<endl;
 ////clog<<"maxn "<<maxn<<endl;
  if (points!=p2){
    //cout<<"points "<<points<<endl;
    //cout<<"p2 "<<p2<<endl;
    mexErrMsgTxt("Matrixs 1 and 2 must be of same length.");
  }
  if (hyp<=2){
    //cout<<"points "<<points<<endl;
    //cout<<"p2 "<<p2<<endl;
    mexErrMsgTxt("Need at least 3 hypothesis");
  }

  if ((mxGetM(in[2])!=p2)||mxGetN(in[2])!=maxn){
    mexErrMsgTxt("second and third arguements must be matrices of same dimension");
  }
  
  conv H(points,hyp,maxn);
  H.un= mxGetPr(in[0]);
  for (unsigned int i = 0; i !=points*hyp; ++i)
    if(H.un[i]<0){
      //clog<<"i: "<<i<<"H.un[i]: "<<H.un[i]<<endl;
      mexErrMsgTxt("hweight[i] must be >= 0");
    }
  H._ncost=pweight;
  
  //H.neigh= mxGetPr(in[1]); 
  for (unsigned int i = 0; i !=points;++i)
    for (unsigned int j  = 0; j !=maxn; ++j){
      H.neigh_p[i*maxn+j]=mxGetPr(in[1])[i*maxn+j]-1;
      if ((H.neigh_p[i*maxn+j]>=0)&&(((unsigned int) H.neigh_p[i*maxn+j]) >=points))
	mexErrMsgTxt("neighbours must be less  or equal to than number of points");
      if (H.neigh_p[i*maxn+j]<-1)
	mexErrMsgTxt("neighbours must be greater than 0");
      // //clog<<H.neigh[i*maxn+j]<<endl;
    } 
  H.solve();
  
  out[0]=mxCreateDoubleMatrix(points,1,mxREAL);   
  H.annotate((double*)mxGetPr(out[0]));
  // for (unsigned int i = 0; i !=points; ++i)
  //   if (H.un[]==thresh)
  //     mxGetPr(out[0])[i]=0;

    
  return ;
}
 
