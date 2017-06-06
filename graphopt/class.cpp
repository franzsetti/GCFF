#ifdef _MSC_VER
    #define isnan(x) _isnan(x)
    
    double _stupid_hack=0;
    double INFINITY=1.0/_stupid_hack;
#endif
#include <string>
//include <iostream>
#include <cmath>
#include <vector>
#include "graph.h"
#include "graph.cpp"
#include "maxflow.cpp"
#include <mex.h>
#define xstr(s) str(s)
#define str(s) #s
#define DEBUG ___
#ifdef DEBUG
  #define mexassert(x) \
    {if(! (x))\
	mexErrMsgTxt("ERROR!! Mexassert: " #x  "failed\n  on line" xstr(__LINE__) "in file " __FILE__ "\n");}
#else
#define mexassert(x)
#endif

using namespace std;

template class Graph<double,double,double>; 

class hypothesis{
public:
  int max_neigh_p;
  int max_neigh_m;
  int nodes, hyp;
  vector<int>  label;
  double *un;

  double lambda;
  double *hypweight;
  int maxnodes;

  vector<int> neigh_p;
  vector<int> neigh_m;
  double* _ncost;

  inline  void add_tweights(int i,double w1,double w2, char id){
    if ((w1==0)&&(w2==0))
      return;
    mexassert(i>=0);
    mexassert(i<maxnodes);
    g-> add_tweights(i,w1,w2);
  }
  inline void add_edge(int i,int j,double w1,double w2,char id){
    //mexPrintf("i %d j %d, id %c\n",i,j,id);
    mexassert(i!=j);
    mexassert(i>=0);
    mexassert(j>=0);
    if (i>=maxnodes ||j >= maxnodes)
        mexPrintf("i %d j %d\n",i,j);
    mexassert(i<maxnodes);
    mexassert(j<maxnodes);
    mexassert(w1>=0.0);
    mexassert(w2>=0.0);
    if ((w1==0)&&(w2==0))
      return;
    g->add_edge(i,j,w1,w2);
  }
  inline double& unary(int node,int lab){
     mexassert(node>=0);
     mexassert(lab>=0);
    mexassert(node<nodes);
    mexassert(lab<hyp);
    return un[node+lab*nodes];
  }

  inline double hweight(int lab){
    mexassert(hypweight);
    mexassert(lab<hyp);
    return hypweight[lab];
  }
  inline int overlap_neighbour(int n,int n2){
    mexassert(n>=0);
    mexassert(n<nodes);
    mexassert(n2<max_neigh_m);
    mexassert(n2>=0);
    mexassert (neigh_m[n*max_neigh_m+n2]<nodes);
    mexassert (neigh_m[n*max_neigh_m+n2]>=-1);
    return neigh_m[n*max_neigh_m+n2];
  }
  int pair_neighbour(int n,int n2){
    mexassert(n<nodes);
    mexassert(n2<max_neigh_p);
    return neigh_p[n*max_neigh_p+n2];
    //return neigh_p[n+nodes*n2];
  }
  double ncost(int n,int n2){
    mexassert(_ncost);
    mexassert(n<nodes);
    mexassert(n2<max_neigh_p);
    mexassert(!isnan(_ncost[n*max_neigh_p+n2]));    
    mexassert(_ncost[n*max_neigh_p+n2]>=0);
    mexassert(_ncost[n*max_neigh_p+n2]!=INFINITY);
    return _ncost[n*max_neigh_p+n2];
  }
  Graph<double,double,double>* g;
  hypothesis(int nodes,int maxn, int hyp,double l,double* w1,int maxn_p){
    max_neigh_p=maxn_p;
    hypweight=w1;
    lambda=l;
    max_neigh_m=maxn;
    this->nodes=nodes;
    int edges=maxn*nodes;
    this->hyp=hyp;
    label.resize(nodes,0);
    neigh_m.resize(nodes*max_neigh_m,-1);
    neigh_p.resize(nodes*max_neigh_p,-1);
    //    max_neigh_p=0;
    maxnodes=2*nodes/*nothing which is alpha*/ + edges/* non interior points*/ +hyp/*MDL prior*/;
    mexPrintf("nodes %d, edges %d hyp %d\n", nodes, edges, hyp);
    g=new Graph<double,double,double> (maxnodes,edges*2+nodes*3+max_neigh_p*nodes);
  }
  ~hypothesis(){
    delete g;
  }

  void construct_multi(const int alpha){
    g->add_node(nodes);
    for (int i = 0; i !=nodes; ++i)	
      if(label[i]!=alpha){
	bool nalpha=false;
	//check if overlap_neighbours contain alpha
	for (int j = 0; (j !=max_neigh_m)&&(overlap_neighbour(i,j)!=-1); ++j)
	  nalpha|=(label[overlap_neighbour(i,j)]==alpha);
	
	if(!nalpha) 
	  // if(nalpha) then boundary set must always contain alpha
	  //and transition to I=alpha is free
	  //and remaining `unary' cost comes in section 2
	  //else
	  {
	    add_tweights(nodes+i,unary(i,alpha)*lambda,0,'d');
	    add_edge(i,nodes+i,0,unary(i,alpha)*lambda,'a');
	    for (int j = 0; (j !=max_neigh_m)&&(overlap_neighbour(i,j)!=-1); ++j)
	      add_edge(overlap_neighbour(i,j),nodes+i,0,lambda*unary(i,alpha),'b');
	  }
      }
    //Half way --- we have added all costs associated with changing node n to alpha 
    // we now need to add the costs to overlap_neighbours of keeping a node fixed

    vector<bool>nclass(hyp);    
    for (int i = 0; i !=nodes; ++i){
      //      //clog<<"i "<<i<<endl;
      for (int j = 0; j !=hyp; ++j)
	nclass[j]=false;
      for (int j = 0; (j !=max_neigh_m)&&(overlap_neighbour(i,j)!=-1); ++j)	
	nclass[label[overlap_neighbour(i,j)]]=true; //check overlap_neighbours
      nclass[label[i]]=true;//check self

      for (int j = 0; j !=hyp; ++j){
	if(nclass[j]&&(j!=alpha)){
	  int temp=g->add_node(1);

	  add_tweights(temp,0,unary(i,j)*lambda,'e');
	  for (int k = 0; (k !=max_neigh_m)&&(overlap_neighbour(i,k)!=-1); ++k)
	    if(label[overlap_neighbour(i,k)]==j)
	      add_edge(overlap_neighbour(i,k),temp,unary(i,j)*lambda,0,'c');
	  if(label[i]==j){
	    add_edge(i,temp,unary(i,j)*lambda,0,'d');
	    add_tweights(i,unary(i,alpha)*(1.0-lambda),unary(i,j)*(1.0-lambda),'b');
	  }
	}
      }
    }
  }
  void construct_mdl_boykov(const int alpha){
    vector<bool>nclass(hyp);    
    for (int j = 0; j !=hyp; ++j)
      nclass[j]=false;
    
    for (int i = 0; i !=nodes; ++i){
      nclass[label[i]]=true;
    }

    
    for (int i = 0; i !=hyp; ++i){
      if ((i!=alpha)&&(nclass[i])&&hweight(i)){
	
	int temp= g->add_node();
	add_tweights(temp,0,hweight(i),'z');
	for (int j = 0; j !=nodes; ++j)
	  if (label[j]==i){
	    add_edge(j,temp,hweight(i),0,'l');
	  }
      }
    }
  }
  void construct_mdl_standard(const int alpha){
    //Not needed see Yuri Boykov's work on MDL priors
    //Left for debugging
    construct_mdl_boykov(alpha);
    vector<bool>nclass(hyp);    
    for (int j = 0; j !=hyp; ++j)
      nclass[j]=false;
    
    for (int i = 0; i !=nodes; ++i){
      nclass[label[i]]=true;
    }
    
    if(!nclass[alpha]&&hweight(alpha)){
      int temp=g->add_node();
      add_tweights(temp,hweight(alpha),0,'a');
      for (int i = 0; i !=nodes; ++i)
	add_edge(i,temp,0,hweight(alpha),'m');
    }  

  }

  void construct_pairwise(const int alpha){//Must be called after construct_multi
    //mexPrintf("Contructing graph for alpha=%d\n",alpha);
    for ( int i = 0; i !=nodes; ++i){	
      //mexPrintf("Node %d labelled %d \n",i,label[i]);
      if(label[i]!=alpha){
	//mexPrintf("not alpha\n");	
	mexassert(label[i]<hyp);
	double cost;
	cost=(1-lambda)*(unary(i,label[i])-unary(i,alpha));
	for ( int j = 0; (j !=max_neigh_p)&&(pair_neighbour(i,j)!=-1); ++j){
	  //mexPrintf(" %dth neighbour is  %d \n",j,pair_neighbour(i,j));	
	  if(label[pair_neighbour(i,j)]==alpha){
	    //mexPrintf("neighbour %d is taking label alpha\n",j);	
	    cost+=ncost(i,j);
	  }
	  else// if (i>pair_neighbour(i,j)) 
	    {
	    if (label[pair_neighbour(i,j)]==label[i]){
	      add_edge(i,pair_neighbour(i,j),ncost(i,j),ncost(i,j),'a');//Tautologically correct
	    }
	    else
	      {
		cost+=ncost(i,j);
		add_edge(i,pair_neighbour(i,j),0,ncost(i,j),'a');
		/*
		  X|1 0
		  ----
		  1|1 1
		  0|1 0
		  = x_1 + x_2 -x_1 x_2 = x_1 + (1-x_1)x_2
		*/
	      }
	  }
	}
	if(label[i]!=alpha){
	  add_tweights(i,max(-cost,0.0),max(cost,0.0),'z'); // Add unaries backward.
	  //mexPrintf("Added unary cost %f\n",cost);	
	}
      }
    }
  }
  
  void expand(const int alpha) {
    //    disp_labels();
    //mexPrintf("Expanding on %d:\n",alpha);
    mexassert(alpha<hyp);
    g->reset();
    g->add_node(nodes);
    
    if(max_neigh_m)
      construct_multi(alpha);
    //if(max_neigh_p)
    construct_pairwise(alpha);
    construct_mdl_boykov(alpha);
    //construct_mdl_standard(alpha);
    
    g->maxflow();	  

    for (int i = 0; i !=nodes; ++i){
      if ((label[i]!=alpha)&&
	  (g->what_segment(i,Graph<double,double,double>::SINK)==Graph<double,double,double>::SINK))
	label[i]=alpha;
    }
  }
  void disp_labels(){
    for (int i = 0; i !=nodes; ++i)
      mexPrintf("%d ",label[i]);
    mexPrintf("\n");
  }

  void annotate (double* out,double thresh){
      mexassert(out);
      for (int i = 0; i !=nodes*hyp; ++i)
          out[i]=0;
      
      for (int i = 0; i !=nodes; ++i){
          for (int j = 0; (j !=max_neigh_m&&overlap_neighbour(i,j)!=-1); ++j)
              out[i+nodes*label[overlap_neighbour(i,j)]]=max(out[i+nodes*label[overlap_neighbour(i,j)]],lambda);
          out[i+nodes*label[i]]=1;
      }
      for (int j=0;j!=hyp;++j)
          for (int i = 0; i !=nodes; ++i)
              if(unary(i,j)==thresh)
                  out[i+nodes*j]=0;
      return ;
  }

  double cost() {
    double cost=0;
    vector<double>lweight(hyp);
    for ( int i = 0; i != nodes ; ++i)
      {
	for ( int j = 0; j != hyp ; ++j)
	  lweight[j]=0;
	
	for (int j = 0; (j !=max_neigh_m&&overlap_neighbour(i,j)!=-1); ++j)
	  lweight[label[overlap_neighbour(i,j)]]=lambda;
	lweight[label[i]]=1;
	
	for (int j = 0; j != hyp ; ++j)
	  if(lweight[j]!=0)
	    cost+=lweight[j]*unary(i,j);
      }

    for ( int i = 0; i != nodes ; ++i)
      	for ( int j = 0; (j !=max_neigh_p)&&(pair_neighbour(i,j)!=-1); ++j)
	  if (label[i]!=label[pair_neighbour(i,j)])
	    cost+=ncost(i,j);

    vector<int>lcost(hyp,0);
    for (int i = 0; i !=nodes; ++i){
      mexassert(label[i]<hyp);
      mexassert(label[i]>=0);
      lcost[label[i]]=1;
    }
    for (int i = 0; i !=hyp; ++i)
      cost+=lcost[i]*hweight(i);
    return cost;
  }

  void debug_solve(){
    double c=cost();
   
    double c2=INFINITY;
    vector<int> oldlabel(nodes);
    for (int i = 0; i !=nodes; ++i)
      oldlabel[i]=label[i];

    int i=0;
    while ((c2>c)||(c==INFINITY)){
      c2=c;
      double temp=0;

      for (int j = 0; j !=hyp; ++j){
	printf("i: %d, j: %d, hyp %d\n",i,j,hyp);
	expand(j);
	//      annotate(grid);
	temp=cost();
	//mexassert(temp==new_cost());
	if (temp>c){
	  char a[80];
	  sprintf(a,"Cost increasing j:%d i:%d c: %f temp:%f",j,i,c,temp);
	  mexErrMsgTxt(a);
	  for (int i = 0; i !=nodes; ++i)
	    label[i]=oldlabel[i];	
	  temp=c;
	}
	else {
	  if (temp<c){
	    printf ("Cost decreasing to %f\n",temp);
	    i=0;
	  }
	  for (int i = 0; i !=nodes; ++i)
	    oldlabel[i]=label[i];
	  c=temp;
	}
	++i;
	if (i==hyp){
	  c2=c;
	  return;
	} 
      }
    }
}

void fast_solve(){
    double c=cost();
   
    double c2=INFINITY;
    vector<int> oldlabel(nodes);
    for (int i = 0; i !=nodes; ++i)
      oldlabel[i]=label[i];

    int i=0;
    while ((c2>c)||(c==INFINITY)){
      c2=c;
      double temp=0;

      for (int j = 0; j !=hyp; ++j){
	//printf("i: %d, j: %d, hyp %d\n",i,j,hyp);
	expand(j);
	//      annotate(grid);
	temp=cost();
	//mexassert(temp==new_cost());
	if (temp>c){
	  // char a[80];
	  // sprintf(a,"Cost increasing j:%d i:%d c: %f temp:%f",j,i,c,temp);
	  // mexErrMsgTxt(a);
	  for (int i = 0; i !=nodes; ++i)
	    label[i]=oldlabel[i];	
	  temp=c;
	}
	else {
	  if (temp<c){
	    printf ("Cost decreasing to %f\n",temp);
	    i=0;
	  }
	  for (int i = 0; i !=nodes; ++i)
	    oldlabel[i]=label[i];
	  c=temp;
	}
	++i;
	if (i==hyp){
	  c2=c;
	  return;
	} 
      }
    }
}
  void solve(){
    fast_solve();
}    
};
