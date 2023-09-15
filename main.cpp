//Algorithm 1 Diversification-Driven Tabu Search (D2TS) for UBQP
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <vector> 
#include <algorithm>
#include <tuple>
#include <utility>
#include <bits/stdc++.h>
#include <chrono>
#include <random>


using namespace std;

tuple<int, int, int> read_qubo(FILE *inFile);
void fill_qubo(FILE *inFile,double **val,int nmin);
void LowerTriangulize(double**qubo,int qubo_size);
double D2ts_tabu(int *best,int *sol,double **qubo,int qubo_size,int R,double beta, int gamma,double max_time,int iter_max,int c,int its);


int main(int argc, char *argv[])
{
  char *inFileName = NULL;
  FILE *inFile = NULL;
  int iter_max=0;
  int time_out=0; //milli seconds
	int c=0; //tabu_tenure constant
	srand(time(0));
  //parse arguments
  for(int i=1;i<argc;++i)
	{
	 // cout << argv[i] << "\n";
	  string s=argv[i];
	  //cout<<s<<endl;
	  if(s.compare("-i")==0)  //argv[i]=="-"&&argv[i][1]=="i") 
	    {
	      inFileName=argv[i+1];
	    }
	  if(s.compare("-it")==0)  //argv[i]=="-"&&argv[i][1]=="i") 
	    {
	      iter_max=atoi(argv[i+1]);
	    }
	
	if(s.compare("-t")==0)  //argv[i]=="-"&&argv[i][1]=="i") 
	    {
	      time_out=atoi(argv[i+1]);
	    }

		//tabu tenure
	if(s.compare("-tt")==0)  //argv[i]=="-"&&argv[i][1]=="i") 
	    {
	      c=atoi(argv[i+1]);
	    }

	} 
   
  //read file
  inFile = fopen(inFileName, "r");
  
  int nmin,nmax,nNodes;
	tie(nmin,nmax,nNodes) = read_qubo(inFile);
	
	cout<<"found nodes "<<nNodes<<"("<<nmin<<","<<nmax<<")"<<endl;
	//fill qubo matrix
	double **val;
	val= new double *[nNodes+1];
	for(int i = 0; i <nNodes+1; i++) val[i] = new double[nNodes+1];
	fill_qubo(inFile,val,nmin);
	fclose(inFile);
//set paramters  default values if not provided
    if(iter_max==0) iter_max=20*nNodes; //T*:=max(500000,5000n) alpha no of moves /* set maximum number of iterations */
	if(c==0) c=max(20,int(nNodes/100)); /* set tabu tenure constant */
	if(time_out==0) time_out=200;
//display paramters	
  cout<<"tabu_search will run using supplied file "<<inFileName <<" and "; 
  cout<<"max_iterations ="<<iter_max<<" or time_out = "<<time_out<<" ms"<<endl;
  ///start of code
	//make qubo matrix lower triangiular
	LowerTriangulize(val,nNodes);
	//initiate solution
	int *sol;
	sol= new int [nNodes+1];
	int *best;
	best= new int [nNodes+1];
 //4: Randomly generate an initial solution S0
	//srand(time(NULL)); 
	for (int i = 0; i < nNodes;i++) 
	{
	  sol[i]=rand()%2;
	}
	//cout<<endl;
	
	//srand(time(0));
///constants values from the paper
	int R=8;  //R constant length of elite solution 
	int gamma=int(nNodes/4);
	double beta=0.3;
	double energy;  //function value
	int its=0;  //iteration count for tabu
	
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now(); //start timer

	energy= D2ts_tabu(best,sol,val,nNodes,R,beta, gamma,time_out,iter_max,c,its);
	
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now(); //end timer

	double etime=float(std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count())/1000.0;  //in milli seconds

	//cout<<"optimal solution found is ";
    //printsol(best,nNodes);
		
	cout<< "Global minimum: "<<energy<<endl;
	cout<<" total time taken"<<etime<<endl;	
		
	return 0;
}

