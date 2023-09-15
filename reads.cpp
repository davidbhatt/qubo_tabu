#include <stdio.h>
#include <iostream>
#include <tuple>
#include <vector> 
#include <algorithm>
using namespace std;
tuple<int, int, int> read_qubo(FILE *inFile)
{
  char a[1000];
  int nmin, nmax, nNodes;
  if (inFile == NULL)
    {
      fprintf(stderr,
	      "\n\t%s error -- no input file (-i option) specified"
	      "\n\n",
	      "qbsolv");
      exit(9);
    }
  fscanf(inFile,"%s",a); //ignore first line
  
  int f1,f2;
  double f3;
  vector <int> entries;
  while (fscanf(inFile, "%d%d%lf\n", &f1, &f2,&f3) == 3)
    {
      //  printf("%d %d %lf\n", f1,f2,f3);
      entries.push_back(f1);
    }
  nmax=*max_element(entries.begin(), entries.end());
  nmin=*min_element(entries.begin(), entries.end());
 nNodes=nmax-nmin+1;
 //cout<<"found nodes "<<nNodes<<"("<<nmin<<","<<nmax<<")"<<endl;
 return make_tuple(nmin, nmax, nNodes);
}
void fill_qubo(FILE *inFile,double **val,int nmin)
{
  char a[100];
  int f1,f2;
  double f3;
  rewind(inFile);
  fscanf(inFile,"%s",a); //ignore first line
  while (fscanf(inFile, "%d%d%lf\n", &f1, &f2,&f3) == 3)
    {
      //  printf("%d %d %lf\n", f1,f2,f3);
      //  *(*(val + f1 - 1) + f2 - 1)=f3
      val[f1-nmin][f2-nmin]=f3/2;
	  val[f2-nmin][f1-nmin]=f3/2;
  }
  
}
