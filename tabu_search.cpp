#include <stdio.h>
#include <iostream>
#include <vector> 
#include <algorithm>
#include <utility>
//local search al;gorith and testiNG
//begin local_search([Xi],V*,t)
using namespace std;

void fchange(double **qubo,int *sol,int bit,int qubo_size,double *dxo);
double local_search(int *sol, int qubo_size, double **qubo, int *flipfreq,double*dxo,double V);
double evaluate(int *sol,double **qubo,int qubo_size);

//local search function
double local_search(int *sol, int qubo_size, double **qubo, int *flipfreq,double*dxo,double V)
{
  bool improve=true;
  //int t=0; //iterations
  
  //improved:=false /* set flag that marks whether an improved solution found or not */
  while(improve)
    {
      improve=false;
      for(int k=0;k<qubo_size;++k) /* examine all variables */
	{
	  //V=evaluate(sol,qubo) //V:= qijXiXj /* evaluate new solution */
	  //calculate the change in objective fun after flipping
	  fchange(qubo,sol,k,qubo_size,dxo);
	  if (dxo[k]<0)   // then /* check for improved solution */
	    {
	      V=V+dxo[k]; /* record improved solution */
	      sol[k]=1-sol[k]; /* make the move for variable k */
	      flipfreq[k]++;  //update flip frequency
	      improve=true; /* set improved flag */
	     // cout<<k<<" flipped "<<V<<endl;
	    }
	}
    }
  return V;
}

double evaluate(int *sol,double **qubo,int qubo_size) //use for first or last time
{
  double energy1=0.0;
  for(int i=0;i<qubo_size;++i)
	{
	  for(int j=0;j<qubo_size;++j)
	    {
	      energy1+=qubo[i][j]*sol[i]*sol[j];
	    }
	}
  return energy1;
 // cout<<energy1<<"in evaluate" <<endl;
}
//function to make matrix lower triangular to make use of fast function calcuations for 1-flip move
void LowerTriangulize(double**qubo,int qubo_size)
{
		//make qubo matrix lower triangiular
	for (int i = 0; i <qubo_size; ++i)
	{
		for (int j = 0; j < qubo_size; ++j)
		{
			if(i<j)
			{
				qubo[i][j]=qubo[i][j]+qubo[j][i];
			}else if(i>j)
			{
			qubo[i][j]=0.0;
			}
		}
	}
}

//Neighboorhood based 1-flip moves
//fast calculations of objective function
void fchange(double **qubo,int *sol,int bit,int qubo_size,double*dxo)
{
  
  //calculate dxo
  dxo[bit]=qubo[bit][bit];
  for (int j = 0; j < qubo_size; ++j)
    {
      if(j!=bit &&sol[j]==1) //ignore zero values instead of multiplication
	{
	  dxo[bit]+=qubo[bit][j]+qubo[j][bit]; //row+column
	}
		
    }
  if(sol[bit]==1) dxo[bit]=-dxo[bit];  //sign change delta
  //flipping
  //sol[bit]=1-sol[bit];
  //update the flipfrequecy
  //flipfreq[bit]+=1;
  
}


//Xi:=0 i=1,...,n /* set the starting solution */
//V*:=0 /* initialise best solution value */
//begin


//basic tabu search follows Beaseley 1998 method
//sol current solution or the starting solution best known
//nTabu  tbau_tenure  L* is the tabu tenure value
//energy V* energy for current solution
//iter_max maximum number of iterations //T*
//TabuK[i]  Li is the tabu value for variable i, where Li=0 if the variable is not tabu

double tabu_search(int *sol,int *best,int qubo_size, double **qubo, double *dxo, int *flipfreq,int iter_max, int *TabuK, int *TabuTenure,double energy,int its)
{
 // cout<<"tabu_started"<<endl;
    //for(int i=0;i<qubo_size;++i) TabuK[i]=0;//Li:=0 i=1,...,n /* initialise tabu values */
  //sol Xi is the value for variable i in the current solution
  //energy V* is the best solution value found so far
  double benergy;//benergy V** is the best solution value associated with a neighbour of [Xi]
  int kmove=0; //flip variable
	energy=evaluate(sol,qubo,qubo_size);
//  for(int i=0;i<qubo_size;++i) TabuK[i]=0;//Li:=0 i=1,...,n /* initialise tabu values */
  benergy=1E100;   //V**:=-∞ /* initialise best neighbour value */
  its=0;
  while(its<iter_max)  /* iter_max* iterations in all */
    {
	for(int k=0;k<qubo_size;++k) /* for all variables */
	{
	  if(TabuK[k]==0)/* examine all non-tabu variables */
	    {
	      (its)++; /* increment iteration counter */
	      fchange(qubo,sol,k,qubo_size,dxo);   // get change in objective function for flipping kth var without flipping
	      energy+=dxo[k];	/* evaluate new solution */
	      if (dxo[k]<0)  //V > V* then /* check for improved solution */
		{
		  sol[k]=1-sol[k]; //flip
	//	  cout<<k<<" th variable flipped"<<energy<<endl;
		 flipfreq[k]++; /* record variable associated with the move */
		  kmove=k;
		 // cout<<"energy before local search"<<energy<<endl;
		  energy=local_search(sol, qubo_size, qubo,flipfreq, dxo,energy);  //make a local search irespective of tabu move
		  break; //go to done:
		}
	     else if(energy < benergy)  /* check for improved neighbouring solution */
		{
		  sol[k]=1-sol[k]; /* make the move */
	//	  cout<<k<<" th variable flipped neighbourhood "<<energy<<endl;
		  flipfreq[k]++;//record
		  kmove=k;
		  benergy=energy; /* record solution value */
		}//end if
	    }
	}//end for and if 
      //done 
      for(int i=0;i<qubo_size;++i) TabuK[i]=max(0,TabuK[i]-1);   // reduce all tabu values by one 
      TabuK[kmove]=TabuTenure[kmove]; /* tabu the chosen variable */
    }//end while
  // copy over the best solution
  for (int i = 0; i < qubo_size; i++) best[i] = sol[i];
  energy=evaluate(best,qubo,qubo_size);
	//cout<<"tabu_ended "<<its<<endl;
  
  return energy;
  
}//end
