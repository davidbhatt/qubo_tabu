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

using namespace std;

bool present(vector <int*>vec,int *key, int n);
void perturbation(int *sol,int*solp,int r,vector <int> elitefreq,int *flipfreq,int qubo_size,double beta, int gamma);
double evaluate(int *sol,double **val,int nNodes); //use for first or last time
double tabu_search(int *sol,int *best,int qubo_size, double **qubo, double *dxo, int *flipfreq,int iter_max, int *TabuK,int *TabuTenure,double energy,int its);
void random_perturbation(int *sol,int*solp,int qubo_size,int gamma);
void printsol(int *sol, int qubo_size);
//Algorithm 1 Diversification-Driven Tabu Search (D2TS) for UBQP
//1: Input: Q matrix  ....qubo
//2: Output: best ->S*: the best solution found so far
double D2ts_tabu(int *best,int *sol,double **qubo,int qubo_size,int R,double beta, int gamma,double max_time,int iter_max,int c,int its)
{
  double energy;
  double etime=0.0;
  int *worst;  //worst solution
  worst= new int [qubo_size+1];
  //3: Set EliteSol = {}, r = 0, EliteFreq(i ) = 0, i = 1, . . . , n
	vector <int*> elite(R); //elisesol
	for(int i = 0; i < R;i++) elite[i]=new int [qubo_size+1];
	
  vector <int> elitefreq(qubo_size,0); //initilizae vector with zero value
  //vector <int> flipfreq(qubo_size,0); //initilizae vector with zero value
  vector <double> elite_energy;
  ///initialization for basic tabu search algorith
  double *dxo;
  dxo= new double [qubo_size+1];
  int *flipfreq;
  flipfreq=new int [qubo_size+1];
  energy=evaluate(sol,qubo,qubo_size);
  //local_search(sol, nNodes, val,flipfreq, dxo,energy);
  
  int *TabuTenure = new int [qubo_size+1];
  for(int i=0;i<qubo_size;++i) TabuTenure[i]=c+rand()%10+1;
  int *TabuK = new int [qubo_size+1];
  for(int i=0;i<qubo_size;++i) TabuK[i]=0;//Li:=0 i=1,...,n /* initialise tabu values */
  /*******************************************/
  //Algorithm start
  //first loop
    int r=0; 
	bool isFound;
	//cout<<"first loop started"<<endl;	

  //5: while r < R do
  int isp;
  while(r<R)
    {
      //6: S* = Tabu_Search(S0)
	 // cout<<r<<" running"<<endl;
	  //printsol(sol,qubo_size);
      energy=tabu_search(sol,best,qubo_size, qubo, dxo, flipfreq,iter_max, TabuK,TabuTenure,energy,its);   // store best in best vector
	  //printsol(sol,qubo_size);
	  //check addition
	 /* cout<<"printing elite solutions"<<endl;
	  for(int i=0;i<r;++i)
	  {
		  cout<<i<<" th solution:";
		  printsol(elite[i],qubo_size);
		  cout<<endl;
	  }	
	  cout<<"printing current solution"<<endl;
	  printsol(best,qubo_size);
	  */
      isFound = present(elite,best,qubo_size);
	  //cout<<"1st loop sol found in elite "<<isFound<<endl;
      if(!isFound) //7: if S* is not in EliteSol then
		{
		//cout<<"1st loop sol not found in elite "<<isFound<<endl;
		//8: Insert Sst into EliteSol: EliteSol = EliteSol + {Sst}
		for(int j = 0; j < qubo_size;++j) elite[r][j]=best[j];
		elite_energy.push_back(energy);//store function v alue also
		r++; //9: r = r + 1
		//10: EliteFreq = EliteFreq + S*
		for (int i = 0; i < qubo_size; ++i) 
			{if(best[i]==1) elitefreq[i]+=best[i];}
		}//11: end if
      //12: Randomly select a solution S’ from EliteSol
      if(r==1)
	  {isp=0;}
	else
	  {
		isp=rand()%r;  //index of sprime
      }
		  //13: S0 = Perturbation_Operator(S’)
	  //use best for storing S'
	  for(int i=0;i<qubo_size;++i)best[i]=elite[isp][i];
	  //cout<<"best sol"<<endl;
	  //printsol(best,qubo_size);
	 // cout<<"selection "<<isp<<" length of elite "<<r<<endl;
     perturbation(sol,best,r,elitefreq,flipfreq,qubo_size,beta, gamma);
	// cout<<"difference in perturbed sol"<<endl;
	 //printsol(best,qubo_size);
	 // for(int i=0;i<qubo_size;++i)cout<<best[i]-sol[i];
	 // cout<<endl;
	//random_perturbation(sol,best,qubo_size,qubo_size);
    energy=evaluate(sol,qubo,qubo_size);  //update function value
	//cout<<"sol after perturbations"<<endl;
	//printsol(sol,qubo_size);
    }//14: end while
  
  //second loop
  //timer start
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  cout<<"second loop started"<<endl;
  //15: while Stop condition is not met do
  int it=0;
  //srand(0);
  while(etime<max_time)
    {
      //16: Randomly select a solution S’ from EliteSol
      int isp=rand()%R;  //index of sprime
	 // cout<<isp<<" selected"<<endl;
      // 17: S0 = Perturbation_Operator(S’)
  	  for(int i=0;i<qubo_size;++i)best[i]=elite[isp][i];
      perturbation(sol,best,R,elitefreq,flipfreq,qubo_size,beta, gamma);
      energy=evaluate(sol,qubo,qubo_size);  //update function value calcuate directly as many flippings in perturbation step
      energy=tabu_search(sol,best,qubo_size, qubo, dxo, flipfreq,iter_max, TabuK,TabuTenure,energy,its);   //18: S* = Tabu_Search(S0)
      //19: Sw = The worst solution in EliteSol in terms of solution quality
      //get the worst solution
	  int maxEI = std::max_element(elite_energy.begin(),elite_energy.end()) - elite_energy.begin();      
      worst=elite[maxEI];
	  double fsw=elite_energy[maxEI];
      //20: if S* is not in EliteSol and f (S*) < f (Sw) then
      if(!present(elite,best,qubo_size) && energy<fsw)
	{
	  //elite[maxEI]=best;  //21: EliteSol = EliteSol + {S∗} − {Sw}
	  for(int j = 0; j < qubo_size;++j) elite[maxEI][j]=best[j];
	  elite_energy[maxEI]=energy;  //function value update
	  //22: EliteFreq = EliteFreq + S∗ − Sw
	  for (int i = 0; i < qubo_size; ++i) 
	    {elitefreq[i]+=best[i]-worst[i];}
	//cout<<maxEI<<" updated"<<endl;
	 //for(int j = 0; j < R;++j)cout<<elite_energy[j]<<" ";
	//cout<<endl;
	}//23: end if
      //end time
      std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
      etime=float(std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count())/1000.0;  //in milli seconds
	  it++;
	  //    int minEI = std::min_element(elite_energy.begin(),elite_energy.end()) - elite_energy.begin();
	 // cout<<it<<" time elasped= "<<etime<<" ms f(best_solution) so far = "<<elite_energy[minEI]<<endl;

	}//24: end while
  //get the best solutions so far found
    int minEI = std::min_element(elite_energy.begin(),elite_energy.end()) - elite_energy.begin();
	//for(int j = 0; j < R;++j)cout<<elite_energy[j]<<" ";
	//cout<<endl;
	for(int j = 0; j < qubo_size;++j) best[j]=elite[minEI][j];
  //calculate the energy for last time 
  energy=evaluate(best,qubo,qubo_size);
  return energy;
}

//function to check if the chosen value solution array key is present in the master list of elite solutions vec
bool present(vector <int*> vec,int *key, int n)
{
  bool ans;
  bool v,vi;
  if(vec.size()==0)
  {ans=false;}
else
{
  for (const auto &item : vec) {
    v=true;
    for(int i=0;i<n;++i)
      {
	if (item[i] == key[i])  vi=true; else vi=false;
	v=v&&vi;
	if(!v) break;
      }
    if(v) {ans=true;break;} else {ans=false;}
  }
}
  return ans;
}



//function for perturbation of solution
//sol the S0 peerturbed solution
//solp S'
//int r current number of solutions recorder 0<r<R
//beta Frequency-related weight in perturbation scoring
//gamma first gamma critical values to be flipped
void perturbation(int *sol,int *solp,int r,vector <int> elitefreq,int *flipfreq,int qubo_size,double beta, int gamma)
{
  vector <double>score(qubo_size,0);
  //determine which bit to perturbation
  //score for all variables
  int maxfreq=*max_element(flipfreq , flipfreq + qubo_size);
  //for(int i=0;i<qubo_size;++i)cout<<elitefreq[i]<<" ";
  //cout<<endl;
  for(int i=0;i<qubo_size;++i)score[i]=double(elitefreq[i])*double(r-elitefreq[i])/double(r*r)+beta*(1-double(flipfreq[i])/double(maxfreq));
	//find maximum score  and variable to be flipped
 // double max_score=*std::max_element(score.begin(),score.end());
  //arrange in descending order of score
  //store the indices in V vector
  vector<int> V(qubo_size);
  std::iota(V.begin(),V.end(),0); //Initializing
  //sort in decreasing value of scores
  sort( V.begin(),V.end(), [&](int i,int j){return score[i]>score[j];} ); //get variable list in descending order of score
   // cout<<"scorer:"<<endl;
  //for(int i=0;i<qubo_size;++i) cout<<score[V[i]]<<" ";
    //cout<<endl;
  //perturbation step
  //cout<<"perturbed indices:"<<endl;
  //for(int i=0;i<gamma;++i) cout<<V[i]<<" "<<elitefreq[i]<<" ";
  //cout<<endl;
  //flip first gamma variables form V list
  for(int i=0;i<gamma;++i)
    {
		  sol[V[i]]=1-solp[V[i]];
		flipfreq[i]++;	  
    }
	/*cout<<"{";
	for(int i=0;i<qubo_size;++i) cout<<elitefreq[i]<<", ";
	cout<<"}"<<endl;
	cout<<"before perturbation :";
	for(int i=0;i<qubo_size;++i) cout<<solp[i];
	cout<<endl;
	
	cout<<"difference in  perturbation in loop:"<<endl;
	for(int i=0;i<qubo_size;++i) cout<<solp[i]-sol[i];
	cout<<endl;
	*/
	
}
void random_perturbation(int *sol,int*solp,int qubo_size,int gamma)
{
  
  //store the indices in V vector
  vector<int> V(qubo_size);
  for(int i=0;i<qubo_size;++i) V.push_back(i);
  std::random_shuffle(V.begin(),V.end());
  
 //perturbation step
  //flip first gamma variables form V list
  for(int i=0;i<gamma;++i)
    {
      sol[V[i]]=1-solp[V[i]];
    }

}
	void printsol(int *sol, int qubo_size)
	{
		for(int i=0;i<qubo_size;++i) cout<<sol[i];
		cout<<endl;
	}
