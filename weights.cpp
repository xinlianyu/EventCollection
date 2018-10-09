/*
  EdgeWeight Functions
  Brian Dean
  1998
  Revised by Song Gao May 2004
  Change weights to be stochastic and time-dependent
*/


#define EDGEWEIGHTS_SOURCE
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parser.h"
#include "network.h"
#include "weights.h"
#include <math.h>
#include "heapSort.h"

/*
void EdgeWeights::calexpected(Network *N, int EV_start, int EV_end, int **event, int start_time) 
{
  	//conditional expected travel times

	int		i, t, m, r, rr;
  	float 	tt;
  	double 	conditionalProb;
  
  	expectedtt = new float*[N->NumLinks()];
  	for (m = 0; m < N->NumLinks(); m++) 
    	expectedtt[m] = new float[_levels];
  
	conditionalProb = 0.0;
	for(i = EV_start; i < EV_end; i++) 
  		conditionalProb += prob[*(*(event + start_time) + i)];
  
  	for(t = start_time; t < _levels; t++)
    	for(m = 0; m < N->NumLinks(); m++)
    		{
      		tt = 0.0;
      		for(r = EV_start; r < EV_end; r++) 
      			{
  				rr = *(*(event + t) + r);
    			tt += *(*(realization + t * N->NumLinks() + m) + rr) * prob[rr];
      			}
      		expectedtt[m][t] = tt/conditionalProb;
    		}
}
*/

void EdgeWeights::calexpectedBP(Network *N) 
{
	long	m, r, links = N->NumLinks();
	int		count, *sorted, i, time;
	bp		x, next;
    float   val;
	//float   *sorting;
	double   *sorting;

	for (m = 0; m < links; m++)
		{
		if (det_flag[m] == 0) //this is a stochastic link
			{
			//calculate the largest possible number of breaking points for expected travel time
			count = 0;
			for (r = 0; r <	_numreal; r++)
				count += BP[m][r].size();
			//sorting = new float[count];
			sorting = new double[count];
			sorted = new int[count];

			//read breaking points from all support points and sort
			count = 0;
			for (r = 0; r < _numreal; r++)
                for (i = 0; i < BP[m][r].size(); i++)
                    {
					sorting[count] = BP[m][r][i].first;
                    sorted[count] = count;
                    count++;
                    }
      		heapSort(sorting, sorted, count); 
			
			//calculate expected travel times at each breaking point
			for (i = 0; i < count; i++)
				{
				if (i > 0 && sorting[sorted[i]] == sorting[sorted[i - 1]])
					continue;

				time = sorting[sorted[i]];
                val = 0;
				for (r = 0; r < _numreal; r++)
					val += weight_link(m, time, N, &BP[m][r]) * prob[r];

                x.first = time;
				x.second = val;
				ExpectedBP[m].push_back(x);
				}
            delete sorting;
            delete sorted;
			}
        }
}

/*
void EdgeWeights::calexpected(Network *N) 
{
  	//unconditional expected tt for all time periods
	int		i, t, m, r, links = N->NumLinks();
  	float 	tt;
  
  	expectedtt = new float*[links];
  	for (m = 0; m < links; m++) 
    	expectedtt[m] = new float[_levels];
  
  	for (t = 0; t < _levels; t++)
    	for(m = 0; m < links; m++)
    		{
      		tt = 0.0;
      		for(r = 0; r < _numreal; r++) 
  				tt += realization[t * links + m][r] * prob[r];
      		expectedtt[m][t] = tt;
    		}
}


void EdgeWeights::calmarginal(Network *N)  //unconditional marginal distribution; different links can have different number of support points
{ 
  	int		m, r, t, i, curr, next, numprob, p;
  	int*	forMarginal = new int[_numreal];
	float*	temp_marginalprob = new float[_numreal];
	float*	temp_marginaltt = new float[_numreal];	
	int		numlinks = N->NumLinks();

  	marginalprob = new float**[numlinks];
	marginaltt 	= new float**[numlinks];
	nummarginalprob = new int*[numlinks];
	for (m = 0; m < numlinks; m++)
		{
		marginalprob[m] = new float*[_levels];
		marginaltt[m] = new float*[_levels];
		nummarginalprob[m] = new int[_levels];
		}

  	for (r = 0; r < _numreal; r++) 
  		forMarginal[r] = r;

  	for (t = 0; t < _levels; t++)
  		{
    	for(m = 0; m < numlinks; m++)
    		{
      		heapSort(realization[t * numlinks + m], forMarginal, _numreal); 

      		//check the sorting result
      		//for(int i=0;i<NUM_REAL;i++) cout<<*(*(realization+t*NUM_ARCS+m)+forMarginal[i])<<"\t";
      		//cout<<endl;

      		for(i = 0; i < _numreal; i++) 
      			temp_marginalprob[i] = 0;
      		curr = 0;
      		next = 0;
      		numprob = 1;
      		for (r = 0; r < _numreal - 1; r++) 
      			{
				if (realization[t * numlinks + m][forMarginal[r]] != realization[t * numlinks + m][forMarginal[r + 1]])
					{
	  				next = r + 1;
	  				numprob++;

	  				for (i = curr; i < next; i++)
	    				temp_marginalprob[numprob - 2] += prob[forMarginal[i]];
	  				temp_marginaltt[numprob - 2] = realization[t * numlinks + m][forMarginal[curr]];
      
	  				curr = next;
					}
      			}

      		for(i = curr; i < _numreal; i++)
				temp_marginalprob[numprob - 1] += prob[forMarginal[i]];
      		temp_marginaltt[numprob - 1] = realization[t * numlinks + m][forMarginal[curr]];
      
      		nummarginalprob[m][t] = numprob;
      		marginaltt[m][t]= new float[numprob];
      		for (p = 0; p < numprob; p++) 
      			{
				marginalprob[m][t][p] = temp_marginalprob[p];
				marginaltt[m][t][p] = temp_marginaltt[p];
      			}
    		}//end of arc loop (m)
  		}//end of time loop (t)

  	delete[] forMarginal;
  	delete[] temp_marginalprob;
  	delete[] temp_marginaltt;
}
 
float *EdgeWeights::pathtimes(long origIndex, long destIndex, long pathIndex, long t, Network *N, long EV_start, long EV_end, long flag, long **event) {

  //flag = 0: unconditional
  
  float *result = new float[2];
  result[0]=0.0; result[1]=0.0;
  
  float *traveltime = new float[EV_end-EV_start];
  float linktime;
  int currtime;
  
  int rr;
  
  vector< int > links = N->path2links(origIndex, destIndex, pathIndex);
  
  //printf("size of links: %d\n", links.size());

  if (EV_start==0 && EV_end==_numreal && flag == 0) {
    //printf("unconditional case!\n");
    for (int r=0;r<_numreal;r++) {
      traveltime[r] = 0.0;
      currtime = t;

      for (vector< int >::iterator it=links.begin();it!=links.end();it++) {
	//printf("currtime: %d\n", currtime);
	linktime = realization[currtime*N->NumLinks()+(*it)][r];
	traveltime[r] += linktime;
	//if (pathIndex == 2)
	//printf("travel time of link %d at time %d: %f\n", (*it), currtime, linktime);
	currtime = (int)(currtime+linktime+0.5);
	currtime = currtime<(_levels-1)?currtime:(_levels-1);
      }
      //printf("current prob: %f\n", prob[r]);
      result[0] += traveltime[r]*prob[r];
      //printf("mean as just computed: %f\n", result[0]);
    }

    for (int r=0;r<_numreal;r++)
      result[1] += pow(traveltime[r-EV_start]-result[0], 2)*prob[r];
  
  }
  else {
    //printf("conditional case!\n");
      for (int r=EV_start;r<EV_end;r++) {
	traveltime[r-EV_start] = 0.0;
	rr = event[t][r];
	currtime = t;

	for (vector< int >::iterator it=links.begin();it!=links.end();it++) {
	  linktime = realization[currtime*N->NumLinks()+(*it)][rr];
	  traveltime[r-EV_start] += linktime;
	  currtime = (int)(currtime+linktime+0.5);
	  currtime = currtime<(_levels-1)?currtime:(_levels-1);
	}
	result[0] += traveltime[r-EV_start]*prob[rr];
      }

      for (int r=EV_start;r<EV_end;r++) {
	rr = event[t][r];
	result[1] += pow(traveltime[r-EV_start]-result[0], 2)*prob[rr];
      }
   
  }
  
  //printf("out of EdgeWeights::pathtimes()...\n");

  return result;
}

float EdgeWeights::pathtime(int origIndex, int destIndex, int pathIndex, int t, Network *N, int r) {

  float tt = 0.0;

  float linktime;
  int currtime;
  
  vector< int > links = N->path2links(origIndex, destIndex, pathIndex);

  currtime = t;

  for (vector< int >::iterator it=links.begin();it!=links.end();it++) {
    linktime = realization[currtime*N->NumLinks()+(*it)][r];
    tt += linktime;

    currtime = (int)(currtime+linktime+0.5);
    currtime = currtime<(_levels-1)?currtime:(_levels-1);
  }
  
  return tt;
}


void EdgeWeights::write_weights(char *filename, Network *_N) 
{
	FILE *fp;
	fp = fopen(filename, "r");
	if (!fp) 
		fp = fopen(filename, "w");
	else 
		{
		fclose(fp);
		fp = fopen(filename, "a");
		}

	if (!fp) 
		{
		printf( "Error: can't open %s\n", filename );
		exit(1);
		}

	int linkIndex;
	//fprintf(fp, "[Policy] : %d\n"); 
	for(int t = 0;t < _levels; t++) 
		for(int m = 0;m < _N->NumLinks(); m++) 
			{
		  	linkIndex = _N->linkIDtoIndex(_N->linkIdbyUser(m));
		  	for (int j = 0;j < _numreal; j++) 
				fprintf(fp, "%f ", this->weight_link(linkIndex, t, j, _N)*_agglength);     
		  	fprintf(fp, "\n");
			}

	fprintf(fp, "\n");

	fclose(fp);
}
*/

