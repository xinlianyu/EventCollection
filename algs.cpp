/*
  Network Algorithms
  Song Gao
  May 2004
  Song Gao, Added NOI, Feburary 2007
  Revised by Jing Ding 2012
  */


#define ALGS_SOURCE

#include <assert.h>
#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <vector>
#include <queue>
#include <iostream>
#include <fstream>
#include "heapSort.h" //heap from Seong
#include "myqueue.h" //queue from Seong
#include "network.h"
#include "weights.h"
#include "algs.h"

#include <algorithm> 

#include <stdlib.h>
#include <time.h>
#include <math.h>

#define EPSILON 0.05
#define alpha 1
#define startingmean 300

#define disutility 1

using namespace std;

float min(float a, float b) {
	return a < b ? a : b;
}

float max(float a, float b){
	return a > b ? a : b;
}

/*
void generate_event_collection(Network *N, EdgeWeights *w, int **event, int **eventTag, int **realizationToEvent, char *filename) {

register int NUM_ARCS, NUM_REAL, TIME_PERIOD;
float 	x, y;
int 	curr, next, r, t, m, k, j;

NUM_ARCS = N->NumLinks();
NUM_REAL = w->numreal();
TIME_PERIOD = w->levels();

float *currRealization = new float[NUM_REAL];
int	*currEvent = new int[NUM_REAL];
int	*sorted = new int[NUM_REAL];

//float **realization = w->realization;
//assert( w->minweight() > 0 );

//////////////////////////////////////////////////////////
//Generate event collections
//////////////////////////////////////////////////////////

for(t = 0; t < TIME_PERIOD; t++)
{
//cout<<"t = "<<t<<endl;

for (r = 0; r < NUM_REAL; r++)
{
if(t == 0)
{
*(*(event + t) + r) = r;
*(*(eventTag+t))=1;
}
else
{
*(*(event + t) + r) = *(*(event + t - 1) + r);
*(*(eventTag + t) + r) = *(*(eventTag + t - 1) + r);
}
}
*(*(eventTag + t) + NUM_REAL) = 1; //mark the end

for (int m = 0; m < NUM_ARCS; m++)
{
if (N->Disabled(m))
continue;
curr = 0;
next = 0;
for (int r = 1; r < NUM_REAL + 1; r++)
{
if (*(*(eventTag + t) + r) == 1)
{
next = r; //determine the range of the current event collection

for (int i = curr; i < next; i++)
{
currEvent[i - curr] = event[t][i];
//currRealization[i - curr] = realization[t*NUM_ARCS + m][event[t][i]];
if (w->det_flag[m])
{
currRealization[i - curr] = w->weight_link(m, t, N, &(w->BP[m][0])) * 60;
}
else
{
currRealization[i - curr] = w->weight_link(m, t, N, &(w->BP[m][event[t][i]])) * 60;
}

sorted[i - curr] = i - curr;
}

heapSort(currRealization, sorted, next - curr);


for (int i = curr; i < next; i++)
event[t][i] = currEvent[sorted[i - curr]];

//heapSort(*(w->realization + t * NUM_ARCS + m) + curr, *(event + t) + curr, next - curr);
for (int j = curr + 1; j < next; j++)
{
//x = *(*(w->realization + t * NUM_ARCS + m) + *(*(event + t) + j - 1));
//y = *(*(w->realization + t * NUM_ARCS + m) + *(*(event + t) + j));

if (w->det_flag[m])
{
x = w->weight_link(m, t, N, &(w->BP[m][0])) * 60;
}
else
{
x = w->weight_link(m, t, N, &(w->BP[m][event[t][j-1]])) * 60;
}

if (w->det_flag[m])
{
y = w->weight_link(m, t, N, &(w->BP[m][0])) * 60;
}
else
{
y = w->weight_link(m, t, N, &(w->BP[m][event[t][j]])) * 60;
}


if (fabs(x - y) / min (x, y) > EPSILON)
*(*(eventTag + t) + j) = 1; //find the breaking points with a tolerance
}
curr = next;
}
}
}//end of arc loop (m)

//Establish the mapping from realization index to event index
curr = 0;
next = 0;
for (r = 1; r < NUM_REAL + 1; r++)
if(*(*(eventTag+t)+r)==1)
{
next = r;
for(j = curr; j < next; j++)
*(*(realizationToEvent + t) + *(*(event + t) + j)) = curr;
curr=next;
}

}//end of time loop (t)

//write_event(N, w, event, filename);
}
*/


void generate_event_collection(Network *N, EdgeWeights *w, int dests, vector<pair<long, int>> orig_dep, int **event, int **eventTag, int **realizationToEvent, char *eventfile, char *filename, int start_time)
{

	register int NUM_ARCS, NUM_REAL, TIME_PERIOD;
	float 	x, y;
	int 	curr, next, m, j, r, t, k;
	//int time;

	NUM_ARCS = N->NumLinks();
	NUM_REAL = w->numreal();
	TIME_PERIOD = w->levels();
	int period_length = w->agglength();

	//float *currRealization = new float[NUM_REAL];
	double *currRealization = new double[NUM_REAL];
	int	*currEvent = new int[NUM_REAL];
	int	*sorted = new int[NUM_REAL];

	float *percentile = new float[NUM_ARCS];
	for (int m = 0; m < NUM_ARCS; m++)
	{
		percentile[m] = 0;
	}

	int	*transTT = new int[NUM_REAL];
	for (int r = 0; r < TIME_PERIOD; r++)
	{
		transTT[r] = 0;
	}

	/*int   **transTT;
	transTT = new int*[NUM_ARCS];
	for (int m = 0; m < NUM_ARCS; m++)
	{
		transTT[m] = new int[NUM_REAL];
	}

	for (int m = 0; m < NUM_ARCS; m++)
	{
		for (int r = 0; r < TIME_PERIOD; r++)
		{
			transTT[m][r] = 0;
		}
	}*/

	//vector<vector<float>> z;
	//z.push_back(z_m)
	//vector<float> z_m;
	//z[m]

	for (int m = 0; m < NUM_ARCS; m++)
	{
		vector<float> z;  //vector z stores travel times for each link
		for (int t = 0; t < TIME_PERIOD; t++)
		{
			for (int r = 0; r < NUM_REAL; r++)
			{
				if (w->det_flag[m] == 1)

				{
					//LinkTT[m][r] = w->weight_link(m, start_time + t* period_length, N, &(w->BP[m][event[t][0]]));
					z.push_back(w->weight_link(m, start_time + t* period_length, N, &(w->BP[m][0])));
				}
				else
				{
					//LinkTT[m][r] = w->weight_link(m, start_time + t* period_length, N, &(w->BP[m][event[t][r]]));
					z.push_back(w->weight_link(m, start_time + t* period_length, N, &(w->BP[m][r])));

				}
			}
		}

		//then calculate the 95 percentile of vector z for each link m.

		sort(z.begin(), z.end());
		percentile[m] = z[int(z.size())*0.99-1];  // the xth percentile

		z.clear();

	}

	/*ofstream out_file("e:\\percentile.dat", std::ofstream::app);

	for (int m = 0; m < NUM_ARCS; m++)
	{
	out_file << percentile[m] << "\r\n";
	}
	out_file << std::endl;
	out_file.close();*/


	//float **realization = w->realization; 
	//assert( w->minweight() > 0 );

	//////////////////////////////////////////////////////////
	//Generate event collections
	//////////////////////////////////////////////////////////

	//time = start_time + (TIME_PERIOD - 1) * period_length;
	//time = start_time + (TIME_PERIOD - 1);

	for (int t = 0; t < TIME_PERIOD; t++)
	{
		//cout<<"t = "<<t<<endl;
		for (int r = 0; r < NUM_REAL; r++)
		{
			event[t][r] = 0;
			eventTag[t][r] = 0;
		}
	}

	for (int t = 0; t < TIME_PERIOD; t++)
	{
		for (int r = 0; r < NUM_REAL; r++)
		{

			if (t == 0)
			{
				event[t][r] = r;
				eventTag[t][0] = 1;
			}
			else
			{
				event[t][r] = event[t - 1][r];
				eventTag[t][r] = eventTag[t - 1][r];
			}
		}
		eventTag[t][NUM_REAL] = 1; //mark the end


		for (int m = 0; m < NUM_ARCS; m++)
		{					
			if (w->det_flag[m] == 0)
			{
				for (int r = 0; r < NUM_REAL; r++)

				{
					if (w->weight_link(m, start_time + t* period_length, N, &(w->BP[m][r])) > percentile[m])
					{
						transTT[r] = 1;
					}
					else
					{
						transTT[r] = 0;
					}
				}

				curr = 0;
				next = 0;

				for (int r = 1; r < NUM_REAL + 1; r++)
				{
					if (*(*(eventTag + t) + r) == 1)
					{
						next = r; //determine the range of the current event collection

						for (int i = curr; i < next; i++)
						{
							currEvent[i - curr] = event[t][i];
							//currRealization[i - curr] = realization[t*NUM_ARCS + m][event[t][i]];
							//currRealization[i - curr] = w->weight_link(m, start_time + t, N, &(w->BP[m][0])) * 60.0;

							currRealization[i - curr] = transTT[i];

							//*(*(link[dest_link].cost_label + t) + r) = w->weight_link(dest_link, start_time + t * period_length, N, &(w->BP[dest_link][0]));

							sorted[i - curr] = i - curr;
						}

						heapSort(currRealization, sorted, next - curr);	//	
						for (int i = curr; i < next; i++)
							event[t][i] = currEvent[sorted[i - curr]];
						for (int j = curr + 1; j < next; j++)
						{
							int x, y;

							x = transTT[event[t][j - 1]];
							y = transTT[event[t][j]];

							if (x != y)
								*(*(eventTag + t) + j) = 1; //find the breaking points with a tolerance
						}
						curr = next;
					}

					/*ofstream out_file("e:\\transTT.dat", std::ofstream::app);

					for (int r = 0; r < NUM_REAL; r++)
					{
					out_file << transTT[r] << " \r\n";
					}
					out_file << std::endl;
					out_file.close();*/

				} //end of r loop

				/*ofstream out_file("e:\\transTT.dat", std::ofstream::app);

				for (int r = 0; r < NUM_REAL; r++)
				{
				out_file << transTT[r] << " ";
				}
				out_file << std::endl;
				out_file.close();*/

			}

		}//end of arc loop (m)
		

		/*ofstream out_file("e:\\eventTag99%.dat", std::ofstream::app);

		for (int r = 0; r < NUM_REAL + 1; r++)
		{
		out_file << *(*(eventTag + t) + r) << " ";
		}
		out_file << std::endl;
		out_file.close();*/


		/*ofstream out_file("e:\\event.dat", std::ofstream::app);

		for (int r = 0; r < NUM_REAL; r++)
		{
			out_file << *(*(event + t) + r) << " ";
		}
		out_file << std::endl;
		out_file.close();*/
		/*for (int m = 0; m < NUM_ARCS; m++)

		{
		delete[] currTT[m];
		}

		delete[] currTT;*/


		//Establish the mapping from realization index to event index
		curr = 0;
		next = 0;
		for (int r = 1; r < NUM_REAL + 1; r++)
			if (*(*(eventTag + t) + r) == 1)
			{
				next = r;
				for (int j = curr; j < next; j++)
					*(*(realizationToEvent + t) + *(*(event + t) + j)) = curr;//????
				curr = next;
			}

	}//end of time loop (t)

	/*ofstream out_file("e:\\eventTag.dat", std::ofstream::app);

	for (t = 0; t < TIME_PERIOD; t++)
	{
	for (int r = 0; r < NUM_REAL + 1; r++)
	{
	out_file << *(*(eventTag + t) + r) << " ";
	}
	}
	out_file << std::endl;
	out_file.close();*/

	delete percentile;
	delete transTT;
	write_event(N, w, dests, orig_dep, eventTag, eventfile);
}
//delete transTT;

//calculate static shortest paths in the static period for each possible realization
//void static_shortest_path( Network *N, EdgeWeights *w, int destIndex, int **eventTag, int **realizationToEvent, node_info *node, int start_time)


void static_shortest_path(Network *N, EdgeWeights *w, int destIndex, int **eventTag, int **realizationToEvent, link_info *link, int start_time)
{
	// Initialization and shortest paths calculation in static area starts

	register int TIME_PERIOD, NUM_NODES, NUM_ARCS, NUM_REAL;

	NUM_NODES = N->NumNodes();
	NUM_ARCS = N->NumLinks();
	NUM_REAL = w->numreal();
	TIME_PERIOD = w->levels();

	assert(destIndex >= 0 && destIndex < NUM_NODES);
	//assert( w->minweight() > 0.5 );

	int n, m, t, r;
	//int indegree; 
	long ilink, jlink, klink;
	int ntime;
	int inode, jnode;
	float temp_cost_label;
	int time;
	//int knode;
	int tnode;
	int T = 0;
	int block_number = 1440;
	long indegree, link_curr, link_next, dest_prev, dest_link;

	//time = start_time + (TIME_PERIOD - 1) * period_length;
	time = start_time + (TIME_PERIOD - 1);

	for (t = 0; t < TIME_PERIOD; t++)
	{
		for (r = 0; r < NUM_REAL; r++)
		{
			if (*(*(eventTag + t) + r) == 1)
			{
				//for(n = 0; n < NUM_NODES; n++)
				for (n = 0; n < NUM_ARCS; n++)
				{
					//if (N->Disabled(n))
					//continue;
					//*(*(node[n].cost_label + t) + r) = INT_INFINITY;
					//*(*(node[n].next_node + t) + r) = -1;
					*(*(link[n].cost_label + t) + r) = INT_INFINITY;
					*(*(link[n].next_node + t) + r) = -1;
				}
				//*(*(node[destIndex].cost_label + t) + r) = 0; 
				//*(*(node[destIndex].next_node + t) + r) = destIndex; //destination

				indegree = N->InDegree(destIndex);
				for (m = 0; m < indegree; m++)
				{
					dest_prev = N->LinkSource(destIndex, m);
					dest_link = N->LinkIndex(dest_prev, destIndex);
					if (N->Disabled(dest_link))
						continue;
					//*(*(link[dest_link].cost_label + t) + r) = 0; 
					*(*(link[dest_link].next_node + t) + r) = destIndex; //destination
					if (w->det_flag[dest_link] == 1)
					{
						//*(*(link[dest_link].cost_label + t) + r) = w->weight_link(dest_link, time, N, &(w->BP[dest_link][0]))*60;
						*(*(link[dest_link].cost_label + t) + r) = w->weight_link(dest_link, time, N, &(w->BP[dest_link][0]));

					}
					else
					{
						//*(*(link[dest_link].cost_label + t) + r) = w->weight_link(dest_link, time, N, &(w->BP[dest_link][r]))*60;
						*(*(link[dest_link].cost_label + t) + r) = w->weight_link(dest_link, time, N, &(w->BP[dest_link][r]));
					}
				}
			}
		}
	}

	for (r = 0; r < NUM_REAL; r++)
	{
		//T=0;
		myqueue SE;
		indegree = N->InDegree(destIndex);
		for (m = 0; m < indegree; m++)
		{
			dest_prev = N->LinkSource(destIndex, m);
			dest_link = N->LinkIndex(dest_prev, destIndex);
			if (N->Disabled(dest_link))
				continue;
			SE.enqueue(dest_link);
		}
		//SE.enqueue(destIndex);
		int eventindex = *(*(realizationToEvent + TIME_PERIOD - 1) + r);
		while (!SE.is_empty())
		{
			//jnode=SE.dequeue();
			jlink = SE.dequeue();
			jnode = N->LinkSource(jlink);
			indegree = N->InDegree(jnode);
			for (m = 0; m < indegree; m++)
			{
				T = 0;
				inode = N->LinkSource(jnode, m);
				link_curr = N->LinkIndex(inode, jnode);
				//if (w->weight_link(link_curr, time, N, &(w->BP[link_curr][r])) > block_number || w->weight_link(link_curr, time, N, &(w->BP[link_curr][r])) < 0)
				//if (w->weight_link(link_curr, time, N, &(w->BP[link_curr][r])) > block_number)
				//{
				//N->SetDisabled(link_curr, 1);
				//}

				if (N->Disabled(link_curr))
					continue;
				//knode = *(*(link[jlink].next_node + TIME_PERIOD - 1) + eventindex);
				//link_next = N->LinkInIndex(jnode, knode);				
				tnode = N->LinkDest(jlink);
				if (inode == tnode)
				{
					T = INT_INFINITY;
				}
				//temp_cost_label = w->weight(inode, jnode, TIME_PERIOD - 1, r, N) + *(*(node[jnode].cost_label + TIME_PERIOD - 1) + eventindex);
				if (w->det_flag[N->LinkIndex(inode, jnode)] == 1)
				{
					//temp_cost_label = w->weight(inode, jnode, time, N, &(w->BP[m][0])) + *(*(node[jnode].cost_label + TIME_PERIOD - 1) + eventindex);
					//temp_cost_label = w->weight(inode, jnode, time, N, &(w->BP[m][0])) + *(*(link[jlink].cost_label + TIME_PERIOD - 1) + eventindex) + T;
					//temp_cost_label = w->weight_link(link_curr, time, N, &(w->BP[link_curr][0]))*60 + *(*(link[jlink].cost_label + TIME_PERIOD - 1) + eventindex) + T;
					temp_cost_label = w->weight_link(link_curr, time, N, &(w->BP[link_curr][0])) + *(*(link[jlink].cost_label + TIME_PERIOD - 1) + eventindex) + T;

				}
				else
				{
					//temp_cost_label = w->weight(inode, jnode, time, N, &(w->BP[m][r])) + *(*(node[jnode].cost_label + TIME_PERIOD - 1) + eventindex);
					//temp_cost_label = w->weight(inode, jnode, time, N, &(w->BP[m][r])) + *(*(link[jlink].cost_label + TIME_PERIOD - 1) + eventindex) + T;
					//temp_cost_label = w->weight_link(link_curr, time, N, &(w->BP[link_curr][r]))*60 + *(*(link[jlink].cost_label + TIME_PERIOD - 1) + eventindex) + T;
					temp_cost_label = w->weight_link(link_curr, time, N, &(w->BP[link_curr][r])) + *(*(link[jlink].cost_label + TIME_PERIOD - 1) + eventindex) + T;

				}

				//*(*(link[dest_link].cost_label + t) + r) = (int)(w->weight_link(dest_link, t, N, &(w->BP[m][0]))*60+0.5);
				//*(*(link[dest_link].cost_label + t) + r) = (int)(w->weight_link(dest_link, t, N, &(w->BP[m][event[t][r]]))*60+0.5);

				//if(temp_cost_label < *(*(node[inode].cost_label + TIME_PERIOD - 1) + eventindex))
				if (temp_cost_label < *(*(link[link_curr].cost_label + TIME_PERIOD - 1) + eventindex))
				{
					//*(*(node[inode].cost_label + TIME_PERIOD - 1) + eventindex) = temp_cost_label;
					//*(*(node[inode].next_node+TIME_PERIOD-1)+eventindex)=jnode;
					*(*(link[link_curr].cost_label + TIME_PERIOD - 1) + eventindex) = temp_cost_label;
					*(*(link[link_curr].next_node + TIME_PERIOD - 1) + eventindex) = tnode;
					if (SE.have_this(link_curr) == 0)
						SE.enqueue(link_curr);
				}
				N->SetDisabled(link_curr, 0);
			}
		}
		SE.clear();
	}
}

/*
void static_shortest_path( Network *N, EdgeWeights *w, int destIndex, int **eventTag, int **realizationToEvent, link_info *link)
{
// Initialization and shortest paths calculation in static area starts

register int TIME_PERIOD, NUM_NODES, NUM_ARCS, NUM_REAL;

NUM_NODES = N->NumNodes();
NUM_ARCS = N->NumLinks();
NUM_REAL = w->numreal();
TIME_PERIOD = w->levels();

assert( destIndex>=0 && destIndex<NUM_NODES );
//assert( w->minweight() > 0.5 );

int n,m,t,r;
//int indegree;
long ilink,jlink,klink;
int ntime;
int inode,jnode;
float temp_cost_label;
//int time;
int knode,tnode;
int T=0;
long indegree,link_curr,link_next,dest_prev,dest_link;

//time = start_time + (TIME_PERIOD - 1) * period_length;
//time = start_time + (TIME_PERIOD - 1);

for(t = 0; t < TIME_PERIOD; t++)
for(r = 0; r < NUM_REAL; r++)
if(*(*(eventTag + t) + r) == 1)
{
//for(n = 0; n < NUM_NODES; n++)
for(n = 0; n < NUM_ARCS; n++)
{
//if (N->Disabled(n))
//continue;
//*(*(node[n].cost_label + t) + r) = INT_INFINITY;
//*(*(node[n].next_node + t) + r) = -1;
*(*(link[n].cost_label + t) + r) = INT_INFINITY;
*(*(link[n].next_node + t) + r) = -1;
}
//*(*(node[destIndex].cost_label + t) + r) = 0;
//*(*(node[destIndex].next_node + t) + r) = destIndex; //destination

indegree = N->InDegree(destIndex);
for(m = 0; m < indegree; m++)
{
dest_prev = N->LinkSource(destIndex, m);
dest_link= N->LinkIndex(dest_prev, destIndex);
if (N->Disabled(dest_link))
continue;
//*(*(link[dest_link].cost_label + t) + r) = 0;
*(*(link[dest_link].next_node + t) + r) = destIndex; //destination
if (w->det_flag[dest_link]==1)
{
*(*(link[dest_link].cost_label + t) + r) = (int)(w->weight_link(dest_link, TIME_PERIOD - 1, N, &(w->BP[m][0]))*60+0.5);

}
else
{
*(*(link[dest_link].cost_label + t) + r) = (int)(w->weight_link(dest_link, TIME_PERIOD - 1, N, &(w->BP[m][r]))*60+0.5);
}
}
}

for(r = 0; r < NUM_REAL; r++)
{
T=0;
myqueue SE;
indegree = N->InDegree(destIndex);
for(m = 0; m < indegree; m++)
{
dest_prev = N->LinkSource(destIndex, m);
dest_link= N->LinkIndex(dest_prev, destIndex);
if (N->Disabled(dest_link))
continue;
SE.enqueue(dest_link);
}
//SE.enqueue(destIndex);
int eventindex = *(*(realizationToEvent + TIME_PERIOD - 1) + r);
while(!SE.is_empty())
{
//jnode=SE.dequeue();
jlink=SE.dequeue();
jnode= N->LinkSource(jlink);
indegree = N->InDegree(jnode);
for(m = 0; m < indegree; m++)
{
inode = N->LinkSource(jnode, m);
link_curr = N->LinkIndex(inode, jnode);
if (N->Disabled(link_curr))
continue;
knode = *(*(link[jlink].next_node + TIME_PERIOD - 1) + eventindex);
tnode=N->LinkDest(jlink);
//link_next = N->LinkInIndex(jnode, knode);
if (inode==tnode)
{
T=INT_INFINITY;
}
//temp_cost_label = w->weight(inode, jnode, TIME_PERIOD - 1, r, N) + *(*(node[jnode].cost_label + TIME_PERIOD - 1) + eventindex);
if (w->det_flag[N->LinkIndex(inode, jnode)]==1)
{
//temp_cost_label = w->weight(inode, jnode, time, N, &(w->BP[m][0])) + *(*(node[jnode].cost_label + TIME_PERIOD - 1) + eventindex);
//temp_cost_label = w->weight(inode, jnode, time, N, &(w->BP[m][0])) + *(*(link[jlink].cost_label + TIME_PERIOD - 1) + eventindex) + T;
temp_cost_label = (int)(w->weight_link(link_curr, TIME_PERIOD - 1, N, &(w->BP[m][0]))*60+0.5) + *(*(link[jlink].cost_label + TIME_PERIOD - 1) + eventindex) + T;

}
else
{
//temp_cost_label = w->weight(inode, jnode, time, N, &(w->BP[m][r])) + *(*(node[jnode].cost_label + TIME_PERIOD - 1) + eventindex);
//temp_cost_label = w->weight(inode, jnode, time, N, &(w->BP[m][r])) + *(*(link[jlink].cost_label + TIME_PERIOD - 1) + eventindex) + T;
temp_cost_label = (int)(w->weight_link(link_curr, TIME_PERIOD - 1, N, &(w->BP[m][r]))*60+0.5) + *(*(link[jlink].cost_label + TIME_PERIOD - 1) + eventindex) + T;
}

//*(*(link[dest_link].cost_label + t) + r) = (int)(w->weight_link(dest_link, t, N, &(w->BP[m][0]))*60+0.5);
//*(*(link[dest_link].cost_label + t) + r) = (int)(w->weight_link(dest_link, t, N, &(w->BP[m][event[t][r]]))*60+0.5);

//if(temp_cost_label < *(*(node[inode].cost_label + TIME_PERIOD - 1) + eventindex))
if(temp_cost_label < *(*(link[link_curr].cost_label + TIME_PERIOD - 1) + eventindex))
{
//*(*(node[inode].cost_label + TIME_PERIOD - 1) + eventindex) = temp_cost_label;
//*(*(node[inode].next_node+TIME_PERIOD-1)+eventindex)=jnode;
*(*(link[link_curr].cost_label + TIME_PERIOD - 1) + eventindex) = temp_cost_label;
*(*(link[link_curr].next_node+TIME_PERIOD-1)+eventindex)=knode;
if (SE.have_this(link_curr) == 0)
SE.enqueue(link_curr);
}
}

}
SE.clear();
}
}
*/

/*
vector<long> foo(){
vector<long> v;
v.push_back(200);
return v;
}
*/
//void CDPI( Network *N, EdgeWeights *w, int dests, node_info *node, char *policyfile, int start_time, vector<pair<long, int>> orig_dep,int *policy) { 
void CDPI(Network *N, EdgeWeights *w, int dests, link_info *link, char *eventfile, char *policyfile, char *policyfile1, int start_time, vector<pair<long, int>> orig_dep, int *policy)
{
	double computation_time = 0.0;
	double computation_time1 = 0.0;
	//double computation_time2 = 0.0;
	double computation_time3 = 0.0;
	clock_t start = clock();

	// the complete dependence perfect online information variant
	int TIME_PERIOD, NUM_NODES, NUM_ARCS, NUM_REAL;
	int t, m, r, inode, jnode;
	double temp_cost_label;
	int tt;
	long i;
	//long outdegree;
	long linkindex;
	int n, knode;
	//int T=0;
	//int block_number=1440;
	int flag = 0;
	int flag2 = 0;
	int flag1 = 0;
	NUM_NODES = N->NumNodes();
	NUM_ARCS = N->NumLinks();
	NUM_REAL = w->numreal();
	TIME_PERIOD = w->levels();
	int period_length = w->agglength();
	vector<long> op_links; //the set of links on the optimal policy
	int **event;
	int **eventTag;
	int **realizationToEvent;
	event = new int*[TIME_PERIOD];
	eventTag = new int*[TIME_PERIOD]; //the starting position of a new event collection at time t
	realizationToEvent = new int*[TIME_PERIOD];
	for (t = 0; t < TIME_PERIOD; t++) {
		event[t] = new int[NUM_REAL];
		eventTag[t] = new int[NUM_REAL + 1]; //note the dimension is NUM_REAL+1
		realizationToEvent[t] = new int[NUM_REAL];
	}

	long ilink, jlink, klink;
	int ntime;
	//int time;
	int tnode;
	long indegree, link_curr, link_next, dest_prev, dest_link;

	clock_t start1 = clock();
	generate_event_collection(N, w, dests, orig_dep, event, eventTag, realizationToEvent, eventfile, policyfile, start_time);

	clock_t ends1 = clock();
	computation_time1 = (double)(ends1 - start1) / CLOCKS_PER_SEC;

	clock_t start3 = clock();
	//time = start_time + (TIME_PERIOD - 1);
	assert(dests >= 0 && dests < NUM_NODES);
	for (t = 0; t < TIME_PERIOD; t++)
	{
		for (r = 0; r < NUM_REAL; r++)
		{
			if (*(*(eventTag + t) + r) == 1)
			{
				for (n = 0; n < NUM_ARCS; n++)
				{
					*(*(link[n].cost_label + t) + r) = INT_INFINITY;
					*(*(link[n].next_node + t) + r) = -1;
				}
				indegree = N->InDegree(dests);
				for (m = 0; m < indegree; m++)
				{
					dest_prev = N->LinkSource(dests, m);
					dest_link = N->LinkIndex(dest_prev, dests);
					//if (N->Disabled(dest_link))
					//continue;
					*(*(link[dest_link].next_node + t) + r) = dests; //destination

					int i3, start3, end3;
					bp start_bp3;
					if (w->det_flag[dest_link] == 1)
					{
						if (w->blocked[dest_link][0][0] == 1)
						{
							continue;
						}
					}
					else
					{
						start_bp3 = w->find_interval(&(w->BP[dest_link][r]), start_time + t * period_length, &i3);
						start3 = i3 - 1;
						end3 = i3;
						if (w->blocked[dest_link][start3][r] == 1 || w->blocked[dest_link][end3][r] == 1)
						{
							continue;
						}
					}

					if (w->det_flag[dest_link] == 1)
					{
						//*(*(link[dest_link].cost_label + t) + r) = w->weight_link(dest_link, start_time + t * period_length, N, &(w->BP[dest_link][0]))*60.0;
						*(*(link[dest_link].cost_label + t) + r) = w->weight_link(dest_link, start_time + t * period_length, N, &(w->BP[dest_link][0]));

					}
					else
					{
						//*(*(link[dest_link].cost_label + t) + r) = w->weight_link(dest_link, start_time + t * period_length, N, &(w->BP[dest_link][r]))*60.0;
						*(*(link[dest_link].cost_label + t) + r) = w->weight_link(dest_link, start_time + t * period_length, N, &(w->BP[dest_link][event[t][r]]));

					}
				}
			}
		}
	}
	myqueue SE;
	indegree = N->InDegree(dests);
	for (m = 0; m < indegree; m++)
	{
		dest_prev = N->LinkSource(dests, m);
		dest_link = N->LinkIndex(dest_prev, dests);
		//if (N->Disabled(dest_link))
		//continue;
		SE.enqueue(dest_link);
	}
	while (!SE.is_empty())
	{
		jlink = SE.dequeue();
		jnode = N->LinkSource(jlink);
		indegree = N->InDegree(jnode);

		//bool dobreak1=false;
		//for(m = 0; m < indegree&&!dobreak1; m++)
		for (m = 0; m < indegree; m++)
		{
			flag = 0;
			inode = N->LinkSource(jnode, m);
			link_curr = N->LinkIndex(inode, jnode);
			//if (N->Disabled(link_curr))
			//continue;
			tnode = N->LinkDest(jlink);
			if (inode == tnode)
				continue;

			//bool dobreak2 = false;
			//for(t =0; t < TIME_PERIOD&&!dobreak2&&!dobreak1; t++)
			for (t = 0; t < TIME_PERIOD; t++)
			{
				int curr = 0; int next = 0;
				int rr, rrr;
				//for(r=1;r<NUM_REAL+1&&!dobreak2&&!dobreak1;r++)
				for (r = 1; r < NUM_REAL + 1; r++)
				{
					flag1 = 0;
					flag2 = 0;
					if (*(*(eventTag + t) + r) == 1)
					{
						temp_cost_label = 0;
						next = r;
						double probEV = 0.0;
						for (rr = curr; rr < next; rr++)
						{

							int	i2, start2, end2;
							bp start_bp2;
							//bp end_bp2;
							if (w->det_flag[link_curr] == 1)
							{
								if (w->blocked[link_curr][0][0] == 1)
								{
									//dobreak2=true;
									flag1 = 1;
									break;

								}
							}
							else
							{
								start_bp2 = w->find_interval(&(w->BP[link_curr][event[t][rr]]), start_time + t * period_length, &i2);
								start2 = i2 - 1;
								end2 = i2;
								//if (w->Blocked(link_curr,start2,r)==1||w->Blocked(link_curr,end2,r)==1)
								if (w->blocked[link_curr][start2][event[t][rr]] == 1 || w->blocked[link_curr][end2][event[t][rr]] == 1)
									//if (w->blocked[link_curr][start2][curr]==1||w->blocked[link_curr][end2][curr]==1)
								{
									//dobreak1=true;
									flag1 = 1;
									break;

								}
							}

							probEV += w->prob[*(*(event + t) + rr)];
						}
						if (flag1 == 1)
						{
							continue;
						}

						if (w->det_flag[link_curr] == 1)
						{
							//tt = (int)(t + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][0]))*60.0/period_length+0.5);
							tt = (int)(t + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][0])) / period_length + 0.5);

						}
						else
						{
							//tt = (int)(t + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][event[t][curr]]))*60.0/period_length+0.5);
							tt = (int)(t + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][event[t][curr]])) / period_length + 0.5);

						}
						//if (tt==t) tt=t+1; //make sure no zero travel time exists
						if (tt > TIME_PERIOD - 1) tt = TIME_PERIOD - 1;
						int curr1 = curr; int next1;
						for (int rr = curr + 1; rr <= next; rr++)
						{
							if (*(*(eventTag + tt) + rr) == 1)
							{
								next1 = rr;
								double probEVPrime = 0.0;
								for (rrr = curr1; rrr < next1; rrr++)
								{

									int i1, start1, end1;
									bp start_bp1;
									if (w->det_flag[jlink] == 1)
									{
										if (w->blocked[jlink][0][0] == 1)
										{
											flag2 = 1;
											break;
										}
									}
									else
									{
										//start_bp1 = w->find_interval(&(w->BP[jlink][event[tt][curr1]]), start_time + tt, &i1);
										start_bp1 = w->find_interval(&(w->BP[jlink][event[tt][rrr]]), start_time + tt * period_length, &i1);
										start1 = i1 - 1;
										end1 = i1;
										//if (w->blocked[jlink][start1][event[tt][curr1]]==1||w->blocked[jlink][end1][event[tt][curr1]]==1)
										if (w->blocked[jlink][start1][event[tt][rrr]] == 1 || w->blocked[jlink][end1][event[tt][rrr]] == 1)
											//if (w->blocked[jlink][start1][curr1]==1||w->blocked[jlink][end1][curr1]==1)
										{
											flag2 = 1;
											break;
										}
									}
									probEVPrime += w->prob[*(*(event + tt) + rrr)];
								}

								temp_cost_label += *(*(link[jlink].cost_label + tt) + curr1)*probEVPrime;
								curr1 = next1;
							}
						}
						if (flag2 == 1)
						{
							continue;
						}

						temp_cost_label /= probEV;
						if (w->det_flag[link_curr] == 1)
						{
							//temp_cost_label = temp_cost_label + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][0]))*60.0;
							temp_cost_label = temp_cost_label + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][0]));
						}
						else
						{
							//temp_cost_label = temp_cost_label + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][event[t][curr]]))*60.0;
							temp_cost_label = temp_cost_label + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][event[t][curr]]));
						}
						if (temp_cost_label < *(*(link[link_curr].cost_label + t) + curr))
							//if(temp_cost_label - *(*(link[link_curr].cost_label+t)+curr) < -EPSILON) 
						{
							*(*(link[link_curr].cost_label + t) + curr) = temp_cost_label;
							*(*(link[link_curr].next_node + t) + curr) = tnode;
							flag = 1;
							//if (SE.have_this(link_curr) == 0) 
							//SE.enqueue(link_curr);
						}

						curr = next;
						//N->SetDisabled(link_curr, 0);
					}
				}
			}
			if (SE.have_this(link_curr) == 0 && flag == 1)
				SE.enqueue(link_curr);
		}
	}
	SE.clear();

	clock_t ends3 = clock();
	computation_time3 = (double)(ends3 - start3) / CLOCKS_PER_SEC;
	clock_t ends = clock();
	computation_time = (double)(ends - start) / CLOCKS_PER_SEC;

	op_links = write_route_choice_result_new(N, w, dests, orig_dep, link, start_time, *policy, policyfile, policyfile1, 1, realizationToEvent, computation_time, computation_time1, computation_time3);
	(*policy)++;


	//generate choice set by link elimination
	//for (i = 0; i < op_links.size(); i+=5)
	//for (i = 0; i < op_links.size(); i+=10)

	for (i = 0; i < op_links.size(); i++)
	{
		N->SetDisabled(op_links[i], 1);

		generate_event_collection(N, w, dests, orig_dep, event, eventTag, realizationToEvent, eventfile, policyfile, start_time);
		assert(dests >= 0 && dests < NUM_NODES);
		for (t = 0; t < TIME_PERIOD; t++)
		{
			for (r = 0; r < NUM_REAL; r++)
			{
				if (*(*(eventTag + t) + r) == 1)
				{
					for (n = 0; n < NUM_ARCS; n++)
					{
						*(*(link[n].cost_label + t) + r) = INT_INFINITY;
						*(*(link[n].next_node + t) + r) = -1;
					}
					indegree = N->InDegree(dests);
					for (m = 0; m < indegree; m++)
					{
						dest_prev = N->LinkSource(dests, m);
						dest_link = N->LinkIndex(dest_prev, dests);
						if (N->Disabled(dest_link))
							continue;
						*(*(link[dest_link].next_node + t) + r) = dests; //destination

						int i3, start3, end3;
						bp start_bp3;
						if (w->det_flag[dest_link] == 1)
						{
							if (w->blocked[dest_link][0][0] == 1)
							{
								continue;
							}
						}
						else
						{
							start_bp3 = w->find_interval(&(w->BP[dest_link][r]), start_time + t * period_length, &i3);
							start3 = i3 - 1;
							end3 = i3;
							if (w->blocked[dest_link][start3][r] == 1 || w->blocked[dest_link][end3][r] == 1)
							{
								continue;
							}
						}

						if (w->det_flag[dest_link] == 1)
						{
							//*(*(link[dest_link].cost_label + t) + r) = w->weight_link(dest_link, start_time + t * period_length, N, &(w->BP[dest_link][0]))*60.0;
							*(*(link[dest_link].cost_label + t) + r) = w->weight_link(dest_link, start_time + t * period_length, N, &(w->BP[dest_link][0]));

						}
						else
						{
							//*(*(link[dest_link].cost_label + t) + r) = w->weight_link(dest_link, start_time + t * period_length, N, &(w->BP[dest_link][r]))*60.0;
							*(*(link[dest_link].cost_label + t) + r) = w->weight_link(dest_link, start_time + t * period_length, N, &(w->BP[dest_link][event[t][r]]));

						}
					}
				}
			}
		}
		myqueue SE;
		indegree = N->InDegree(dests);
		for (m = 0; m < indegree; m++)
		{
			dest_prev = N->LinkSource(dests, m);
			dest_link = N->LinkIndex(dest_prev, dests);
			if (N->Disabled(dest_link))
				continue;
			SE.enqueue(dest_link);
		}
		while (!SE.is_empty())
		{
			jlink = SE.dequeue();
			jnode = N->LinkSource(jlink);
			indegree = N->InDegree(jnode);

			//bool dobreak1=false;
			//for(m = 0; m < indegree&&!dobreak1; m++)
			for (m = 0; m < indegree; m++)
			{
				flag = 0;
				inode = N->LinkSource(jnode, m);
				link_curr = N->LinkIndex(inode, jnode);
				if (N->Disabled(link_curr))
					continue;
				tnode = N->LinkDest(jlink);
				if (inode == tnode)
					continue;

				//bool dobreak2 = false;
				//for(t =0; t < TIME_PERIOD&&!dobreak2&&!dobreak1; t++)
				for (t = 0; t < TIME_PERIOD; t++)
				{
					int curr = 0; int next = 0;
					int rr, rrr;
					//for(r=1;r<NUM_REAL+1&&!dobreak2&&!dobreak1;r++)
					for (r = 1; r < NUM_REAL + 1; r++)
					{
						flag1 = 0;
						flag2 = 0;
						if (*(*(eventTag + t) + r) == 1)
						{
							temp_cost_label = 0;
							next = r;
							double probEV = 0.0;
							for (rr = curr; rr < next; rr++)
							{

								int	i2, start2, end2;
								bp start_bp2;
								//bp end_bp2;
								if (w->det_flag[link_curr] == 1)
								{
									if (w->blocked[link_curr][0][0] == 1)
									{
										//dobreak2=true;
										flag1 = 1;
										break;

									}
								}
								else
								{
									start_bp2 = w->find_interval(&(w->BP[link_curr][event[t][rr]]), start_time + t * period_length, &i2);
									start2 = i2 - 1;
									end2 = i2;
									//if (w->Blocked(link_curr,start2,r)==1||w->Blocked(link_curr,end2,r)==1)
									if (w->blocked[link_curr][start2][event[t][rr]] == 1 || w->blocked[link_curr][end2][event[t][rr]] == 1)
										//if (w->blocked[link_curr][start2][curr]==1||w->blocked[link_curr][end2][curr]==1)
									{
										//dobreak1=true;
										flag1 = 1;
										break;

									}
								}

								probEV += w->prob[*(*(event + t) + rr)];
							}
							if (flag1 == 1)
							{
								continue;
							}

							if (w->det_flag[link_curr] == 1)
							{
								//tt = (int)(t + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][0]))*60.0/period_length+0.5);
								tt = (int)(t + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][0])) / period_length + 0.5);

							}
							else
							{
								//tt = (int)(t + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][event[t][curr]]))*60.0/period_length+0.5);
								tt = (int)(t + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][event[t][curr]])) / period_length + 0.5);

							}
							//if (tt==t) tt=t+1; //make sure no zero travel time exists
							if (tt > TIME_PERIOD - 1) tt = TIME_PERIOD - 1;
							int curr1 = curr; int next1;
							for (int rr = curr + 1; rr <= next; rr++)
							{
								if (*(*(eventTag + tt) + rr) == 1)
								{
									next1 = rr;
									double probEVPrime = 0.0;
									for (rrr = curr1; rrr < next1; rrr++)
									{

										int i1, start1, end1;
										bp start_bp1;
										if (w->det_flag[jlink] == 1)
										{
											if (w->blocked[jlink][0][0] == 1)
											{
												flag2 = 1;
												break;
											}
										}
										else
										{
											//start_bp1 = w->find_interval(&(w->BP[jlink][event[tt][curr1]]), start_time + tt, &i1);
											start_bp1 = w->find_interval(&(w->BP[jlink][event[tt][rrr]]), start_time + tt * period_length, &i1);
											start1 = i1 - 1;
											end1 = i1;
											//if (w->blocked[jlink][start1][event[tt][curr1]]==1||w->blocked[jlink][end1][event[tt][curr1]]==1)
											if (w->blocked[jlink][start1][event[tt][rrr]] == 1 || w->blocked[jlink][end1][event[tt][rrr]] == 1)
												//if (w->blocked[jlink][start1][curr1]==1||w->blocked[jlink][end1][curr1]==1)
											{
												flag2 = 1;
												break;
											}
										}
										probEVPrime += w->prob[*(*(event + tt) + rrr)];
									}

									temp_cost_label += *(*(link[jlink].cost_label + tt) + curr1)*probEVPrime;
									curr1 = next1;
								}
							}
							if (flag2 == 1)
							{
								continue;
							}

							temp_cost_label /= probEV;
							if (w->det_flag[link_curr] == 1)
							{
								//temp_cost_label = temp_cost_label + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][0]))*60.0;
								temp_cost_label = temp_cost_label + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][0]));
							}
							else
							{
								//temp_cost_label = temp_cost_label + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][event[t][curr]]))*60.0;
								temp_cost_label = temp_cost_label + w->weight_link(link_curr, start_time + t * period_length, N, &(w->BP[link_curr][event[t][curr]]));
							}
							if (temp_cost_label < *(*(link[link_curr].cost_label + t) + curr))
								//if(temp_cost_label - *(*(link[link_curr].cost_label+t)+curr) < -EPSILON) 
							{
								*(*(link[link_curr].cost_label + t) + curr) = temp_cost_label;
								*(*(link[link_curr].next_node + t) + curr) = tnode;
								flag = 1;
								//if (SE.have_this(link_curr) == 0) 
								//SE.enqueue(link_curr);
							}

							curr = next;
							//N->SetDisabled(link_curr, 0);
						}
					}
				}
				if (SE.have_this(link_curr) == 0 && flag == 1)
					SE.enqueue(link_curr);
			}
		}
		SE.clear();
		N->SetDisabled(op_links[i], 0);
		write_route_choice_result_new(N, w, dests, orig_dep, link, start_time, *policy, policyfile, policyfile1, 0, realizationToEvent, computation_time, computation_time1, computation_time3);
		(*policy)++;
	}



	//write_result(N, w, *it, eventTag, node, realizationToEvent, policyfile);
	//write_result_CDPI(N, w, dests, event, eventTag, node, realizationToEvent, policyfile);
	//free spaces
	for (t = 0; t < TIME_PERIOD; t++)
	{
		delete[] event[t];
		delete[] eventTag[t]; //note the dimension is NUM_REAL+1
		delete[] realizationToEvent[t];
	}
	delete[] event;
	delete[] eventTag;
	delete[] realizationToEvent;
}
/*
//calculate static shortest paths in the static period where link travel costs are expected values
void static_shortest_path_noi(Network *N, EdgeWeights *w, int destIndex, node_info *node, int start_time)
{

// Initialization and shortest paths calculation in static area for the noi heuristic
long 	TIME_PERIOD, NUM_NODES, NUM_ARCS, NUM_REAL, link, period_length;
int		i, n, m, t, r, time;
float	val;
long 	inode, jnode, indegree;
double 	temp_cost_label;
bp		next;

NUM_NODES = N->NumNodes();
TIME_PERIOD = w->levels();
period_length = w->agglength();
time = start_time + (TIME_PERIOD - 1) * period_length;

assert(destIndex >= 0 && destIndex < NUM_NODES);
//assert( w->minweight() > 0.5 );


for(t = 0; t < TIME_PERIOD; t++)
{
for(n = 0; n < NUM_NODES; n++)
{
node[n].cost_label_noi[t] = INT_INFINITY;
node[n].next_node_noi[t] = -1;
}
node[destIndex].cost_label_noi[t] = 0;
}


myqueue SE;
SE.enqueue(destIndex);

while(!SE.is_empty())
{
jnode = SE.dequeue();

indegree = N->InDegree(jnode);
for(m = 0; m < indegree; m++)
{
inode = N->LinkSource(jnode, m);
link = N->LinkInIndex(jnode, m);
if (N->Disabled(link))
continue;

if (w->det_flag[link])
val = w->weight_link(link, time, N, &(w->BP[link][0])) * 60;
else
val = w->weight_link(link, time, N, &(w->ExpectedBP[link])) * 60;
temp_cost_label = val + node[jnode].cost_label_noi[TIME_PERIOD - 1];
if (temp_cost_label < node[inode].cost_label_noi[TIME_PERIOD - 1])
{
node[inode].cost_label_noi[TIME_PERIOD - 1] = temp_cost_label;
node[inode].next_node_noi[TIME_PERIOD - 1] = jnode;
if (!SE.have_this(inode))
SE.enqueue(inode);
}
}

}
}
*/
/*
void static_shortest_path_noi_StepTT(Network *N, EdgeWeights *w, int destIndex, node_info *node, int start_time)
{

// Initialization and shortest paths calculation in static area for the noi heuristic
long 	TIME_PERIOD, NUM_NODES, NUM_ARCS, NUM_REAL, link, period_length;
int		i, n, m, t, r, time;
float	val;
long 	inode, jnode, indegree, timeslice, p;
double 	temp_cost_label;

NUM_NODES = N->NumNodes();
TIME_PERIOD = w->levels();
period_length = w->agglength();

assert(destIndex >= 0 && destIndex < NUM_NODES);
//assert( w->minweight() > 0.5 );


for(t = 0; t < TIME_PERIOD; t++)
{
for(n = 0; n < NUM_NODES; n++)
{
node[n].cost_label_noi[t] = DOUBLE_MAX;
node[n].next_node_noi[t] = -1;
}
if (disutility)
node[destIndex].cost_label_noi[t] = exp((double)-startingmean);
else
node[destIndex].cost_label_noi[t] = 0;
}


myqueue SE;
SE.enqueue(destIndex);

timeslice = w->find_interval(0, start_time + (TIME_PERIOD - 1) * period_length);
while(!SE.is_empty())
{
jnode = SE.dequeue();

indegree = N->InDegree(jnode);
for(m = 0; m < indegree; m++)
{
inode = N->LinkSource(jnode, m);
link = N->LinkInIndex(jnode, m);
if (N->Disabled(link))
continue;

val = .0;
for (p = 0; p < w->numsupportpoints[link][timeslice]; p++)
val += w->supportpoints[link][timeslice][p] * w->probability[link][timeslice][p];

if (disutility)
{
temp_cost_label = exp((double)(alpha * val - startingmean)) * node[jnode].cost_label_noi[TIME_PERIOD - 1];
if (exp((double)(log(temp_cost_label) + startingmean)) < node[inode].cost_label_noi[TIME_PERIOD - 1])
{
node[inode].cost_label_noi[TIME_PERIOD - 1] = exp((double)(log(temp_cost_label) + startingmean));
node[inode].next_node_noi[TIME_PERIOD - 1] = jnode;
if (!SE.have_this(inode))
SE.enqueue(inode);
}
}
else
{
temp_cost_label = val + node[jnode].cost_label_noi[TIME_PERIOD - 1];
if (temp_cost_label < node[inode].cost_label_noi[TIME_PERIOD - 1])
{
node[inode].cost_label_noi[TIME_PERIOD - 1] = temp_cost_label;
node[inode].next_node_noi[TIME_PERIOD - 1] = jnode;
if (!SE.have_this(inode))
SE.enqueue(inode);
}
}
}

}
}
*/
/*
void CE( Network *N, EdgeWeights *w, vector<int> dests, node_info *node, char *policyfile, int start_time, int flag) {
//flag indicated whether this is a standalone CE or a component of OLFCE
//flag = 0: CE flag = 1: OLFCE

register int TIME_PERIOD, NUM_NODES, NUM_ARCS, NUM_REAL;
register int t, m, r, inode, jnode, tt;
float temp_cost_label;

NUM_NODES = N->NumNodes();
NUM_ARCS = N->NumLinks();
NUM_REAL = w->numreal();
TIME_PERIOD = w->levels();

//assert( w->minweight() > 0.5 );

if (flag == 0)
w->calexpected(N); //calculate the expected time-dependent link travel times
//if (flag == 1) expected tt has been calculated in OLFCE

for (vector<int>::iterator it = dests.begin(); it != dests.end(); it++)
{
static_shortest_path_ce(N, w, *it, node);

//begin the main loop
for (t = TIME_PERIOD - 2; t >= start_time; t--)
for(m = 0; m < NUM_ARCS; m++)
{
inode=N->LinkSource(m);
jnode=N->LinkDest(m);

tt = (int)(t + w->expectedtt[m][t] + 0.5);
if (tt == t)
tt = t + 1; //make sure no zero travel time exists
if( tt > TIME_PERIOD - 1)
tt = TIME_PERIOD-1;

temp_cost_label = w->expectedtt[m][t] + *(node[jnode].cost_label_ce + tt);

if(temp_cost_label < *(node[inode].cost_label_ce + t))
{
*(node[inode].cost_label_ce + t) = temp_cost_label;
*(node[inode].next_node_ce + t) = jnode;
}
}

if (flag == 0)
write_result_simple(N, w, *it,  node, policyfile);
}
}
*/
/*
void NOI(Network *N, EdgeWeights *w, odpairs ods, node_info *node, char *policyfile, int start_time, int flag, int *policy)
{
//flag indicated whether this is a stand-alone NOI or a component of OLFNOI
//flag = 0: NOI flag = 1: open-loop-feedback with NOI

long TIME_PERIOD, NUM_NODES, NUM_ARCS, NUM_REAL;
long i, t, m, r, inode, jnode, tt, p, n, period_length;
float temp_cost_label, weight;
char    x[256];
vector<long> op_links; //the set of links on the optimal policy

NUM_NODES = N->NumNodes();
NUM_ARCS = N->NumLinks();
NUM_REAL = w->numreal();
TIME_PERIOD = w->levels();
period_length = w->agglength();

if (flag == 0)
w->calexpectedBP(N);

for (odpairs::iterator it = ods.begin(); it != ods.end(); it++)
{
static_shortest_path_noi(N, w, (*it).first, node, start_time);

//begin the main loop
for (t = TIME_PERIOD - 2; t >= 0; t--)
for (m = 0; m < NUM_ARCS; m++)
{
if (N->Disabled(m))
continue;

inode = N->LinkSource(m);
if (inode < 0 || inode >= NUM_NODES)
cin >> x;
jnode = N->LinkDest(m);
if (jnode < 0 || jnode >= NUM_NODES)
cin >> x;

if (w->det_flag[m]==1)
{
weight = w->weight_link(m, t * period_length, N, &(w->BP[m][0])) * 60;
tt = (int)(t + weight / period_length + 0.5); //travel time is in minutes, while period length is in seconds
if (tt == t)
tt = t + 1; //make sure arrival at least next interval
if (tt > TIME_PERIOD - 1)
tt = TIME_PERIOD - 1;

temp_cost_label = weight + node[jnode].cost_label_noi[tt];
}
else
{
temp_cost_label = .0;
for (p = 0; p < w->numreal(); p++)
{
weight = w->weight_link(m, t * period_length, N, &(w->BP[m][p])) * 60;
tt = (int)(t + weight / period_length + 0.5);
if (tt == t)
tt = t + 1; //make sure arrival at least next interval
if (tt > TIME_PERIOD - 1)
tt = TIME_PERIOD - 1;

temp_cost_label += (weight + node[jnode].cost_label_noi[tt]) * w->prob[p];
}
}

if (temp_cost_label < node[inode].cost_label_noi[t])
{
node[inode].cost_label_noi[t] = temp_cost_label;
node[inode].next_node_noi[t] = jnode;
}
}

op_links = write_route_choice_result(N, w, (*it).first, (*it).second, node, start_time, *policy, policyfile, 1);
(*policy)++;

//generate choice set by link elimination
for (i = 0; i < op_links.size(); i++)
{
N->SetDisabled(op_links[i], 1);

static_shortest_path_noi(N, w, (*it).first, node, start_time);

//begin the main loop
for (t = TIME_PERIOD - 2; t >= 0; t--)
for (m = 0; m < NUM_ARCS; m++)
{
if (N->Disabled(m))
continue;

inode = N->LinkSource(m);
if (inode < 0 || inode >= NUM_NODES)
cin >> x;
jnode = N->LinkDest(m);
if (jnode < 0 || jnode >= NUM_NODES)
cin >> x;

if (w->det_flag[m]==1)
{
weight = w->weight_link(m, t * period_length, N, &(w->BP[m][0])) * 60;
tt = (int)(t + weight / period_length + 0.5); //travel time is in minutes, while period length is in seconds
if (tt == t)
tt = t + 1; //make sure arrival at least next interval
if (tt > TIME_PERIOD - 1)
tt = TIME_PERIOD - 1;

temp_cost_label = weight + node[jnode].cost_label_noi[tt];
}
else
{
temp_cost_label = .0;
for (p = 0; p < w->numreal(); p++)
{
weight = w->weight_link(m, t * period_length, N, &(w->BP[m][p])) * 60;
tt = (int)(t + weight / period_length + 0.5);
if (tt == t)
tt = t + 1; //make sure arrival at least next interval
if (tt > TIME_PERIOD - 1)
tt = TIME_PERIOD - 1;

temp_cost_label += (weight + node[jnode].cost_label_noi[tt]) * w->prob[p];
}
}

if (temp_cost_label < node[inode].cost_label_noi[t])
{
node[inode].cost_label_noi[t] = temp_cost_label;
node[inode].next_node_noi[t] = jnode;
}
}
N->SetDisabled(op_links[i], 0);
write_route_choice_result(N, w, (*it).first, (*it).second, node, start_time, *policy, policyfile, 0);
(*policy)++;
}
}
}
*/
/*
void remove_cycles(Network *N, EdgeWeights *w, long destIdx, vector<pair<long, int>> orig_dep, node_info *node, int start_time)
{
//for each entry in the queue, store the sequence of (node, time) and the probability
//check whether the new next node is already in the sequence, if so, cycle
//to remove cycle, go backwards from current node, and divert at the first possible node

vector<pair<vector<pair<long, long>>, float>>	reachables_queue;
//a vector of realized paths so far, where each path is represented by a sequence of (node, time) pairs and
//the associated probability
//this is an exact method to obtain the realized paths and probabilities
pair<long, long>			node_time;
vector<pair<long, long>>	path;
float						prob;
pair<vector<pair<long, long>>, float>	path_prob;

long	i, j, orig, nodes, currnode, nextnode, jnode, currtime, nexttime, timeslice, curr, next, link,
t1, t2, TIME_PERIOD = w->levels(), m;
int		dep, period_length = w->agglength(), flag, outdegree, p, tt;
char	x[256];
double  weight, temp_cost_label;

nodes = N->NumNodes();

for (i = 0; i < orig_dep.size(); i++)
{

//breadth search to find all reacheable (node, time) pairs
t1 = (long)((orig_dep[i].second - start_time) / period_length + 0.5);
node_time.first = orig_dep[i].first;
node_time.second = t1;
path.clear();
path.push_back(node_time);
path_prob.first = path;
path_prob.second = 1;
reachables_queue.clear();
reachables_queue.push_back(path_prob);  //put the origin in the queue

while (!reachables_queue.empty())
{
//pop out a path-prob pair
path_prob = reachables_queue.front();
reachables_queue.erase(reachables_queue.begin());
path = path_prob.first;
prob = path_prob.second;
node_time = path.front();
currnode = node_time.first;
if (currnode == destIdx)
continue;
t1 = node_time.second;
timeslice = w->find_interval(0, start_time + t1 * period_length);

//add outgoing (node,time) pairs from the current node-time pair
nextnode = node[currnode].next_node_noi[t1];

//check whether a cycle is formed

link = N->LinkIndex(currnode, nextnode);

node_time.first = nextnode;
for (j = 0; j < w->numsupportpoints[link][timeslice]; j++)
{
if (!(w->probability[link][timeslice][j] > 0))
continue;
weight = w->supportpoints[link][timeslice][j];
t2 = (long)((t1 * period_length + weight)/ period_length + 0.5);
if (t2 == t1)
t2++;
if (t2 > TIME_PERIOD - 1)
t2 = TIME_PERIOD - 1;
node_time.second = t2;
//check whether the element already exists in the queue
flag = 0;
for (vector<pair<long,long>>::iterator x = reachables_queue.begin(); x != reachables_queue.end(); x++)
if (*x == node_time)
{
flag = 1;
break;
}
if (!flag)
reachables_queue.push_back(node_time);
flag = 0;
}
}
}
}

*/


/*
void NOI_StepTT(Network *N, EdgeWeights *w, odpairs ods, node_info *node, char *policyfile, int start_time, int flag, int *policy)
{
//flag indicated whether this is a stand-alone NOI or a component of OLFNOI
//flag = 0: NOI flag = 1: open-loop-feedback with NOI

long TIME_PERIOD, NUM_NODES, NUM_ARCS, NUM_REAL;
long i, t, m, r, inode, jnode, tt, p, n, period_length, timeslice;
double temp_cost_label, weight;
char    x[256];
vector<long> op_links; //the set of links on the optimal policy

NUM_NODES = N->NumNodes();
NUM_ARCS = N->NumLinks();
NUM_REAL = w->numreal();
TIME_PERIOD = w->levels();
period_length = w->agglength();


for (odpairs::iterator it = ods.begin(); it != ods.end(); it++)
{
static_shortest_path_noi_StepTT(N, w, (*it).first, node, start_time);

//begin the main loop
for (t = TIME_PERIOD - 2; t >= 0; t--)
for (m = 0; m < NUM_ARCS; m++)
{
if (N->Disabled(m))
continue;

inode = N->LinkSource(m);
if (inode == (*it).first)
continue;
if (inode < 0 || inode >= NUM_NODES)
cin >> x;
jnode = N->LinkDest(m);
if (jnode < 0 || jnode >= NUM_NODES)
cin >> x;

timeslice = w->find_interval(0, start_time + t * period_length);
temp_cost_label = 0.;
for (p = 0; p < w->numsupportpoints[m][timeslice]; p++)
{
weight = w->supportpoints[m][timeslice][p];
if (w->probability[m][timeslice][p] > 0)
{
tt = (int)((t * period_length + weight) / period_length + 0.5);
if (tt == t)
tt = t + 1; //make sure arrival at least next interval
if (tt > TIME_PERIOD - 1)
tt = TIME_PERIOD - 1;

if (disutility)
temp_cost_label += (exp((double)(alpha * weight - startingmean)) * node[jnode].cost_label_noi[tt]) * w->probability[m][timeslice][p];
else
temp_cost_label += (weight + node[jnode].cost_label_noi[tt]) * w->probability[m][timeslice][p];
}
}

if (disutility)
{
if (exp((double)(log(temp_cost_label) + startingmean)) < node[inode].cost_label_noi[t])
{
node[inode].cost_label_noi[t] = exp((double)(log(temp_cost_label) + startingmean));
node[inode].next_node_noi[t] = jnode;
}
}
else
{
if (temp_cost_label < node[inode].cost_label_noi[t])
{
node[inode].cost_label_noi[t] = temp_cost_label;
node[inode].next_node_noi[t] = jnode;
}
}
}

//remove_cycles(N, w, (*it).first, (*it).second, node, start_time);

op_links = Policy2Path_NOI(N, w, (*it).first, (*it).second, node, start_time, *policy, policyfile, 1);
(*policy)++;

//op_links.clear();

//generate choice set by link elimination
for (i = 0; i < op_links.size(); i++)
{
N->SetDisabled(op_links[i], 1);

static_shortest_path_noi_StepTT(N, w, (*it).first, node, start_time);

//begin the main loop
for (t = TIME_PERIOD - 2; t >= 0; t--)
for (m = 0; m < NUM_ARCS; m++)
{
if (N->Disabled(m))
continue;

inode = N->LinkSource(m);
if (inode == (*it).first)
continue;
if (inode < 0 || inode >= NUM_NODES)
cin >> x;
jnode = N->LinkDest(m);
if (jnode < 0 || jnode >= NUM_NODES)
cin >> x;


timeslice = w->find_interval(0, start_time + t * period_length);
temp_cost_label = .0;
for (p = 0; p < w->numsupportpoints[m][timeslice]; p++)
{
weight = w->supportpoints[m][timeslice][p];
if (w->probability[m][timeslice][p] > 0)
{
tt = (int)((t * period_length + weight + 0.5) / period_length);
if (tt == t)
tt = t + 1; //make sure arrival at least next interval
if (tt > TIME_PERIOD - 1)
tt = TIME_PERIOD - 1;

if (disutility)
temp_cost_label += (exp((double)(alpha * weight - startingmean)) * node[jnode].cost_label_noi[tt]) * w->probability[m][timeslice][p];
else
temp_cost_label += (weight + node[jnode].cost_label_noi[tt]) * w->probability[m][timeslice][p];
}
}

if (disutility)
{
if (exp((double)(log(temp_cost_label) + startingmean)) < node[inode].cost_label_noi[t])
{
node[inode].cost_label_noi[t] = exp((double)(log(temp_cost_label) + startingmean));
node[inode].next_node_noi[t] = jnode;
}
}
else
{
if (temp_cost_label < node[inode].cost_label_noi[t])
{
node[inode].cost_label_noi[t] = temp_cost_label;
node[inode].next_node_noi[t] = jnode;
}
}
}
N->SetDisabled(op_links[i], 0);
//remove_cycles(N, w, (*it).first, (*it).second, node, start_time);
Policy2Path_NOI(N, w, (*it).first, (*it).second, node, start_time, *policy, policyfile, 0);
(*policy)++;
}
}
}
*/

/*
/*
void OLFCE( Network *N, EdgeWeights *w, vector<int> dests, node_info *node, char *policyfile)
{
register int TIME_PERIOD, NUM_NODES, NUM_ARCS, NUM_REAL;
register int t, m, r, inode, jnode, curr, next, n;
double temp_cost_label;

NUM_NODES = N->NumNodes();
NUM_ARCS = N->NumLinks();
NUM_REAL = w->numreal();
TIME_PERIOD = w->levels();

//assert( w->minweight() > 0.5 );

int **event;
int **eventTag;
int **realizationToEvent;

event = new int*[TIME_PERIOD];
eventTag = new int*[TIME_PERIOD]; //the starting position of a new event collection at time t
realizationToEvent=new int*[TIME_PERIOD];

for(t=0;t<TIME_PERIOD;t++)
{
event[t]=new int[NUM_REAL];
eventTag[t]=new int[NUM_REAL+1]; //note the dimension is NUM_REAL+1
realizationToEvent[t]=new int[NUM_REAL];
}

generate_event_collection(N, w, event, eventTag, realizationToEvent, policyfile);

for (vector<int>::iterator it = dests.begin(); it != dests.end(); it++)
{

static_shortest_path(N, w, *it,  eventTag, realizationToEvent, node);

vector<int> dest;
dest.push_back(*it);

//begin the main loop

for(t = 0;t < TIME_PERIOD - 1; t++)
{
curr=0;
next=0;
for(r = 1;r < NUM_REAL + 1; r++)
{
if(*(*(eventTag+t)+r)==1)
{
next=r;
w->calexpected(N, curr, next, event, t);
CE(N, w, dest, node, policyfile, t, 1);

for(n=0;n<NUM_NODES;n++)
{
*(*(node[n].cost_label+t)+curr)=node[n].cost_label_ce[t];
*(*(node[n].next_node+t)+curr)=node[n].next_node_ce[t];
}
curr=next;
}
}//for r
}//for t

write_result(N, w, *it, eventTag, node, realizationToEvent, policyfile);
dest.erase(dest.begin(), dest.end());
}//for it

//free spaces
for(t=0;t<TIME_PERIOD;t++)
{
delete[] event[t];
delete[] eventTag[t]; //note the dimension is NUM_REAL+1
delete[] realizationToEvent[t];
}
delete[] event;
delete[] eventTag;
delete[] realizationToEvent;
}

void OLF(Network *N, EdgeWeights *w, vector<int> dests, node_info *node, char *policyfile)
{

register int TIME_PERIOD, NUM_NODES, NUM_ARCS, NUM_REAL;
register int t, m, r, inode, jnode;
double temp_cost_label;

NUM_NODES = N->NumNodes();
NUM_ARCS = N->NumLinks();
NUM_REAL = w->numreal();
TIME_PERIOD = w->levels();

//assert( w->minweight() > 0.5 );

int **event;
int **eventTag;
int **realizationToEvent;

event = new int*[TIME_PERIOD];
eventTag = new int*[TIME_PERIOD]; //the starting position of a new event collection at time t
realizationToEvent=new int*[TIME_PERIOD];

for(t=0;t<TIME_PERIOD;t++) {
event[t]=new int[NUM_REAL];
eventTag[t]=new int[NUM_REAL+1]; //note the dimension is NUM_REAL+1
realizationToEvent[t]=new int[NUM_REAL];
}

generate_event_collection(N, w, event, eventTag, realizationToEvent, policyfile);

for (vector<int>::iterator it = dests.begin(); it != dests.end(); it++) {

static_shortest_path(N, w, *it,  eventTag, realizationToEvent, node);

//begin the main loop

for(t = 0; t < TIME_PERIOD - 1; t++) {
int curr=0; int next=0;

for(r=1;r<NUM_REAL+1;r++)
if(*(*(eventTag+t)+r)==1){
next=r;

//find the minimum expected travel time path conditional on the current realizations
float minpathtime, currpathtime;
int minpathindex;

for(int n=0;n<NUM_NODES;n++) {
if (n!=*it) {
minpathtime = INT_INFINITY;
for (int i=0;i<N->numpaths(n, *it);i++) {
currpathtime = (w->pathtimes(n, *it, i, t, N, curr, next, 1, event))[0];
if (currpathtime < minpathtime) {
minpathtime = currpathtime;
minpathindex = i;
}
}


//*(*(node[n].cost_label+t)+curr)=node[n].cost_label_ce[t];
*(*(node[n].next_node+t)+curr)=N->LinkDest((N->path2links(n, *it, minpathindex))[0]);
}
}
curr=next;
}
}

write_result(N, w, *it, eventTag, node, realizationToEvent, policyfile);

}
//free spaces
for(t=0;t<TIME_PERIOD;t++)
{
delete[] event[t];
delete[] eventTag[t]; //note the dimension is NUM_REAL+1
delete[] realizationToEvent[t];
}
delete[] event;
delete[] eventTag;
delete[] realizationToEvent;
}


//It does not make any sense to have the event alone.  Eithere we have both the events and the event tags.  Or we have none.
*/
void write_event(Network *N, EdgeWeights *w, long destIdx, vector<pair<long, int>> orig_dep, int **eventtag, char *filename) {

FILE *fp;
register int TIME_PERIOD, NUM_REAL;

fp = fopen( filename, "a");
if ( !fp ) {
printf( "Error: can't open %s\n", filename );
exit(1);
}

NUM_REAL = w->numreal();
TIME_PERIOD = w->levels();

//assert( w->minweight() > 0 );
long orig, destID;
int dep;
destID = N->NodeId(destIdx);


for (int i = 0; i < orig_dep.size(); i++)
{
	orig = orig_dep[i].first;
	dep = orig_dep[i].second;

		
	//write event collection
	for (int i = 0; i < TIME_PERIOD; i++)
	{
		fprintf(fp, "%16d%16d%16d ", N->NodeId(orig), destID, dep);
		for (int j = 0; j < NUM_REAL; j++)
			
			fprintf(fp, "%d ", eventtag[i][j]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	fclose(fp);

}

}

/*
void write_result_CDPI(Network *N, EdgeWeights *w, int destIndex, int **event,
int **eventTag, node_info *node, int **realizationToEvent, char *filename)
{
register int TIME_PERIOD, NUM_NODES, NUM_ARCS, NUM_REAL, agglength;

NUM_NODES = N->NumNodes();
NUM_ARCS = N->NumLinks();
NUM_REAL = w->numreal();
TIME_PERIOD = w->levels();

agglength = w->agglength();

int n,t,r,currArc,currNode,currTime,cost;
double **aggregateCost=new double*[NUM_NODES]; //expected cost from each node at each departure time
double aggregateProb; //probability of each event collection
FILE *fpDOT;
float	sum = 0.;

float *prob = w->prob;

fpDOT = fopen(filename, "a");
if (!fpDOT)
{
printf("Cannot open file %s!", filename);
exit(1);
}


for (n = 0; n < NUM_NODES; n++)
{
aggregateCost[n] = new double[TIME_PERIOD];
for (t = 0; t < TIME_PERIOD; t++)
aggregateCost[n][t] = 0;
}

fprintf(fpDOT,"*** RESULTS FROM DOT-SPI ALGORITHM ***\n");
fprintf(fpDOT,"NUM_NODES: %d NUM_ARCS: %d TIME_PERIOD: %d NUM_REAL: %d\n", NUM_NODES,NUM_ARCS,TIME_PERIOD,NUM_REAL);

for(t = 0; t < TIME_PERIOD;t++)
{
fprintf(fpDOT,"\n*** At Time %3d ***: ",t);

int curr=0;
int next=0;
for(int r=1;r<NUM_REAL+1;r++)
if(*(*(eventTag+t)+r)==1)
{
next=r;
//fprintf(fpDOT, "%\\psframe[dimen=middle](%d,%d)(%d,%d)\n",curr,t,next,t+1);

aggregateProb=0.0;

//fprintf(fpDOT,"\n*** For Event Collection (");
for(int rr=curr;rr<next;rr++)
{
//fprintf(fpDOT,"%d ",*(*(event+t)+rr));
aggregateProb+=prob[*(*(event+t)+rr)];
}
//fprintf(fpDOT,") ***\n");

for(n=0;n<NUM_NODES;n++)
{
//fprintf(fpDOT,"\nFrom Node %d: ",n);
//fprintf(fpDOT,"   Cost = %.5f ",*(*(node[n].cost_label+t)+curr));
//fprintf(fpDOT,"   Next_Node = %3d",*(*(node[n].next_node+t)+curr));

aggregateCost[n][t]+=*(*(node[n].cost_label+t)+curr)*aggregateProb;
}
curr=next;
//fprintf(fpDOT,"\n");
}

fprintf(fpDOT,"\n*** Aggregated Expected Cost at Time %3d ***\n",t);
for(n=0;n<NUM_NODES;n++)
{
//fprintf(fpDOT,"\nFrom Node %d: ",n);
//fprintf(fpDOT,"   Cost = %.5f ",aggregateCost[n][t]);
if (n == 0 && t >= 30 && t < 60) //from origin only
sum += aggregateCost[n][t] * agglength;
//fprintf(fpDOT,"%.5f\n",aggregateCost[n][t] * agglength);
}
//fprintf(fpDOT,"\n\n");
}
fprintf(fpDOT, "\n");
fprintf(fpDOT, "%.5f\n", sum/30.0);


//the following is to get the expected tt from joint realizations
int **next_arc=new int*[NUM_NODES*TIME_PERIOD];

for(t=0;t<TIME_PERIOD;t++)
for(n=0;n<NUM_NODES-1;n++)
{
//cout<<"t: "<<t<<" n: "<<n<<endl;
next_arc[t*NUM_NODES+n]=new int[NUM_REAL];

int curr=0;int next=0;
for(r=1;r<NUM_REAL+1;r++)
if(eventTag[t][r]==1)
{
next=r;

int nextNode=*(*(node[n].next_node+t)+curr);
//cout<<"curr: "<<curr<<" event collection starting with: "<<*(*(event+t)+curr)<<endl;
//cout<<"nextNode: "<<nextNode<<endl;
for(int m=node[nextNode].arc_start;m<node[nextNode+1].arc_start;m++)
if(arc[m].inode==n)
*(*(next_arc+t*NUM_NODES+n)+curr)=m;

curr=next;
}
}

for(t = 0;t < TIME_PERIOD; t++)
for (n = 0; n < NUM_NODES; n++)
node[n].true_cost_label[t] = 0;

for(r=0;r<NUM_REAL;r++)
for(n=0;n<NUM_NODES-1;n++)
for(t=0;t<TIME_PERIOD;t++)
{
currNode=n;currTime=t;
do
{
currArc=*(*(next_arc+currTime*NUM_NODES+currNode)+*(*(realizationToEvent+currTime)+r));
cost=*(*(realization+currTime*NUM_ARCS+currArc)+r);
node[n].true_cost_label[t] += cost * prob[r];
currNode = *(*(node[currNode].next_node+currTime)+*(*(realizationToEvent+currTime)+r));
currTime += cost;
if(currTime>TIME_PERIOD-1)
currTime=TIME_PERIOD-1;
}while(currNode!=NUM_NODES-1); //the destination node
}

for(t=0;t<TIME_PERIOD;t++)
{
for(n=0;n<NUM_NODES;n++)
{
fprintf(fpDOT,"\nFrom Node %d: ",n);
fprintf(fpDOT,"   Claimed Cost = %.5f  True Cost = %.5f",aggregateCost[n][t], node[n].true_cost_label[t]);
//fprintf(fpDOT,"%.5f ",node[n].true_cost_label[t]);
}
fprintf(fpDOT,"\n\n");
}

fclose(fpDOT);
delete[] aggregateCost;
}
*/
void write_result(Network *N, EdgeWeights *w, int destIndex, int **eventTag, node_info *node, int **realizationToEvent, char *filename) {

	FILE *fp;
	register int TIME_PERIOD, NUM_NODES, NUM_ARCS, NUM_REAL, event, n, i, j;

	fp = fopen(filename, "a");
	if (!fp)
	{
		printf("Error: can't open %s\n", filename);
		exit(1);
	}

	NUM_NODES = N->NumNodes();
	NUM_ARCS = N->NumLinks();
	NUM_REAL = w->numreal();
	TIME_PERIOD = w->levels();

	assert(destIndex >= 0 && destIndex < NUM_NODES);
	//assert( w->minweight() > 0 );

	int nextnode;

	//write the destination node
	fprintf(fp, "%d\n", N->NodeId(destIndex));

	//write next node for each node at each time period for each realization (r)
	//note the next nodes are sorted by realization index, not by event collection.  
	for (n = 0; n < NUM_NODES; n++)
	{
		for (i = 0; i < TIME_PERIOD; i++)
		{
			for (j = 0; j < NUM_REAL; j++)
			{
				event = realizationToEvent[i][j];
				nextnode = (node[n].next_node)[i][event];
				if (nextnode >= 0)
					fprintf(fp, "%d ", N->NodeId(nextnode)); //next node is in node index, while link travel times are in order of link ID
				//fprintf(fp,"%f ", (node[n].cost_label)[i][event]*w->agglength()); //next node is in node index, while output link travel times are in order of link ID
				else
					fprintf(fp, "%d ", nextnode);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}


vector<long> write_route_choice_result(Network *N, EdgeWeights *w, long destIdx, vector<pair<long, int>> orig_dep, node_info *node, int start_time, int policy, char *filename, int op_flag)
{
	//calculate the true expected travel time and obtain corresponding paths (for now, only set of links)
	long	orig, nodes, destID, curr, next, link, TIME_PERIOD = w->levels();
	int		dep, period_length = w->agglength(), i, r, t1, t2, *sorted, opsize;
	float   exp_weight, sum_weight, weight, time;
	double *links;
	FILE	*fp;
	vector<long>	op_links;

	nodes = N->NumNodes();
	destID = N->NodeId(destIdx);

	if (policy == 0)
		fp = fopen(filename, "w");
	else
		fp = fopen(filename, "a");

	if (!fp)
	{
		printf("Error: cannot open file %s.\n", filename);
		exit(1);
	}

	//sg	if (policy == 0)
	//sg		fprintf(fp, "%4s%20s%20s%20s%20s%20s\n", "//", "Origin", "Destination", "Departure", "Policy#", "Mean");

	for (i = 0; i < orig_dep.size(); i++)
	{
		orig = orig_dep[i].first;
		dep = orig_dep[i].second;

		exp_weight = 0;
		for (r = 0; r < w->numreal(); r++)
		{
			fprintf(fp, "%10d%10d%10d", N->NodeId(orig), destID, dep);

			curr = orig;
			time = dep;
			t2 = (int)((time - start_time) / period_length + 0.5);
			//if (t2 == t1)
			//t2++;
			if (t2 > TIME_PERIOD - 1)
				t2 = TIME_PERIOD - 1;
			next = node[curr].next_node_noi[t2];
			t1 = t2;

			sum_weight = 0;
			while (next != -1)
			{
				link = N->LinkIndex(curr, next);
				if (link == -1)
				{
					printf("Wrong routing policy (node %d, time period %d : next node %d) to destination %d.\n",
						N->NodeId(curr), t2, N->NodeId(next), destID);
					exit(1);
				}
				fprintf(fp, "%10d", N->LinkId(link));
				if (op_flag)
					op_links.push_back(link);

				if (w->det_flag[link])
					//weight = w->weight_link(link, time, N, &(w->BP[link][0])) * 60.0;	  //travel time is in minutes
					weight = w->weight_link(link, time, N, &(w->BP[link][0]));
				else
					//weight = w->weight_link(link, time, N, &(w->BP[link][r])) * 60.0;
					weight = w->weight_link(link, time, N, &(w->BP[link][r]));
				sum_weight += weight;

				curr = next;
				time += weight;
				t2 = (int)((time - start_time) / period_length + 0.5);
				if (t2 == t1)
					t2++;
				if (t2 > TIME_PERIOD - 1)
					t2 = TIME_PERIOD - 1;
				next = node[curr].next_node_noi[t2];
				t1 = t2;
			}
			if (curr != destIdx)
			{
				printf("No path between origin %d and destination %d for departure time %d on day %d.\n", N->NodeId(orig), destID, dep, r);
				exit(1);
			}
			exp_weight += sum_weight * w->prob[r];
			fprintf(fp, "%12.2f\n", sum_weight);
		}
		//sg		fprintf(fp, "%4s%20d%20d%20d%20d%20.6f\n", " ", N->NodeId(orig), destID, dep, policy, exp_weight);
	}
	if (op_flag)
	{
		opsize = op_links.size();
		//links = new float[opsize];
		links = new double[opsize];
		sorted = new int[opsize];
		for (i = 0; i < opsize; i++)
		{
			links[i] = op_links[i];
			sorted[i] = i;
		}
		heapSort(links, sorted, opsize);

		op_links.clear();
		op_links.push_back(links[sorted[0]]);
		for (i = 1; i < opsize; i++)
			if (fabs(links[sorted[i]] - links[sorted[i - 1]]) > EPSILON)
				op_links.push_back(links[sorted[i]]);

		delete links;
		delete sorted;
	}

	fclose(fp);
	return(op_links);
}

//x = {x1, x2, ...}, combination returns a vector, where each element y = (y1, y2, ...) 
//is a vector with the same length as x, and y_j = 0, 1, ..., x_j-1
vector<vector<int>> combination(vector<int> x)
{
	vector<vector<int>> result, result2;
	vector<int>			tmp, y;
	vector<int>::iterator	it, it3;
	vector<vector<int>>::iterator it2;
	unsigned int i, j;

	if (x.size() == 1)
	{
		for (i = 0; i < x.at(0); i++)
		{
			tmp.push_back(i);
			result.push_back(tmp);
			tmp.clear();
		}
	}
	else
	{ //use recursive call
		y = x;
		it = y.end();
		it--;
		y.erase(it);

		result2 = combination(y);

		for (it2 = result2.begin(); it2 != result2.end(); it2++)
		{
			it3 = x.end();
			it3--;
			for (j = 0; j < (*it3); j++)
			{
				tmp = (*it2);
				tmp.push_back(j);
				result.push_back(tmp);
			}
		}
	}
	return result;
}

/*
vector<long> Policy2Path_NOI(Network *N, EdgeWeights *w, long destIdx, vector<pair<long, int>> orig_dep, node_info *node, int start_time, int policy, char *filename, int op_flag)
{
//This function convert a routing policy to a set of paths with their associated travel times and probabilities in a network
//	represented by marginal distributions of time-dependent travel time random variables
//The current implementation might not be optimal, as it simply enumerates all possible combinations of marginal realizations along the policy

//key data structures
vector<pair<long, long>>	reachables_queue;
int	**reachable_link;  //flag to denote whether a dynamic link is reachable by a policy AND the current support point index
vector<pair<long, long>>		reachables;
vector<int>					reachables_numsupports;
vector<double>				path_probs;
vector<vector<int>>			result;
pair<long, long>			node_time;
pair<long, long>			link_time;

vector<pair<long,long>>     reachable_nodes; //reachable node-time pairs

long	i, j, orig, nodes, destID, currnode, nextnode, currtime, nexttime, timeslice, curr, next, link, numlinks,
t1, t2, TIME_PERIOD = w->levels();
int		dep, period_length = w->agglength(), *sorted, opsize, flag;
float   exp_weight, sum_weight, weight, *links;
FILE	*fp;
vector<long>	op_links;

if (policy == 0)
fp = fopen(filename, "w");
else
fp = fopen(filename, "a");

if (!fp)
{
printf("Error: cannot open file %s.\n", filename);
exit(1);
}

nodes = N->NumNodes();
numlinks = N->NumLinks();

destID = N->NodeId(destIdx);

reachable_link = new int*[numlinks];
for (i = 0; i < numlinks; i++)
{
reachable_link[i] = new int[w->numtimeslices()];
for (j = 0; j < w->numtimeslices(); j++)
reachable_link[i][j] = 0;
}


for (i = 0; i < orig_dep.size(); i++) //for each origin-departure pair
{

//breadth search to find all reacheable (node, time) pairs  (here time is the index into time slices of travel time step functions)
t1 = (long)((orig_dep[i].second - start_time) / period_length + 0.5);
node_time.first = orig_dep[i].first;
node_time.second = t1;
reachables_queue.clear();
//reachable_nodes.clear();
reachables_queue.push_back(node_time);  //put the origin in the queue
//reachable_nodes.push_back(node_time);
while (!reachables_queue.empty())
{
//pop out a node-time pair
node_time = reachables_queue.front();
reachables_queue.erase(reachables_queue.begin());
currnode = node_time.first;
if (currnode == destIdx)
continue;
t1 = node_time.second;
timeslice = w->find_interval(0, start_time + t1 * period_length);

//add outgoing (node,time) pairs from the current node-time pair
nextnode = node[currnode].next_node_noi[t1];
link = N->LinkIndex(currnode, nextnode);
reachable_link[link][timeslice] = 1;

node_time.first = nextnode;
for (j = 0; j < w->numsupportpoints[link][timeslice]; j++)
{
if (!(w->probability[link][timeslice][j] > 0))
continue;
weight = w->supportpoints[link][timeslice][j];
t2 = (long)((t1 * period_length + weight)/ period_length + 0.5);
if (t2 == t1)
t2++;
if (t2 > TIME_PERIOD - 1)
t2 = TIME_PERIOD - 1;
node_time.second = t2;
//check whether the element already exists in the queue
flag = 0;
for (vector<pair<long,long>>::iterator x = reachables_queue.begin(); x != reachables_queue.end(); x++)
if (*x == node_time)
{
flag = 1;
break;
}
if (!flag)
reachables_queue.push_back(node_time);
flag = 0;
//for (vector<pair<long,long>>::iterator x = reachable_nodes.begin(); x != reachable_nodes.end(); x++)
//	if (*x == node_time)
//	{
//		flag = 1;
//		break;
//	}
//if (!flag)
//	reachable_nodes.push_back(node_time);
}
}


//vector<pair<long,long>> next_nodes;
//for (vector<pair<long,long>>::iterator x = reachable_nodes.begin(); x != reachable_nodes.end(); x++)
//{
//	node_time = *x;
//	currnode = node_time.first;
//	t1 = node_time.second;
//	node_time.second = node[currnode].next_node_noi[t1];
//	next_nodes.push_back(node_time);
//}

//////////////////////////
//build a vector to store the reachable dynamic links
reachables.clear();
reachables_numsupports.clear();
for (link = 0; link < numlinks; link++)
for (timeslice = 0; timeslice < w->numtimeslices(); timeslice++)
{
if (reachable_link[link][timeslice] && !w->det_flag[link][timeslice])
{
link_time.first = link;
link_time.second = timeslice;
reachables.push_back(link_time);
reachables_numsupports.push_back(w->numsupportpoints[link][timeslice]);
}
}

vector<int>::iterator it2;
vector<pair<long, long>>::iterator it3;
long	max_numsupports = 0;
for (it2 = reachables_numsupports.begin(); it2 != reachables_numsupports.end(); it2++)
max_numsupports = (long) max(max_numsupports, *it2);

//get the combination of support point indices
//result = combination(reachables_numsupports);

//random sample realizations + selective samples (to include extreme usually low prob realizations)

//find the corresponding path for each sample
orig = orig_dep[i].first;
dep = orig_dep[i].second;

exp_weight = 0;
//for (vector<vector<int>>::iterator it = result.begin(); it != result.end(); it++)
srand((unsigned)time( NULL ));

vector<long> path;
int samplesize = 100;  //number of joint realizations sampled
double p;
for (i = 0; i < numlinks; i++)
for (j = 0; j < w->numtimeslices(); j++)
reachable_link[i][j] = 0;

for (i = 0; i < samplesize; i++)
{
//set the support point index for the current sample
for (it2 = reachables_numsupports.begin(), it3 = reachables.begin(); it2 != reachables_numsupports.end(); it2++, it3++)
{
if (i >= max_numsupports) //end of selective samples
p = rand()/(double)(RAND_MAX + 1);
else
p = (double) i / (double) max_numsupports;
long	prob_idx = floor(p * (*it2));
if (prob_idx >= *it2)
prob_idx = *it2 - 1;
reachable_link[(*it3).first][(*it3).second] = prob_idx;
}

fprintf(fp, "%10d%10d%10d", N->NodeId(orig), destID, dep);

curr = orig;
currtime = dep;
t2 = (long)((currtime - start_time) / period_length + 0.5);
if (t2 == t1)
t2++;
if (t2 > TIME_PERIOD - 1)
t2 = TIME_PERIOD - 1;
next = node[curr].next_node_noi[t2];

sum_weight = 0;
while (next != -1)
{
link = N->LinkIndex(curr, next);
if (link == -1)
{
printf("Wrong routing policy (node %d, time period %d : next node %d) to destination %d.\n",
N->NodeId(curr), t2, N->NodeId(next), destID);
exit(1);
}
fprintf(fp, "%10d", N->LinkId(link));
if (op_flag)
op_links.push_back(link);

timeslice = w->find_interval(0, currtime);
weight = w->supportpoints[link][timeslice][reachable_link[link][timeslice]];
sum_weight += weight;

curr = next;
currtime += weight;
t2 = (long)((currtime - start_time) / period_length + 0.5);
if (t2 == t1)
t2++;
if (t2 > TIME_PERIOD - 1)
t2 = TIME_PERIOD - 1;
next = node[curr].next_node_noi[t2];
t1 = t2;
}
if (curr != destIdx)
{
printf("No path between origin %d and destination %d for departure time %d.\n", N->NodeId(orig), destID, dep);
exit(1);
}
//exp_weight += sum_weight * prob;
//fprintf(fp, "%12.2f\n", sum_weight);
fprintf(fp, "\n");
//path_probs.push_back(prob);
}
//sg		fprintf(fp, "%4s%20d%20d%20d%20d%20.6f\n", " ", N->NodeId(orig), destID, dep, policy, exp_weight);
}
if (op_flag)
{
opsize = op_links.size();
links = new float[opsize];
sorted = new int[opsize];
for (i = 0; i < opsize; i++)
{
links[i] = op_links[i];
sorted[i] = i;
}
heapSort(links, sorted, opsize);

op_links.clear();
op_links.push_back(links[sorted[0]]);
for (i = 1; i < opsize; i++)
if (fabs(links[sorted[i]] - links[sorted[i - 1]]) > EPSILON)
op_links.push_back(links[sorted[i]]);

delete links;
delete sorted;
}

fclose(fp);
return(op_links);
}
*/


void write_result_simple(Network *N, EdgeWeights *w, int destIndex, node_info *node, char *filename)
{
	//need to calculate adaptive_cost_label_ce and true_cost_label_noi

	FILE *fp;
	register int TIME_PERIOD, NUM_NODES, NUM_ARCS, NUM_REAL, event, n, i;
	int nextnode;

	fp = fopen(filename, "a");
	if (!fp)
	{
		printf("Error: can't open %s\n", filename);
		exit(1);
	}

	NUM_NODES = N->NumNodes();
	NUM_ARCS = N->NumLinks();
	NUM_REAL = w->numreal();
	TIME_PERIOD = w->levels();

	assert(destIndex >= 0 && destIndex < NUM_NODES);
	//assert( w->minweight() > 0 );

	//write the destination node
	fprintf(fp, "%d\n", N->NodeId(destIndex));

	//write next node for each node at each time period   
	for (n = 0; n < NUM_NODES; n++)
	{
		for (i = 0; i < TIME_PERIOD; i++)
		{
			nextnode = (node[n].next_node_noi)[i];
			if (nextnode >= 0)
				fprintf(fp, "%d\n", N->NodeId(nextnode));
			else
				fprintf(fp, "%d\n", nextnode);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}

/*
vector<long> Policy2Path_CDPI(Network *N, EdgeWeights *w, long destIdx, vector<pair<long, int>> orig_dep, node_info *node, int start_time, int policy, char *filename, int op_flag)
{
//This function convert a routing policy to a set of paths with their associated travel times and probabilities in a network
//	represented by marginal distributions of time-dependent travel time random variables
//The current implementation might not be optimal, as it simply enumerates all possible combinations of marginal realizations along the policy

//key data structures
vector<pair<long, long>>	reachables_queue;
int	**reachable_link;  //flag to denote whether a dynamic link is reachable by a policy AND the current support point index
vector<pair<long, long>>		reachables;
vector<int>					reachables_numsupports;
vector<double>				path_probs;
vector<vector<int>>			result;
pair<long, long>			node_time;
pair<long, long>			link_time;

vector<pair<long,long>>     reachable_nodes; //reachable node-time pairs

long	i, j, orig, nodes, destID, currnode, nextnode, currtime, nexttime, timeslice, curr, next, link, numlinks,
t1, t2, TIME_PERIOD = w->levels();
int		dep, period_length = w->agglength(), *sorted, opsize, flag;
float   exp_weight, sum_weight, weight, *links;
FILE	*fp;
vector<long>	op_links;

if (policy == 0)
fp = fopen(filename, "w");
else
fp = fopen(filename, "a");

if (!fp)
{
printf("Error: cannot open file %s.\n", filename);
exit(1);
}

nodes = N->NumNodes();
numlinks = N->NumLinks();

destID = N->NodeId(destIdx);

reachable_link = new int*[numlinks];
for (i = 0; i < numlinks; i++)
{
reachable_link[i] = new int[w->numtimeslices()];
for (j = 0; j < w->numtimeslices(); j++)
reachable_link[i][j] = 0;
}


for (i = 0; i < orig_dep.size(); i++) //for each origin-departure pair
{

//breadth search to find all reacheable (node, time) pairs  (here time is the index into time slices of travel time step functions)
t1 = (long)((orig_dep[i].second - start_time) / period_length + 0.5);
node_time.first = orig_dep[i].first;
node_time.second = t1;
reachables_queue.clear();
//reachable_nodes.clear();
reachables_queue.push_back(node_time);  //put the origin in the queue
//reachable_nodes.push_back(node_time);
while (!reachables_queue.empty())
{
//pop out a node-time pair
node_time = reachables_queue.front();
reachables_queue.erase(reachables_queue.begin());
currnode = node_time.first;
if (currnode == destIdx)
continue;
t1 = node_time.second;
timeslice = w->find_interval(0, start_time + t1 * period_length);

//add outgoing (node,time) pairs from the current node-time pair
nextnode = node[currnode].next_node_noi[t1];

link = N->LinkIndex(currnode, nextnode);
reachable_link[link][timeslice] = 1;

node_time.first = nextnode;
for (j = 0; j < w->numsupportpoints[link][timeslice]; j++)
{
if (!(w->probability[link][timeslice][j] > 0))
continue;
weight = w->supportpoints[link][timeslice][j];
t2 = (long)((t1 * period_length + weight)/ period_length + 0.5);
if (t2 == t1)
t2++;
if (t2 > TIME_PERIOD - 1)
t2 = TIME_PERIOD - 1;
node_time.second = t2;
//check whether the element already exists in the queue
flag = 0;
for (vector<pair<long,long>>::iterator x = reachables_queue.begin(); x != reachables_queue.end(); x++)
if (*x == node_time)
{
flag = 1;
break;
}
if (!flag)
reachables_queue.push_back(node_time);
flag = 0;
//for (vector<pair<long,long>>::iterator x = reachable_nodes.begin(); x != reachable_nodes.end(); x++)
//	if (*x == node_time)
//	{
//		flag = 1;
//		break;
//	}
//if (!flag)
//	reachable_nodes.push_back(node_time);
}
}


//vector<pair<long,long>> next_nodes;
//for (vector<pair<long,long>>::iterator x = reachable_nodes.begin(); x != reachable_nodes.end(); x++)
//{
//	node_time = *x;
//	currnode = node_time.first;
//	t1 = node_time.second;
//	node_time.second = node[currnode].next_node_noi[t1];
//	next_nodes.push_back(node_time);
//}

//////////////////////////
//build a vector to store the reachable dynamic links
reachables.clear();
reachables_numsupports.clear();
for (link = 0; link < numlinks; link++)
for (timeslice = 0; timeslice < w->numtimeslices(); timeslice++)
{
if (reachable_link[link][timeslice] && !w->det_flag[link][timeslice])
{
link_time.first = link;
link_time.second = timeslice;
reachables.push_back(link_time);
reachables_numsupports.push_back(w->numsupportpoints[link][timeslice]);
}
}

vector<int>::iterator it2;
vector<pair<long, long>>::iterator it3;
long	max_numsupports = 0;
for (it2 = reachables_numsupports.begin(); it2 != reachables_numsupports.end(); it2++)
max_numsupports = (long) max(max_numsupports, *it2);

//get the combination of support point indices
//result = combination(reachables_numsupports);

//random sample realizations + selective samples (to include extreme usually low prob realizations)

//find the corresponding path for each sample
orig = orig_dep[i].first;
dep = orig_dep[i].second;

exp_weight = 0;
//for (vector<vector<int>>::iterator it = result.begin(); it != result.end(); it++)
srand((unsigned)time( NULL ));

vector<long> path;
int samplesize = 100;  //number of joint realizations sampled
double p;
for (i = 0; i < numlinks; i++)
for (j = 0; j < w->numtimeslices(); j++)
reachable_link[i][j] = 0;

for (i = 0; i < samplesize; i++)
{
//set the support point index for the current sample
for (it2 = reachables_numsupports.begin(), it3 = reachables.begin(); it2 != reachables_numsupports.end(); it2++, it3++)
{
if (i >= max_numsupports) //end of selective samples
p = rand()/(double)(RAND_MAX + 1);
else
p = (double) i / (double) max_numsupports;
long	prob_idx = floor(p * (*it2));
if (prob_idx >= *it2)
prob_idx = *it2 - 1;
reachable_link[(*it3).first][(*it3).second] = prob_idx;
}

fprintf(fp, "%10d%10d%10d", N->NodeId(orig), destID, dep);

curr = orig;
currtime = dep;
t2 = (long)((currtime - start_time) / period_length + 0.5);
if (t2 == t1)
t2++;
if (t2 > TIME_PERIOD - 1)
t2 = TIME_PERIOD - 1;
next = node[curr].next_node_noi[t2];

sum_weight = 0;
while (next != -1)
{
link = N->LinkIndex(curr, next);
if (link == -1)
{
printf("Wrong routing policy (node %d, time period %d : next node %d) to destination %d.\n",
N->NodeId(curr), t2, N->NodeId(next), destID);
exit(1);
}
fprintf(fp, "%10d", N->LinkId(link));
if (op_flag)
op_links.push_back(link);

timeslice = w->find_interval(0, currtime);
weight = w->supportpoints[link][timeslice][reachable_link[link][timeslice]];
sum_weight += weight;

curr = next;
currtime += weight;
t2 = (long)((currtime - start_time) / period_length + 0.5);
if (t2 == t1)
t2++;
if (t2 > TIME_PERIOD - 1)
t2 = TIME_PERIOD - 1;
next = node[curr].next_node_noi[t2];
t1 = t2;
}
if (curr != destIdx)
{
printf("No path between origin %d and destination %d for departure time %d.\n", N->NodeId(orig), destID, dep);
exit(1);
}
//exp_weight += sum_weight * prob;
//fprintf(fp, "%12.2f\n", sum_weight);
fprintf(fp, "\n");
//path_probs.push_back(prob);
}
//sg		fprintf(fp, "%4s%20d%20d%20d%20d%20.6f\n", " ", N->NodeId(orig), destID, dep, policy, exp_weight);
}
if (op_flag)
{
opsize = op_links.size();
links = new float[opsize];
sorted = new int[opsize];
for (i = 0; i < opsize; i++)
{
links[i] = op_links[i];
sorted[i] = i;
}
heapSort(links, sorted, opsize);

op_links.clear();
op_links.push_back(links[sorted[0]]);
for (i = 1; i < opsize; i++)
if (fabs(links[sorted[i]] - links[sorted[i - 1]]) > EPSILON)
op_links.push_back(links[sorted[i]]);

delete links;
delete sorted;
}

fclose(fp);
return(op_links);
}
*/
//vector<long> write_route_choice_result_new(Network *N, EdgeWeights *w, long destIdx, vector<pair<long, int>> orig_dep, node_info *node, int start_time, int policy, char *filename, int op_flag, int **realizationToEvent)

/*
* searches for a value in sorted array
*   arr is an array to search in
*   value is searched value
*   left is an index of left boundary
*   right is an index of right boundary
* returns position of searched value, if it presents in the array
* or -1, if it is absent
*/
//int binarySearch(float arr[], int value, int left, int right) {
int binarySearch(double arr[], int value, int left, int right) {
	while (left <= right) {
		int middle = (left + right) / 2;
		if (arr[middle] == value)
			return middle;
		else if (arr[middle] > value)
			right = middle - 1;
		else
			left = middle + 1;
	}
	return -1;
}


vector<long> write_route_choice_result_new(Network *N, EdgeWeights *w, long destIdx, vector<pair<long, int>> orig_dep, link_info *link, int start_time, int policy, char *filename, char *filename1, int op_flag, int **realizationToEvent, double runtime, double runtime1, double runtime3)
{
	//calculate the true expected travel time and obtain corresponding paths (for now, only set of links)
	long	orig, nodes, destID, curr, next, TIME_PERIOD = w->levels();
	//long linkindex;
	int		dep, period_length = w->agglength(), i, r, t1, t2, *sorted, opsize;
	//float   exp_weight, sum_weight, weight, time, *links;
	float   exp_weight, sum_weight, weight, time;
	//float *links2, *links;
	double *links2, *links;
	int temp = 0, exist = 0;
	FILE	*fp;
	FILE	*fp1;
	vector<long>	op_links2;
	long orig_next, curr_link, outdegree, orig_link;
	long optimal_next, optimal_link;
	int m, j, k, n, l, q;
	//float **sum_weight_all; 
	//int *optimal_m;
	int blocked_index = 0;
	//int blocked_day=0;

	//op_links2.push_back(100);
	//return op_links2;

	nodes = N->NumNodes();
	destID = N->NodeId(destIdx);

	if (policy == 0)
	{
		fp = fopen(filename, "w");
		fp1 = fopen(filename1, "w");
	}
	else
	{
		fp = fopen(filename, "a");
		fp1 = fopen(filename1, "a");
	}

	if (!fp)
	{
		printf("Error: cannot open file %s.\n", filename);
		exit(1);
	}


	if (!fp1)
	{
		printf("Error: cannot open file %s.\n", filename1);
		exit(1);
	}


	//sg	if (policy == 0)
	//sg		fprintf(fp, "%4s%20s%20s%20s%20s%20s\n", "//", "Origin", "Destination", "Departure", "Policy#", "Mean");

	/*
	if (op_flag==1)
	{
	fprintf(fp,"ORP Total Runtime:\n%f\n",runtime);
	fprintf(fp,"ORP Generate Event Collection Runtime:\n%f\n",runtime1);
	//fprintf(fp,"ORP Static Shortest Running Time:\n%f\n",runtime2);
	fprintf(fp,"ORP Label Correcting Runtime:\n%f\n",runtime3);
	}
	*/

	for (i = 0; i < orig_dep.size(); i++)
	{
		orig = orig_dep[i].first;
		dep = orig_dep[i].second;
		//orig_next = N->LinkDest(orig, optimal_m);
		//orig_link = N->LinkIndex(orig, orig_next);					
		exp_weight = 0;
		for (r = 0; r < w->numreal(); r++)
		{
			blocked_index = 0;
			sum_weight = 0;
			fprintf(fp, "%16d%16d%16d", N->NodeId(orig), destID, dep);
			fprintf(fp1, "%16d%16d%16d", N->NodeId(orig), destID, dep);
			outdegree = N->OutDegree(orig);
			for (m = 0; m < outdegree; m++)
			{
				orig_next = N->LinkDest(orig, m);
				orig_link = N->LinkIndex(orig, orig_next);
				if (N->Disabled(orig_link))
					continue;
				int i1, start1, end1;
				bp start_bp1;
				int dep_time = (int)((dep - start_time) / period_length + 0.5);
				if (dep_time > TIME_PERIOD - 1)
					dep_time = TIME_PERIOD - 1;
				if (w->det_flag[orig_link] == 1)
				{
					if (w->blocked[orig_link][0][0] == 1)
					{
						continue;
					}
				}
				else
				{
					//start_bp1 = w->find_interval(&(w->BP[orig_link][realizationToEvent[dep_time][r]]), dep, &i1);
					start_bp1 = w->find_interval(&(w->BP[orig_link][r]), dep, &i1);
					start1 = i1 - 1;
					end1 = i1;
					//if (w->blocked[orig_link][start1][realizationToEvent[dep_time][r]]==1||w->blocked[orig_link][end1][realizationToEvent[dep_time][r]]==1)
					if (w->blocked[orig_link][start1][r] == 1 || w->blocked[orig_link][end1][r] == 1)
					{
						continue;
					}
				}
				optimal_next = N->LinkDest(orig, m);
				optimal_link = N->LinkIndex(orig, optimal_next);
				break;
			}
			for (m = 0; m < outdegree; m++)
			{
				//blocked_index=0;
				orig_next = N->LinkDest(orig, m);
				orig_link = N->LinkIndex(orig, orig_next);
				if (N->Disabled(orig_link))
					continue;
				int i1, start1, end1;
				bp start_bp1;
				int dep_time = (int)((dep - start_time) / period_length + 0.5);
				if (dep_time > TIME_PERIOD - 1)
					dep_time = TIME_PERIOD - 1;
				if (w->det_flag[orig_link] == 1)
				{
					if (w->blocked[orig_link][0][0] == 1)
					{
						continue;
					}
				}
				else
				{
					//start_bp1 = w->find_interval(&(w->BP[orig_link][realizationToEvent[dep_time][r]]), dep, &i1);
					start_bp1 = w->find_interval(&(w->BP[orig_link][r]), dep, &i1);
					start1 = i1 - 1;
					end1 = i1;
					//if (w->blocked[orig_link][start1][realizationToEvent[dep_time][r]]==1||w->blocked[orig_link][end1][realizationToEvent[dep_time][r]]==1)
					if (w->blocked[orig_link][start1][r] == 1 || w->blocked[orig_link][end1][r] == 1)
					{
						continue;
					}
				}
				if (m > 0)
				{
					if (*(*(link[orig_link].cost_label + dep_time) + realizationToEvent[dep_time][r]) < *(*(link[optimal_link].cost_label + dep_time) + realizationToEvent[dep_time][r]))
					{
						optimal_link = orig_link;
					}
				}


			}
			orig_link = optimal_link;
			if (N->Disabled(orig_link))
				continue;
			int i1, start1, end1;
			bp start_bp1;
			int dep_time = (int)((dep - start_time) / period_length + 0.5);

			if (dep_time > TIME_PERIOD - 1)
				dep_time = TIME_PERIOD - 1;

			if (w->det_flag[orig_link] == 1)
			{
				if (w->blocked[orig_link][0][0] == 1)
				{
					continue;
				}
			}
			else
			{
				//start_bp1 = w->find_interval(&(w->BP[orig_link][realizationToEvent[dep_time][r]]), dep, &i1);
				start_bp1 = w->find_interval(&(w->BP[orig_link][r]), dep, &i1);
				start1 = i1 - 1;
				end1 = i1;
				//if (w->blocked[orig_link][start1][realizationToEvent[dep_time][r]]==1||w->blocked[orig_link][end1][realizationToEvent[dep_time][r]]==1)
				if (w->blocked[orig_link][start1][r] == 1 || w->blocked[orig_link][end1][r] == 1)
				{
					continue;
				}
			}
			curr = orig;
			curr_link = orig_link;
			time = dep;
			t1 = (int)((time - start_time) / period_length + 0.5);
			//if (t2 == t1)
			//t2++;
			if (t1 > TIME_PERIOD - 1)
				t1 = TIME_PERIOD - 1;
			//next = node[curr].next_node_noi[t2];
			//next = node[curr].next_node[t1][realizationToEvent[t1][r]];			
			//orig_next = N->LinkDest(orig, 0);
			//curr_link = N->LinkInIndex(orig, orig_next);
			next = link[curr_link].next_node[t1][realizationToEvent[t1][r]];


			if (next == -1)
			{
				fprintf(fp, " No path between origin %d and destination %d for departure time %d on day %d.\n", N->NodeId(orig), destID, dep, r);
				fprintf(fp1, " No path between origin %d and destination %d for departure time %d on day %d.\n", N->NodeId(orig), destID, dep, r);
				continue;
			}


			// t1 = t2;
			while (curr != N->nodeIDtoIndex(destID))
				//while (next != -1)
			{
				/*
				if (next == -1)
				{
				fprintf(fp, " No path between origin %d and destination %d for departure time %d on day %d.\n", N->NodeId(orig), destID, dep, r);
				blocked_day=1;
				break;
				}
				*/


				//linkindex = N->LinkIndex(curr, next);
				//if (N->Disabled(linkindex))
				if (curr_link == -1)
				{
					printf("Wrong routing policy (node %d, time period %d : next node %d) to destination %d.\n",
						N->NodeId(curr), t2, N->NodeId(next), destID);
					exit(1);
				}
				if (N->Disabled(curr_link))
					continue;


				int i2, start2, end2;
				bp start_bp2;
				if (w->det_flag[curr_link] == 1)
				{
					if (w->blocked[curr_link][0][0] == 1)
					{
						blocked_index = 1;
						break;
					}
				}
				else
				{
					//start_bp2 = w->find_interval(&(w->BP[curr_link][realizationToEvent[t1][r]]), time, &i2);
					start_bp2 = w->find_interval(&(w->BP[curr_link][r]), start_time + t1 * period_length, &i2);
					start2 = i2 - 1;
					end2 = i2;
					//if (w->blocked[curr_link][start2][realizationToEvent[t1][r]]==1||w->blocked[curr_link][end2][realizationToEvent[t1][r]]==1)
					if (w->blocked[curr_link][start2][r] == 1 || w->blocked[curr_link][end2][r] == 1)
					{
						blocked_index = 1;
						break;
					}
				}
				if (op_flag == 1)
				{
					if (r == temp)
					{
						op_links2.push_back(curr_link);
					}
					else
					{
						opsize = op_links2.size();
						//links = new float[opsize];
						links = new double[opsize];
						sorted = new int[opsize];
						//links2 = new float[opsize];
						links2 = new double[opsize];
						for (q = 0; q < opsize; q++)
						{
							links[q] = op_links2[q];
							sorted[q] = q;
						}
						heapSort(links, sorted, opsize);
						for (q = 0; q < opsize; q++)
						{
							links2[q] = links[sorted[q]];
						}
						int left = 0;
						int right = opsize - 1;
						exist = binarySearch(links2, curr_link, left, right);
						if (exist == -1)
						{
							op_links2.push_back(curr_link);
						}
						delete links;
						delete sorted;
						delete links2;
					}
					if (N->LinkDest(curr_link) == destIdx)
					{
						temp = r;
					}
				}

				if (w->det_flag[curr_link] == 1)
					//weight = w->weight_link(curr_link, time, N, &(w->BP[curr_link][0])) * 60.0;	  //travel time is in minutes
					weight = w->weight_link(curr_link, start_time + t1 * period_length, N, &(w->BP[curr_link][0]));
				else
					//weight = w->weight_link(curr_link, time, N, &(w->BP[curr_link][realizationToEvent[t1][r]])) * 60.0;
					//weight = w->weight_link(curr_link, time, N, &(w->BP[curr_link][realizationToEvent[t1][r]]));
					weight = w->weight_link(curr_link, start_time + t1 * period_length, N, &(w->BP[curr_link][r]));
				sum_weight += weight;

				fprintf(fp, "%16d", N->LinkId(curr_link));
				fprintf(fp1, "%16.2f", weight);

				//curr = N->LinkDest(curr_link);			
				//time += weight;  
				//t2 = (int)((time - start_time) / period_length + 0.5); //si she wu ru
				t2 = (int)(t1 + weight / period_length + 0.5);
				//if (t2 == t1)
				//t2++;
				if (t2 > TIME_PERIOD - 1)
					t2 = TIME_PERIOD - 1;
				t1 = t2;
				//next = node[curr].next_node_noi[t2];
				curr = N->LinkDest(curr_link);
				if (curr == next)
					continue;
				curr_link = N->LinkIndex(curr, next);
				next = link[curr_link].next_node[t1][realizationToEvent[t1][r]];
			}
			if (curr != destIdx)
			{
				//if(blocked_day==1)
				//continue;
				if (blocked_index == 0)
				{
					printf("No path between origin %d and destination %d for departure time %d on day %d.\n", N->NodeId(orig), destID, dep, r);
					exit(1);
				}
				else
					continue;
			}
			exp_weight += sum_weight * w->prob[r];
			fprintf(fp, "%24.2f\n", sum_weight);
			fprintf(fp1, "%24.2f\n", sum_weight);
		}
		//sg		fprintf(fp, "%4s%20d%20d%20d%20d%20.6f\n", " ", N->NodeId(orig), destID, dep, policy, exp_weight);
		//fprintf(fp, "Expected %4s%20d%20d%20d%20d%20.6f\n", " ", N->NodeId(orig), destID, dep, policy, exp_weight);
	}

	fclose(fp);
	fclose(fp1);
	return op_links2;

}
