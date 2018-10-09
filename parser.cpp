/*
  Network IO Routines
  Brian Dean
  1997
  Revised by Song Gao May 2004
  Added a function to read OD demand
  Revised the function to read stochastic time-dependent weights
*/


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include "parser.h"
#include "network.h"
#include "od.h"
#include "HeapSort.h"

Network *readnetwork_simple(char *filename)
{
  	char s[256], s2[256];
  	FILE *fp;
  	long Nnodes=0, Nlinks=0, x, y, z, temp, ID, index1, index2, linkIndex;
  	float x1, y1, x2, y2;
  	float t;

  	fp = fopen( filename, "rt" );
  	if( !fp ) 
  		{
    	printf( "Error: can't open %s\n", filename );
    	exit(1);
  		}

  	while( fgets( s, 255, fp ) ) 
    	if( sscanf( s, " [Nodes] : %d", &Nnodes ) == 1 )   //   ???  Nnodes
      		break;
  	while( fgets( s, 255, fp ) ) 
    	if( sscanf( s, " [Links] : %d", &Nlinks ) == 1 ) 
      		break;
  	fclose( fp );
  	if( !Nnodes || !Nlinks ) 
  		{
    	printf( "Invalid network file: %s\n", filename );
    	exit(1);
  		}

  	Network *N = new Network(0, Nnodes, Nlinks);
  	fp = fopen( filename, "rt" );

  	while( fgets( s, 255, fp ) ) 
    	if( sscanf( s, " [Nodes] : %d", &Nnodes ) == 1 ) //???
      		break;
  	for( x=0; x<Nnodes; x++ )
    	while( fgets( s, 255, fp ) ) 
      		if( sscanf( s, " { %d }", &ID) == 1 ) 
      			{
				N->AddNode( ID );// N 是结构体变量  指向 ？
				//printf("node ID: %d\n", ID);
				break;	
      			}

  	while( fgets( s, 255, fp ) ) 
    	if( sscanf( s, " [Links] : %d", &Nlinks ) == 1 ) 
      		break;
  	for( x = 0; x < Nlinks; x++ )
    	while( fgets( s, 255, fp ) ) 
			if( sscanf( s, " { %d %d %d }", &ID, &y, &z) == 3 )
				{
				index1 = N->nodeIDtoIndex(y);//  what‘s N??
				index2 = N->nodeIDtoIndex(z);

				//printf("ID: %d\t UpperNode: %d(%d)\t DownNode: %d(%d)\n", ID, y, index1, z, index2);
				//for( index1=0; N->NodeId(index1) != y && index1 < Nnodes; index1++ );
				//for( index2=0; N->NodeId(index2) != z && index2 < Nnodes; index2++ );
				assert( index1 < Nnodes && index2 < Nnodes);//？？？
				//if( !N->LinkExists(index1, index2) )
				linkIndex = N->AddLink(index1, index2, ID);//？？
				//else
				//printf("Link between nodes %d and %d already exists in the network.\n", y, z);
				break;	
				}
	  
	N->IndexNetwork();  
	fclose( fp );
	return N;

}

Network *readnetwork_mitsim(char *filename) 
{
  	char s[256], s2[256];
  	FILE *fp;
  	long Nnodes=0, Nlinks=0, x, y, z, temp, ID, index1, index2, linkIndex;
  	float x1, y1, x2, y2;
  	float t;

  	fp = fopen( filename, "rt" );
  	if( !fp ) 
  		{
    	printf( "Error: can't open %s\n", filename );
    	exit(1);
  		}

  	while( fgets( s, 255, fp ) ) 
    	if( sscanf( s, " [Nodes] : %d", &Nnodes ) == 1 ) 
      		break;
  	while( fgets( s, 255, fp ) ) 
    	if( sscanf( s, " [Links] : %d", &Nlinks ) == 1 ) 
      		break;
  	fclose( fp );
  	if( !Nnodes || !Nlinks ) 
  		{
    	printf( "Invalid network file: %s\n", filename );
    	exit(1);
  		}

  	Network *N = new Network(0, Nnodes, Nlinks);
  	fp = fopen( filename, "rt" );

  	while( fgets( s, 255, fp ) ) 
    	if( sscanf( s, " [Nodes] : %d", &Nnodes ) == 1 ) 
      		break;
  	for( x=0; x<Nnodes; x++ )
    	while( fgets( s, 255, fp ) ) 
      		if( sscanf( s, " { %d %d \"%s\"}", &ID, &temp, s2 ) == 3 ) 
      			{
				N->AddNode( ID );
				//printf("node ID: %d\n", ID);
				break;	
      			}

  	while( fgets( s, 255, fp ) ) 
    	if( sscanf( s, " [Links] : %d", &Nlinks ) == 1 ) 
      		break;
  	for( x = 0; x < Nlinks; x++ )
    	while( fgets( s, 255, fp ) ) 
      		if( !(sscanf( s, " { %f %f %f %f %f %s", &t, &t, &t, &t, &t, s2 ) == 6 && s2[0]=='}' ) )
				if( sscanf( s, " { %d %d %d %d %d", &ID, &temp, &y, &z, &temp ) == 5 ) 
					{
					index1 = N->nodeIDtoIndex(y);
					index2 = N->nodeIDtoIndex(z);

					//printf("ID: %d\t UpperNode: %d(%d)\t DownNode: %d(%d)\n", ID, y, index1, z, index2);
					//for( index1=0; N->NodeId(index1) != y && index1 < Nnodes; index1++ );
					//for( index2=0; N->NodeId(index2) != z && index2 < Nnodes; index2++ );
					assert( index1 < Nnodes && index2 < Nnodes);
					//if( !N->LinkExists(index1, index2) )
					linkIndex = N->AddLink(index1, index2, ID);
					//else
					//printf("Link between nodes %d and %d already exists in the network.\n", y, z);
					while( fgets( s, 255, fp ) )
						if( sscanf( s, " { %f %f %f %f %f %s", &x1, &y1, &t, &x2, &y2, s2 ) == 6 && s2[0]=='}' ) 
							{
					  		N->SetLinkLength(linkIndex, sqrt(pow(x1-x2,2)+pow(y1-y2,2)));
					  		break;
							}
					break;	
					}
	  
	N->IndexNetwork();  
	fclose( fp );
	return N;

}

EdgeWeights *readlinktimes(Network *N, char *filename, long levels, long agglength, long numreal) 
{
  
  FILE *fp;
  long r, t, m;
  float tt;
  long temp; //to store unused data
  long y;
  EdgeWeights *w;

  fp = fopen( filename, "rt" );
  if( !fp ) {
    printf( "Error: can't open %s\n", filename );
    exit(1);
  }

  //fscanf( fp, "%d %d %d", &levels, &agglength, &numreal );
  w = new EdgeWeights( N, levels, agglength, numreal);

	w->realization = new float*[levels * N->NumLinks()];
	w->prob = new float[numreal];
	for(long x = 0; x < levels * N->NumLinks(); x++)
		{
		w->realization[x] = new float[numreal];
		for (long y = 0; y < numreal; y++)
			w->realization[x][y] = 0.0;
		}
		
	for( r=0; r<numreal; r++) {
    fscanf( fp, "%d %f", &temp, w->prob+r);
    //printf("prob[%d]: %f\n", r, w->prob[r]);
    for( t=0; t<levels; t++ )
      	{
		for( m=0; m<N->NumLinks(); m++ ) 
			{
	  		y = N->linkIDtoIndex(N->linkIdbyUser(m));
            fscanf( fp, " %f", &tt );
	  
	  		w->setweight_link( y, t, r, tt/agglength, N);

			}
     	}
  }
  fclose( fp );
  return w;
}

/*
EdgeWeights *ReadStepLinkTimes(Network *N, char *filename, long levels, long agglength) //agglength is in number of seconds 
{
	char s[256], s2[256];
  	FILE *fp;
  	long i, x, t, UserID, index, numpoints, j, k, numtimeslices = 0, numlinks = 0, numsto = 0, starttime, tmp;
	float sum, tt;
  	EdgeWeights *w;

  	fp = fopen(filename, "rt" );
  	if(!fp) 
  		{
    	printf( "Error: can't open %s\n", filename );
    	exit(1);
  		}

  	w = new EdgeWeights(N, levels, agglength);

	//read start times of time slices
  	while (fgets(s, 255, fp)) 
    	if (sscanf(s, " [Dynamic Dividing Points] : %d", &tmp) == 1 ) 
      		break;
	numtimeslices = tmp - 1;
	if (numtimeslices >= 1)
		{
		w->slicestart = new long[numtimeslices+1];
		w->setnumslices(numtimeslices);

		for (i = 0; i < numtimeslices; i++)
			while (fgets(s, 255, fp))
				if (sscanf(s, " { %d }", &w->slicestart[i]) == 1)
					break;
		}

	//read free flow travel times
	while (fgets(s, 255, fp))
		if (sscanf(s, " [Static Travel Times] : %d", &numlinks) == 1)
			break;

	if (numlinks != N->NumLinks())
	{
		printf("Wrong number of links in section [Free flow Travel Time]!");
		exit(1);
	}

	w->supportpoints = new long**[numlinks];
	w->probability = new float**[numlinks];
	w->avgtime = new float[numlinks];
	w->numsupportpoints = new long*[numlinks];
	w->det_flag = new short*[numlinks];
	for (i = 0; i < numlinks; i++)
	{
		w->det_flag[i] = new short[numtimeslices];
		for (j = 0; j < numtimeslices; j++)
			w->det_flag[i][j] = 1;
	}

	for (i = 0; i < numlinks; i++)
	{
		while (fgets(s, 255, fp))
			if (sscanf(s, " {%d %f}", &UserID, &tt) == 2) //travel time in seconds
			{
				index = N->linkIDtoIndex(UserID);
				if (index < 0)
					{
					printf("Link ID %d does not exist in the network!\n", N->linkIdbyUser(UserID));
					exit(1);
					}
				w->avgtime[index] = tt > 1 ? tt : 1;
				w->numsupportpoints[index] = new long[numtimeslices];
				w->supportpoints[index] = new long*[numtimeslices];
				w->probability[index] = new float*[numtimeslices];
				break;
			}
	}	

	//read stochastic dynamic links
	while (fgets(s, 255, fp))
		if (sscanf(s, " [Stochastic Dynamic Links] : %d", &numsto) == 1)
			break;

	for (i = 0; i < numsto; i++)
		while (fgets(s, 255, fp))
			if (sscanf(s, " {%d", &UserID) == 1)
				{
				index = N->linkIDtoIndex(UserID);

				if (index < 0)
					{
					printf("Link ID %d does not exist in the network!\n", N->linkIdbyUser(UserID));
					exit(1);
					}
				
				for (x = 0; x < numtimeslices; x++)
				{
					while (fgets(s, 255, fp))
						if (s[0] == '}')
							goto out_of_loop;
						else if (sscanf(s, " [%d,%d) : %d ", &starttime, &tmp, &numpoints) == 3)
							break;
					t = w->find_interval(x, starttime);

					w->probability[index][t] = new float[numpoints];
					w->numsupportpoints[index][t] = numpoints;
					w->supportpoints[index][t] = new long[numpoints];

					for (j = 0; j < numpoints; j++)
					{
						while (fgets(s, 255, fp))
							if (sscanf(s, " {%d %f}", &tmp, &w->probability[index][t][j]) == 2)
							{
								w->supportpoints[index][t][j] = tmp > 1 ? tmp : 1;
								break;
							}
					}
					w->det_flag[index][t] = 0;
				}
out_of_loop:
				break;
			   	}

	//set numsupportpoints and supportpoints for deterministic link-time
	for (i = 0; i < numlinks; i++)
		for (j = 0; j < numtimeslices; j++)
		{
			if (w->det_flag[i][j])
			{
				w->numsupportpoints[i][j] = 1;
				w->supportpoints[i][j] = new long[1];
				w->supportpoints[i][j][0] = w->avgtime[i];
				w->probability[i][j] = new float[1];
				w->probability[i][j][0] = 1.;
			}
		}

  	fclose(fp);
  	return w;
}
*/

/*
int ***blocked;

int Blocked(long linkIndex, int t, int r)
{
	return blocked[linkIndex][t][r];
}
void SetBlocked(long linkIndex, int t, int r, int state)
{
	blocked[linkIndex][t][r] = state;
}
*/

EdgeWeights *ReadPiecewiseLinkTimes(Network *N, char *filename, long levels, long agglength) //agglength is in number of seconds 
{
	char s[256], s2[256];
  	FILE *fp;
    int default_flag;
  	long i, UserID, index, numpoints, j, k, numdefaults = 0, numdet = 0, numsto = 0, numreal = 0;
	long m,n,p;
	float sum;
  	EdgeWeights *w;
	bp	breakingpoint;
	
  	fp = fopen(filename, "rt" );
  	if(!fp) 
  		{
    	printf( "Error: can't open %s\n", filename );
    	exit(1);
  		}

  	w = new EdgeWeights(N, levels, agglength);

	//read default breaking points
  	while (fgets(s, 255, fp)) 
    	if (sscanf(s, " [Default Breaking Points] : %d", &numdefaults) == 1 ) 
      		break;
	if (numdefaults >= 1)
		{
		w->defaults = new long[numdefaults];
		w->setnumdefaults(numdefaults);

		for (i = 0; i < numdefaults; i++)
			while (fgets(s, 255, fp))
				if (sscanf(s, " { %d }", &w->defaults[i]) == 1)
					break;
		}

	//read deterministic links
	while (fgets(s, 255, fp))
		if (sscanf(s, " [Deterministic Links] : %d", &numdet) == 1)
			break;
	while (fgets(s, 255, fp))
        if (sscanf(s, " [Stochastic Links] : %d", &numsto) == 1)
			break;
  	fclose(fp);
	if (numdet + numsto != N->NumLinks())
		{
		printf("The sum of numbers of deterministic and stochastic links is not equal to the number of links.\n");
		exit(1);
		}

    w->det_flag = new short[N->NumLinks()];
	
	for (i = 0; i < N->NumLinks(); i++)
	{
			w->det_flag[i] = -1;
	}

  	fp = fopen(filename, "rt");
	while (fgets(s, 255, fp))
		if (sscanf(s, " [Deterministic Links] : %d", &numdet) == 1)
			break;
	for (i = 0; i < numdet; i++)
		while (fgets(s, 255, fp))
			if (sscanf(s, " {%d", &UserID) == 1)
				{
				index = N->linkIDtoIndex(UserID);
				if (index < 0)
					{
					printf("Link ID %d does not exist in the network!\n", N->linkIdbyUser(UserID));
					exit(1);
					}
				else if (w->det_flag[index] != -1)
					{
					printf("Duplicate data for link ID %d.\n", UserID);
					exit(1);
					}
				w->det_flag[index] = 1;

				while (fgets(s, 255, fp))
					if (sscanf(s, "     {%d", &numpoints) == 1)
						break;

				default_flag = 0;

				if (numpoints <= 0 && numdefaults >= 1)
					{
					numpoints = numdefaults;  //use the default breaking points
					default_flag = 1;
					}
				if (numpoints <= 0)
					{
					printf("Number of breaking points must be at least 1.\n");
					exit(1);
					}
				w->BP[index] = new bps[1];
				w->blocked[index]=new int*[1];
				w->blocked[index][0]=new int[1];
				w->blocked[index][0][0]=0;
				
				for (j = 0; j < numpoints; j++)
					while (fgets(s, 255, fp))
						if (!default_flag && sscanf(s, "         { %d %f }", &breakingpoint.first, &breakingpoint.second) == 2)
						{
							if(breakingpoint.second==-1)
							{
								w->blocked[index][0][0]=1;
							}
							w->BP[index][0].push_back(breakingpoint);
							break;
			   			}
						else if (default_flag && sscanf(s, "         {  %f }", &breakingpoint.second) == 1)
						{
							if(breakingpoint.second==-1)
							{
								w->blocked[index][0][0]=1;
							}
							breakingpoint.first = w->defaults[j];
							w->BP[index][0].push_back(breakingpoint);
							break;
			   			}					
				break;
			   	}	

	//read support point probabilities
	while (fgets(s, 255, fp))
		if (sscanf(s, " [Support Points] : %d", &numreal) == 1)
			break;

	if (numreal >= 1)
		{
		w->setnumreal(numreal);
		w->prob = new float[numreal];
		sum = 0.;
		for (i = 0; i < numreal; i++)
			while (fgets(s, 255, fp))
				if (sscanf(s, " %f ", w->prob + i) == 1)
					{
					sum += w->prob[i];
					break; 
					}
		for (i = 0; i < numreal; i++)
			w->prob[i] /= sum;
		}

	//read stochastic links
	while (fgets(s, 255, fp))
		if (sscanf(s, " [Stochastic Links] : %d", &numsto) == 1)
			break;
	for (i = 0; i < numsto; i++)
		while (fgets(s, 255, fp))
			if (sscanf(s, " {%d", &UserID) == 1)
				{
				index = N->linkIDtoIndex(UserID);

				if (index < 0)
					{
					printf("Link ID %d does not exist in the network!\n", N->linkIdbyUser(UserID));
					exit(1);
					}
				else if (w->det_flag[index] != -1)
					{
					printf("Duplicate data for link ID %d.\n", UserID);
					exit(1);
					}
				w->det_flag[index] = 0;
				w->BP[index] = new bps[numreal];

				int count=0;
				for (j = 0; j < numreal; j++)
					{
					while (fgets(s, 255, fp))
						if (sscanf(s, "     {%d", &numpoints) == 1)
							break;

					default_flag = 0;
					if (numpoints <= 0 && numdefaults >= 1)
						{
						numpoints = numdefaults;  //use the default breaking points
						default_flag = 1;
						}
					if (numpoints <= 0)
						{
						printf("Number of breaking points must be at least 1.\n");
						exit(1);
						}
					
					/*
					int ***blocked;
					blocked=new int**[N->NumLinks()];
					for (m=0;m<N->NumLinks();m++)
					{
						blocked[m]=new int*[numpoints];
					}

					for (m=0;m<N->NumLinks();m++)
					{
						for(n=0;n<numpoints;n++)
						{
							blocked[m][n]=new int[numreal];
						}
					}
					for (m=0;m<N->NumLinks();m++)
					{
						for(n=0;n<numpoints;n++)
						{
							for(p=0;n<numreal;n++)
							{
								blocked[m][n][p]=0;
							}
						}
					}
					*/

					if(count==0)
					{
						w->blocked[index]=new int*[numpoints];
						for(n=0;n<numpoints;n++)
						{
							w->blocked[index][n]=new int[numreal];
						}
						for(n=0;n<numpoints;n++)
						{
							for(p=0;p<numreal;p++)
							{
								w->blocked[index][n][p]=0;
							}
						}
					}
					count=1;

					for (k = 0; k < numpoints; k++)
						while (fgets(s, 255, fp))
							if (!default_flag && sscanf(s, "         { %d %f }", &breakingpoint.first, &breakingpoint.second) == 2)
							{
								if(breakingpoint.second==-1)
								{
									w->blocked[index][k][j]=1;
								}
								w->BP[index][j].push_back(breakingpoint);
								break;
				   			}
							else if (default_flag && sscanf(s, "         {  %f }", &breakingpoint.second) == 1)
							{
								if(breakingpoint.second==-1)
								{	
									w->blocked[index][k][j]=1;
								}
								breakingpoint.first = w->defaults[k];
								w->BP[index][j].push_back(breakingpoint);
								break;
				   			}					
					}
				break;
			   	}

  	fclose(fp);
  	return w;
}


odpairs ReadDynamicODPairs(Network *N, char *filename)  //for Borlange data
{
	FILE 					*fp;
	odpairs					ods;
	pair<long, int> 		orig_dep;
	vector<pair<long, int>>	orig_dep_bydest;	
	long					destID, destIdx, *orig;
	int						tmp, count, *dep, *sorted, i, curr;
    //float                   *dest;
	double                   *dest;
    char                    s[256];

	fp = fopen(filename, "rt");
	if (!fp)
		{
		printf("Error: cannot open %s.\n", filename);
		exit(1);
		}

	count = 0;
	while (fgets(s, 255, fp))
		count++;
	fclose(fp);

	if (count == 0)
		{
		printf("Error: the od file %s is empty.\n", filename);
		exit(1);
		}

	orig = new long[count];
	//dest = new float[count];
	dest = new double[count];
	dep = new int[count];
	sorted = new int[count];

	fp = fopen(filename, "rt");
	count = 0;
	while (fgets(s, 255, fp))
		{
		if (sscanf(s, "%d %lf %d %d\n", &orig[count], &dest[count], &tmp, &dep[count]) != 4)
			{
			printf("Error reading %d.\n", filename);
			exit(1);
			}
		count++;
		}
	for (i = 0; i < count; i++)
		sorted[i] = i;

	heapSort(dest, sorted, count);

	for (i = 0; i < count; i++)
		{
		curr = sorted[i];
		orig_dep.first = N->nodeIDtoIndex(orig[curr]);
		orig_dep.second = dep[curr];
		orig_dep_bydest.push_back(orig_dep);
		if (i == count - 1 || dest[curr] != dest[sorted[i + 1]])
			{
			destIdx = N->nodeIDtoIndex(dest[curr]);
			ods[destIdx] = orig_dep_bydest;
			orig_dep_bydest.clear();
			}
		}
    return ods;
}


ODDemand *readdemand(char *filename, long levels, long stepsize, Network *N) {

  FILE *fp;

  ODDemand *ods = new ODDemand(levels);
  
  char s[256], s2[256];
  long starttime, time, orig, dest, origIndex, destIndex, demand=0, x, thisdemand;
  float scale;
  vector<long> times;

  fp = fopen( filename, "rt" );
  if( !fp ) {
    printf( "Error: can't open %s\n", filename );
    exit(1);
  }
  
  while( fgets( s, 255, fp ) )
    if( sscanf( s, "%d %d %f", &time, &x, &scale ) == 3 )
      times.push_back(time);
 
  fclose(fp);

  assert( times.back() - times.front() <= levels*stepsize);
  starttime = times[0];
  times.push_back(levels*stepsize+starttime); //one more time point to mark the end
  
  ods->setStart(starttime);

  fopen( filename, "rt" );
  
  long counter = 0;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, "%d %d %f", &x, &x, &scale ) == 3 ) {
      while ( fgets( s, 255, fp) ) {
	if ( s[0] == '}') break;
	if ( sscanf( s, " { %d %d %d %s", &orig, &dest, &thisdemand, s2 ) == 4 && s2[0] == '}' ) {
	  origIndex = N->nodeIDtoIndex(orig);
      if (origIndex == -1)
      {
          printf("Origin node ID %d in demand file %s not found in the network.", orig, filename);
          exit(1);
      }
      
	  destIndex = N->nodeIDtoIndex(dest);
      if (destIndex == -1)
      {
          printf("Destination node ID %d in demand file %s not found in the network.", dest, filename);
          exit(1);
      }

	  for (long i=(times[counter]-starttime)/stepsize;i<(times[counter+1]-starttime)/stepsize;i++) {

	    ods->addOD(origIndex, destIndex, i);
	    ods->addDemand(thisdemand, i);
	  }
	  demand += (long)(ceil(scale*thisdemand*(times[counter+1]-times[counter])/3600));
	  
	}

      } 
      counter++;
    }

  fclose(fp);
  


  ods->find_list_of_dest();
  ods->setNumVeh(demand);

  //ods->printOD();

  return ods;
}
  

