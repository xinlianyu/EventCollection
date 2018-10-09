/*
  Network Class
  Brian Dean
  6/27/97
  Revised by Song Gao May 2004
  Add functions related to path calculation
*/

#include <vector>
#include <utility>
#include <assert.h>
#include <stdio.h>
#include "network.h"


#define BADID "4294967295"

void Network::DebugOut( void ) {
  long x;
  printf( "(V %d) (E %d) (Indexed %s)\n", Nnodes, Nlinks, indexed ? "Yes" : "No" );
  printf( "\nNode# Indeg OutDeg   ID firstA firstB\n" );
  for( x=0; x<Nnodes; x++ )
    printf( "%5d%6d%7d%5d%7d%7d\n", x, indegree[x], outdegree[x], Nids[x], firstA[x], firstB[x] );
  printf( "\nLink# Source Dest   ID destsorted\n" );
  for( x=0; x<Nlinks; x++ )
    printf( "%5d%7d%5d%5d%11d\n", x, source[x], dest[x], Lids[x], indexed ? destsorted[x] : 0 );
}

#define iswap( x, y ) { t=x; x=y; y=t; }
void Network::RowSwap( long i1, long i2 )
{
  long t;
  iswap( source[i1], source[i2] );
  iswap( dest[i1], dest[i2] );
  iswap( Lids[i1], Lids[i2] );
}

void Network::IndexRows( void ) {
  long *temp = new long[Nnodes], x, y;
  firstB[0] = temp[0] = 0;
  for( x=1; x<Nnodes; x++ )
    firstB[x] = temp[x] = indegree[x-1] + firstB[x-1];
  for( x=0; x<Nlinks; x++ ) {
    y = dest[x];
    destsorted[temp[y]++] = x;
  }  
  indexed = 1;

  //deal with links with no outgoing links and therefore firstA is not set
  firstA[Nnodes] = Nlinks;
  for (x = Nnodes - 1; x >= 0; x--)
	  if (firstA[x] == -1)
		  firstA[x] = firstA[x + 1];

}

#define pcomp( s1,d1,s2,d2 ) (s1>s2 || (s1==s2 && d1>d2))

long Network::LinkIndex(long n1, long n2) 
{
    long i;

//4_10  	//printf("entering LinkIndex(long, long)...\n");
//4_10  	long first = 0;
//4_10  	//printf("after first line of LinkIndex(long, long)...\n");
//4_10  	//printf("Nlinks: %d\n", Nlinks);
//4_10  
//4_10  	long last = Nlinks-1;
//4_10  	// printf("after second line of LinkIndex(long, long)...\n");
//4_10  	long mid;
//4_10
//4_10  	while( first <= last ) 
//4_10  		{
//4_10    	mid = (first + last)/2;
//4_10    	if( n1==source[mid] && n2==dest[mid] )
//4_10      		return mid;
//4_10    	if( pcomp( n1,n2,source[mid],dest[mid] ) )
//4_10      		first = mid+1;
//4_10    	else
//4_10      		last = mid-1;
//4_10  		}
//4_10  	return -1;
	firstA[Nnodes] = Nlinks;
	for (i = firstA[n1]; i < firstA[n1 + 1]; i++)
		{
		if (n2 == dest[i])
			return i;
		}
	return -1;		
}

Network::Network( long _Nnodes, long _maxnodes=-1, long _maxlinks=-1 ) {
  	long x;
  	if( _maxnodes == -1 )
    	_maxnodes = _Nnodes;
  	if( _maxlinks == -1 )
    	_maxlinks = _Nnodes * _Nnodes;
  	maxnodes = _maxnodes;
  	maxlinks = _maxlinks;
  	Nlinks = Nnodes = indexed = 0;
  	source =     new long[maxlinks];
 	dest =       new long[maxlinks];
  	destsorted = new long[maxlinks];
  	Lids =       new long[maxlinks];
  	Lids_byuser = new long[maxlinks];
  	linklength = new float[maxlinks];
  	disabled	= new int[maxlinks];
  
  	for (x = 0; x < maxlinks; x++)
		disabled[x] = 0;

  	firstA     = new long[maxnodes + 1];
  	firstB     = new long[maxnodes];
  	indegree   = new long[maxnodes];
  	outdegree  = new long[maxnodes];
  	Nids       = new long[maxnodes];
  	for( x=0; x<_Nnodes; x++ )
    	AddNode(-1);
}

Network::~Network() {
  delete source;
  delete dest;
  delete destsorted;
  delete Lids;
  delete Lids_byuser;
  delete firstA;
  delete firstB;
  delete indegree;
  delete outdegree;
  delete Nids;
}

void Network::AddNode( long Nid=-1 ) {
  assert( Nnodes < maxnodes );
  indegree[Nnodes] = outdegree[Nnodes] = 0;
  firstA[Nnodes] = firstB[Nnodes] = -1;
  Nids[Nnodes] = Nid;
  Nnodes++;
  indexed = 0;
}

void Network::RemoveNode( long index ) {
  long x;
  assert( index>=0 && index<Nnodes );
  for( x=0; x<Nlinks; x++ ) 
    if( source[x]==index || dest[x]==index ) {
      RemoveLink( source[x], dest[x] );
      x--;
    }
  for( x=0; x<Nlinks; x++ ) {
    if( source[x] > index ) 
      source[x]--;
    if( dest[x] > index )
      dest[x]--;
  }
  for( x=index; x<Nnodes-1; x++ ) {
    indegree[x] = indegree[x+1];
    outdegree[x] = outdegree[x+1];
    Nids[x] = Nids[x+1];
    firstA[x] = firstA[x+1];
  }
  Nnodes--;
  indexed = 0;
}

long Network::AddLink(long indexSource, long indexDest, long Lid = -1) 
{
  	long i, x;
  	assert( Nlinks < maxlinks );
  	//assert( !LinkExists( indexSource, indexDest ) );
  	outdegree[indexSource]++;
  	indegree[indexDest]++;
  	source[Nlinks] 	= indexSource;
  	dest[Nlinks]   	= indexDest;
  	Lids[Nlinks]   	= Lid;
	Lids_byuser[Nlinks] = Lid;
  	i = Nlinks;
  	while( i>=1 && pcomp( source[i-1], dest[i-1], source[i], dest[i] ) ) 
  		{
    	i--;
    	RowSwap( i, i+1 );
  		}
  	Nlinks++;
  	for( x=0; x<Nnodes; x++ )
    	if( firstA[x] >= i && x != indexSource )
      		firstA[x]++;
  	if( firstA[indexSource] == -1 )
    	firstA[indexSource] = i;
  	indexed = 0;
  
  	return i;
}

void Network::RemoveLink( long indexSource, long indexDest ) {
  long i, x;
  assert( LinkExists( indexSource, indexDest ) );
  i = LinkIndex( indexSource, indexDest );
  outdegree[indexSource]--;
  indegree[indexDest]--;
  if( outdegree[indexSource] == 0 )
    firstA[indexSource] = -1;
  for( x=0; x<Nnodes; x++ )
    if( firstA[x] > i )
      firstA[x]--;
  while( i < Nlinks-1 ) {
    RowSwap( i, i+1 );
    i++;
  }
  Nlinks--;
  indexed = 0;
}

Network *Network::Copy( void ) {
  long x;

  Network *N = new Network(0, maxnodes, maxlinks );
  for( x=0; x<Nnodes; x++ )
    N->AddNode( Nids[x] );
  for( x=0; x<Nlinks; x++ )
    N->AddLink( LinkSource(x), LinkDest(x), Lids[x] );
  return N;
}

Network *Network::Inverse( void ) {
  long x, y;

  Network *N = new Network(0, Nnodes, Nnodes*Nnodes-Nlinks );
  for( x=0; x<Nnodes; x++ )
    N->AddNode( Nids[x] );
  for( x=0; x<Nnodes; x++ )
    for( y=0; y<Nnodes; y++ )
      if( !LinkExists( x, y ) )
	N->AddLink( x, y );
  return N;
}

Network *Network::Dual( void ) {
  long x, y, z, Tlinks=0;
  
  for( x=0; x<Nnodes; x++ )
    Tlinks += indegree[x] * outdegree[x];
  Network *N = new Network(0, Nlinks, Tlinks );
  for( x=0; x<Nlinks; x++ )
    N->AddNode( Lids[x] );
  for( x=0; x<Nnodes; x++ )
    for( y=0; y<indegree[x]; y++ )
      for( z=0; z<outdegree[x]; z++ ) {
	N->AddLink( LinkIndex(LinkSource(y,x),x), LinkIndex(x,LinkDest(z,x)) );
      }
  return N;  
}

//the following three functions: generatecompletepaths only apply in acyclic networks

void Network::generatecompletepaths(char *pathfile, vector<long> destIndex) 
{
  
	long linkIndex, jnode, nextlink, i, k, l;

	FILE *fp;
	fp = fopen(pathfile, "wt");
	if (!fp) 
		{
		printf("Error: Cannot open file %s!\n", pathfile);
		exit(1);
		}

	fprintf(fp, "%d\n", destIndex.size());
	for (vector<long>::iterator it = destIndex.begin(); it != destIndex.end(); it++) 
		{
		fprintf(fp, "%d\n", *it);
		fprintf(fp, "%d\n", NumLinks());
		for (i = 0;i < NumLinks(); i++) 
			{
			linkIndex = linkIDtoIndex(linkIdbyUser(i));
			fprintf(fp, "%d\n", numpaths(linkIndex, *it));
			jnode = LinkDest(linkIndex);
			if (jnode == *it) 
				fprintf(fp, "%s %s\n", BADID, BADID);
			else 
				for (k = 0; k < OutDegree(jnode); k++) 
					{
					nextlink = firstA[jnode] + k;
					for (l = 0; l < numpaths(nextlink, *it); l++)
						fprintf(fp, "%d %d\n", LinkId(nextlink), l);
					}
			}
		}

	fclose(fp);
}
 
/*
long Network::numpaths(long linkIndex, long destIndex) {
  //we use a recursive call to compute the number of paths out of any given link

  long num = 0;
  long jnode = LinkDest(linkIndex);

  if (jnode == destIndex)
    return 1;
  else {
    for (long i=0;i<OutDegree(jnode);i++) 
      num += numpaths(firstA[jnode]+i, destIndex);
    return num;
  }

}
*/

long Network::numpaths(long origIndex, long destIndex) {
 
	long currlink, counter = 0, i;

	//printf("Origin node: %d\t Origin link: %d\n", origin, LinkId(links[0]));

	PT::iterator itt;

	for (itt = PathTables.begin(); itt != PathTables.end(); itt++) 
		{
		//printf("destination: %d\n", (*itt).first);
		if ((*itt).first == destIndex) {
	  	//printf("dest: %d\n", (*itt).first);
		  	for (i = 0; i < OutDegree(origIndex); i++) 
		  		{
				currlink = firstA[origIndex] + i;
				//printf("current link ID: %d\n", LinkId(currlink));
				counter += (*itt).second[LinkId(currlink)].size();
		  		}
		  	break;
			}
		}
	return counter;
}


long Network::numpaths_byID(long origID, long destID) {
  return numpaths(nodeIDtoIndex(origID), nodeIDtoIndex(destID));
}

/*
long Network::links2path(vector<long> links, long destIndex) {
  long counter=0;
  long currnode, currlink;

  for (vector<long>::iterator it=links.begin();it!=links.end();it++) {
    currnode = LinkSource(*it);

    for (long i=0;i<OutDegree(currnode);i++) {
      currlink = firstA[currnode]+i;
      if (currlink==*it) break;
      counter += numpaths(currlink, destIndex);
    }
  }

  assert(counter < numpaths_node(LinkSource(links[0]), destIndex));
  return counter;
}
*/

void Network::buildpathtable(char *pathfile) 
{
 
	vector< vector< pair< long, long > > > PathTable;
	pair<long, long> SubPath;
	vector< pair< long, long > > LinkTable;

	long numdest, currdest, numlinks, subpathsize, currlink, nextlink;
	long counter = 0, maxnumpath = 0, tmp, i, j, k;


	FILE *fp;
	fp = fopen(pathfile, "rt");
	if(!fp) 
		{
		printf("Error: Cannot open file %s!\n", pathfile);
		exit(1);
		}

	fscanf(fp, "%d", &numdest); 
  	for (k = 0; k < numdest; k++) 
  		{
    	fscanf(fp, "%d", &currdest);
    	fscanf(fp, "%d", &numlinks);
    	assert(numlinks == Nlinks);
    	for (i = 0; i < numlinks; i++) 
    		{
      		fscanf(fp, "%d", &subpathsize);
			//printf("subpathsize: %d\n", subpathsize);
      		for (j = 0; j < subpathsize; j++) 
      			{
				fscanf(fp, "%d %d", &SubPath.first, &SubPath.second);
				//printf("[nextlink nextpath]: %d %d\n", SubPath.first, SubPath.second);
				LinkTable.push_back(SubPath);
      			}
      		//printf("LinkTable size: %d\n", LinkTable.size());
      		PathTable.push_back(LinkTable);
      		LinkTable.erase(LinkTable.begin(), LinkTable.end());
    		}
    	PathTables[currdest] = PathTable;
    	PathTable.erase(PathTable.begin(), PathTable.end());
  		}
  	fclose(fp);

}

//get the path index given a set of link indices
long Network::links2path(vector<long> links, long destIndex) 
{
	//links are in User ID
	//my program does not have system ID
	//DynaMIT order it by system ID
	//So make sure that User ID and system IDs are the same in DynaMIT network file
	long subpathsize, currlink, nextlink;
	long counter = 0, i, origin = LinkSource(links[0]), appnum;
  	//printf("Origin node: %d\t Origin link: %d\n", origin, LinkId(links[0]));
  	PT::iterator itt;
	vector<long>::iterator it;

	for (itt = PathTables.begin(); itt != PathTables.end(); itt++) 
		{
		//printf("destination: %d\n", (*itt).first);
		if ((*itt).first == destIndex) 
			{
	  		//printf("dest: %d\n", (*itt).first);
	  		for (i = 0; i < OutDegree(origin); i++) 
	  			{
				currlink = firstA[origin] + i;
				//printf("current link ID: %d\n", LinkId(currlink));
				if (LinkId(currlink) < LinkId(links[0]))
	  				counter += (*itt).second[LinkId(currlink)].size();
	  			}
	  		break;
			}
		}

	for (it = links.begin(); it != links.end(); it++) 
		{
		currlink = *it;
		if (LinkDest(currlink) == destIndex) 
			break;

		subpathsize = (*itt).second[LinkId(currlink)].size();
		for (i = 0; i < subpathsize; i++) 
			{
	  		nextlink = (*itt).second[LinkId(currlink)][i].first;
	  		if (nextlink == LinkId(*(it+1)))
				break;
	  		counter++;
			}
		if (subpathsize == 0 || i >= subpathsize) {//what if the path is not found???
		  	//find an approximate path??
		  	printf("Error: path not found!\n");
		  	printf("Return an approximate path instead.\n");
		  	appnum = counter - (*itt).second[LinkId(currlink)].size();
		  	//appnum = appnum<(numlinks-1)?appnum:numlinks-1;
		  	return appnum;
			}
		}
	return counter;
}

vector< long >  Network::path2links(long origIndex, long destIndex, long pathIndex) {

  if (numpaths(origIndex, destIndex) < 1) {
    printf("no paths between %d and %d!\n", origIndex, destIndex);
    exit(1);
  }

  assert(origIndex!=destIndex);
  assert(pathIndex>=0 && pathIndex < numpaths(origIndex, destIndex));

 

  vector< long > links;
  long currlink, nextlink, currsize, localindex;
  long linkcounter = 0;

  PT::iterator itt;

  //printf("origin: %d destIndex: %d pathindex: %d\n", origIndex, destIndex, pathIndex);

  for (itt=PathTables.begin();itt!=PathTables.end();itt++) {
    //printf("destination: %d\n", (*itt).first);
    if ((*itt).first==destIndex) {
      currsize=0;
      //printf("dest: %d\n", (*itt).first);
      for (linkcounter=0;linkcounter<(*itt).second.size();linkcounter++) {
	if (LinkSource(linkIDtoIndex(linkcounter))==origIndex) {
	  currsize += (((*itt).second)[linkcounter]).size();
	  //printf("linkcounter: %d\n", linkcounter);
	}
	if (currsize > pathIndex) 
	  break;
      }
      
      //printf("just outside the loop...\n");

      if (linkcounter == (*itt).second.size())
	linkcounter--;

      currsize -= (*itt).second[linkcounter].size();
      localindex = pathIndex - currsize;
      
      currlink = linkIDtoIndex(linkcounter);
      links.push_back(currlink);
      
      //printf("currsize: %d\t localindex: %d\t currlink: %d currlinkID: %d\n", currsize, localindex, currlink, linkcounter);

      while (LinkDest(currlink)!=destIndex) {
	nextlink = linkIDtoIndex((*itt).second[linkcounter][localindex].first);
	localindex = (*itt).second[linkcounter][localindex].second;

	currlink = nextlink;
	links.push_back(currlink);

	linkcounter = LinkId(nextlink);
      }
      break;
    }
  }
  
  return links;
}

void Network::PrintPath(long origIndex, long destIndex, long pathIndex) {
  
  vector< long > links = path2links(origIndex, destIndex, pathIndex);
  for (vector< long >::iterator it=links.begin();it!=links.end();it++)
    printf("%d ", LinkId(*it));
  printf("\n");

}

void Network::PrintPath() {
  PT::iterator it;
  vector< vector< pair< long, long > > >::iterator it2;
  vector< pair< long, long > >::iterator it3;
  long counter=0;
  
  for (it=PathTables.begin();it!=PathTables.end();it++) {
    printf("Destination: %d\n", (*it).first);
    counter=0;
    for (it2=(*it).second.begin();it2!=(*it).second.end();it2++) {
      printf("***From Link %d***\n", counter++);
      for (it3=(*it2).begin();it3!=(*it2).end();it3++)
	printf("Next Link: %d\t Next Path: %d\n", (*it3).first, (*it3).second);	
      printf("\n");
    }
    printf("\n");
  }
}
