/*
  Network Class
  Brian Dean
  6/27/97
  Revised by Song Gao May 2004
  Add functions related to path calculation
*/
#pragma comment(lib,"ws2_32.lib")
#ifndef NETWORK
#define NETWORK
#include "assert.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <map>

using namespace std;

typedef map< long, vector< vector< pair< long, long > > > >  PT;

class Network {

//I need to rewrite some of the network structure to distinguish between User ID, System ID and system index.
  
private:

  long 	maxnodes, maxlinks, Nnodes, Nlinks, *destsorted, *Nids, *Lids, *Lids_byuser, *source, *dest, *firstA, *firstB, *indegree, *outdegree;
  int	indexed, *disabled;
  float *linklength;
  void 	RowSwap( long i1, long i2 );
  void 	IndexRows( void );
  PT 	PathTables; //a very important feature of PathTables: the links are in IDs, not in indices, the nodes are in indicies
public:
  
  Network( long _Nnodes, long _maxnodes, long _maxlinks);
  ~Network();

//long  NumNodes   ( void );
  void AddNode    ( long );
  void RemoveNode ( long index );
//long  NodeId     ( long index );
//long  InDegree   ( long index );
//long  OutDegree  ( long index );

//long  NumLinks   ( void );
  long AddLink    ( long indexSource, long indexDest, long LinkId );
  void RemoveLink ( long indexSource, long indexDest );
//long  LinkId     ( long indexSource, long indexDest );
//long  LinkId     ( long index );
//long  LinkExists ( long indexSource, long indexDest );
//long  LinkSource ( long i, long nodeIndex=-1 );
//long  LinkDest   ( long i, long nodeIndex=-1 );
  long  LinkIndex   ( long indexSource, long indexDest );
//long  LinkOutIndex( long nodeIndex, long i );
//long  LinkInIndex ( long nodeIndex, long i );

  
  Network *Copy( void );
  Network *Inverse( void );
  Network *Dual( void );

  void DebugOut(void);  
//void IndexNetwork( void );

  void buildpathtable(char *pathfile);
  void generatecompletepaths(char *pathfile, vector<long> destIndex);
  long numpaths(long origIndex, long destIndex);
  long numpaths_byID(long origId, long destID);
  long links2path(vector<long> links, long destIndex);
  vector< long > path2links(long origIndex, long destIndex, long pathIndex);
  void PrintPath( void );
  void PrintPath(long origIndex, long destIndex, long pathIndex);

  // Inline functions:

  inline void IndexNetwork( void ) {
    IndexRows();
  }

  inline long LinkOutIndex( long nodeIndex, long i ) {
    assert( indexed );
    return firstA[nodeIndex]+i;
  }
  
  inline long LinkInIndex( long nodeIndex, long i ) {
    assert( indexed );
    return destsorted[firstB[nodeIndex]+i];
  }

  inline long Network::NumNodes( void ) {
    return Nnodes;
  }

  inline long Network::InDegree( long index ) {
    assert( index>=0 && index<Nnodes );
    return indegree[index];  
  }

  inline long Network::OutDegree( long index ) {
    assert( index>=0 && index<Nnodes );
    return outdegree[index];
  }

  inline long Network::NumLinks( void ) {
    return Nlinks;
  }

  inline long Network::LinkExists( long indexSource, long indexDest ) {
    assert( indexSource>=0 && indexSource<Nnodes );
    assert( indexDest>=0 && indexDest<Nnodes );
    return LinkIndex( indexSource, indexDest ) != -1;
  }

  inline long Network::NodeId( long index ) {
    assert( index>=0 && index<Nnodes );
    return Nids[index];
  }

  inline long Network::LinkId( long indexSource, long indexDest ) {
    assert( LinkExists( indexSource, indexDest ) );
    return Lids[LinkIndex(indexSource,indexDest)];
  }

  inline long Network::LinkId( long index ) {
    assert( index>=0 && index<Nlinks );
    return Lids[index];
  }

  	inline long Network::nodeIDtoIndex( long ID) 
  		{
    	long index;
    	for (index = 0; index < Nnodes && NodeId(index) != ID; index++);
    	if (index < Nnodes)
        	return index;
    	else
        	return -1;
  		}

  	inline long Network::linkIDtoIndex( long ID) 
  		{
    	long index;
    	for (index = 0; index < Nlinks && LinkId(index) != ID; index++);
    	if (index < Nlinks) 
        	return index;
    	else
        	return -1;
  		}

	inline long Network::linkIdbyUser(long user_idx)
		{
		assert(user_idx >= 0 && user_idx < Nlinks);
		return Lids_byuser[user_idx];
		}

  inline long Network::LinkSource( long nodeIndex, long i=-1 ) {
    if( i==-1 ) {
      assert( nodeIndex<Nlinks );
      return source[nodeIndex];
    }
    assert( nodeIndex >=0 && nodeIndex<Nnodes );
    assert( indexed );
    assert( i>=0 );
    assert( i<indegree[nodeIndex] );
    return source[destsorted[firstB[nodeIndex]+i]];
  }
  
  inline long Network::LinkDest( long nodeIndex, long i=-1 ) {
    if( i==-1 ) {
      assert( nodeIndex<Nlinks );
      return dest[nodeIndex];
    }
    assert( nodeIndex >=0 && nodeIndex<Nnodes );
    assert( indexed );
    assert( i>=0 );
    assert( i<outdegree[nodeIndex] );
    return dest[firstA[nodeIndex]+i];
  }

	
  inline void Network::SetLinkLength(long linkIndex, float x) {
    linklength[linkIndex] = x;
  }

  inline float Network::GetLinkLength(long linkIndex) {
    return linklength[linkIndex];
  }

	inline int Network::Disabled(long linkIndex)
		{
		return disabled[linkIndex];
		}
	inline void Network::SetDisabled(long linkIndex, int state)
		{
		disabled[linkIndex] = state;
		}

};

#endif



