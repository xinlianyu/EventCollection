/*
  Network IO Routines
  Brian Dean
  1997
  Revised by Song Gao May 2004
  Added a function to read OD demand
  Revised the function to read stochastic time-dependent weights
*/

#pragma comment(lib,"ws2_32.lib")
#ifndef PARSER
#define PARSER

#include <stdlib.h>
#include "network.h"
#include "weights.h"
#include "od.h"

Network 	*readnetwork_mitsim( char *filename );
Network 	*readnetwork_simple( char *filename );
EdgeWeights *readlinktimes(Network *N, char *filename, long levels, long agglength, long numreal);
EdgeWeights *ReadPiecewiseLinkTimes(Network *N, char *filename, long levels, long agglength); 
//EdgeWeights *ReadStepLinkTimes(Network *N, char *filename, long levels, long agglength); 
ODDemand 	*readdemand(char *filename, long levels, long stepsize, Network *N);
odpairs 	ReadDynamicODPairs(Network *N, char *filename);

#endif









