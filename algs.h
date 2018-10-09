/*
  Network Algorithms
  Song Gao
  May 2004
  Revised by Jing Ding 2012
*/
#pragma comment(lib,"ws2_32.lib")
#ifndef ALGS
#define ALGS

#include <vector>
#include "network.h"
#include "weights.h"
#include "od.h"

class node_info
{
public:
  	//for CDPI
  	float 	**cost_label; //at each node, optimal cost is dependent on time t and event collection EV
  	int 	**next_node; 

  	//for CE
//4_10  	float 	*cost_label_ce; //nominal optimal cost computed from CE
//4_10  	int 	*next_node_ce;
//4_10  	float	*adaptive_cost_label_ce; //cost of adaptive CE evaluated over joint distribution of link travel times  

	//for NI
	double 	*cost_label_noi; //nominal optimal cost computed from NOI
	int 	*next_node_noi;
//4_10	float	*true_cost_label_noi; //true cost of NOI evaluated over joint distribution of link times

  	void node_info::init(int _levels, int numreal) 
  		{
		int t, r;
		//should initialize all values to zero
    	cost_label = new float*[_levels];
    	next_node = new int*[_levels];
//4_10
//4_10    	cost_label_ce = new float[_levels];
//4_10    	next_node_ce = new int[_levels];
//4_10    	adaptive_cost_label_ce = new float[_levels];
    	
    	cost_label_noi = new double[_levels];
    	next_node_noi = new int[_levels];
//4_10    	true_cost_label_noi = new float[_levels];

    	for (t = 0; t < _levels; t++)
    		{
			//cost_label_ce[t] = INT_INFINITY;
			//adaptive_cost_label_ce[t] = INT_INFINITY;
			cost_label_noi[t] = INT_INFINITY;
//4_10			true_cost_label_noi[t] = INT_INFINITY;

      		cost_label[t] = new float[numreal];
     		next_node[t] = new int[numreal];
			for (r = 0; r < numreal; r++)
				cost_label[t][r] = INT_INFINITY;
    		}
  		}

  	node_info() {}  //a blank constructor.  the true initialization is in init()

  	~node_info() 
  		{
    	delete[] cost_label;
    	delete[] next_node;
//4_10    	delete[] cost_label_ce;
//4_10    	delete[] next_node_ce;
//4_10    	delete[] adaptive_cost_label_ce;
    	delete[] cost_label_noi;
    	delete[] next_node_noi;
//4_10    	delete[] true_cost_label_noi;
  		}
};

class link_info
{
public:
  	//for CDPI
  	//float 	**cost_label; //at each node, optimal cost is dependent on time t and event collection EV
	double 	**cost_label;
  	int 	**next_node; 

  	//for CE
//4_10  	float 	*cost_label_ce; //nominal optimal cost computed from CE
//4_10  	int 	*next_node_ce;
//4_10  	float	*adaptive_cost_label_ce; //cost of adaptive CE evaluated over joint distribution of link travel times  

	//for NI
	//double 	*cost_label_noi; //nominal optimal cost computed from NOI
	//int 	*next_node_noi;
//4_10	float	*true_cost_label_noi; //true cost of NOI evaluated over joint distribution of link times

	void link_info::init(int _levels, int numreal) 
  		{
		int t, r;
		//should initialize all values to zero
    	//cost_label = new float*[_levels];
		cost_label = new double*[_levels];
    	next_node = new int*[_levels];
//4_10
//4_10    	cost_label_ce = new float[_levels];
//4_10    	next_node_ce = new int[_levels];
//4_10    	adaptive_cost_label_ce = new float[_levels];
    	
    	//cost_label_noi = new double[_levels];
    	//next_node_noi = new int[_levels];
//4_10    	true_cost_label_noi = new float[_levels];

    	for (t = 0; t < _levels; t++)
    		{
			//cost_label_ce[t] = INT_INFINITY;
			//adaptive_cost_label_ce[t] = INT_INFINITY;
			//cost_label_noi[t] = INT_INFINITY;
//4_10			true_cost_label_noi[t] = INT_INFINITY;

      		//cost_label[t] = new float[numreal];
			cost_label[t] = new double[numreal];
     		next_node[t] = new int[numreal];
			for (r = 0; r < numreal; r++)
				cost_label[t][r] = INT_INFINITY;
    		}
  		}

  	link_info() {}  //a blank constructor.  the true initialization is in init()

  	~link_info() 
  		{
    	delete[] cost_label;
    	delete[] next_node;
//4_10    	delete[] cost_label_ce;
//4_10    	delete[] next_node_ce;
//4_10    	delete[] adaptive_cost_label_ce;
    	//delete[] cost_label_noi;
    	//delete[] next_node_noi;
//4_10    	delete[] true_cost_label_noi;
  		}
};



//void CDPI( Network *N, EdgeWeights *w, int dests, node_info *node, char *filename, int start_time, vector<pair<long, int>> orig_dep, int *policy);
void CDPI( Network *N, EdgeWeights *w, int dests, link_info *link, char *eventfile, char *filename, char *filename1, int start_time, vector<pair<long, int>> orig_dep, int *policy);


//void CE( Network *N, EdgeWeights *w, vector<int> dests, node_info *node, char *filename, int start, int flag); //the transform of the data can be done in the network class
//expected travel times are a property of the network
//add start in the argument list, and then CE can be used in OLFCE
//both nominal and true expected cost for each pair (j, t) are calculated

//void OLFCE( Network *N, EdgeWeights *w, vector<int> dests, node_info *node, char *filename);
//the calculation of conditional expected travel times can be implemented as a member function of network class
//call the function only in OLFCE

//void OLF(Network *N, EdgeWeights *w, vector<int> dests, node_info *node, char *policyfile); 
//path enumeration is used to compute minimum expetected travel time path

//void NOI(Network *N, EdgeWeights *w, odpairs ods, node_info *node, char *policyfile, int start_time, int flag, int *policy);
//No-online-information approximation to CDPI
//both nominal and true expected cost for each (j, t) pair are calculated
 
//void NOI_StepTT(Network *N, EdgeWeights *w, odpairs ods, node_info *node, char *policyfile, int start_time, int flag, int *policy);

//void generate_event_collection(Network *N, EdgeWeights *w, int start_time);
//void generate_event_collection(Network *N, EdgeWeights *w);
void generate_event_collection(Network *N, EdgeWeights *w, int dests, vector<pair<long, int>> orig_dep, int **event, int **eventTag, int **realizationToEvent, char* eventfile, char *filename, int start_time);
//void static_shortest_path( Network *N, EdgeWeights *w, int destIndex,  node_info *node, int start_time);
void static_shortest_path( Network *N, EdgeWeights *w, int destIndex,  link_info *link, int start_time);
//void static_shortest_path( Network *N, EdgeWeights *w, int destIndex,  link_info *link);
//void static_shortest_path_noi( Network *N, EdgeWeights *w, int destIndex, node_info *node);
//void write_result_CDPI( Network *N, EdgeWeights *w, int destIndex, int **event, int **eventTag, node_info *node, int **realizationToEvent, char *filename); //write joint-realization based results
//void write_result( Network *N, EdgeWeights *w, int destIndex, int **eventTag, node_info *node, int **realizationToEvent, char *filename); //write joint-realization based results
void write_event(Network *N, EdgeWeights *w, long destIdx, vector<pair<long, int>> orig_dep, int **event, char *filename);
//void write_result_simple( Network *N, EdgeWeights *w, int destIndex, node_info *node, char *filename); //write results based on arrival time only (CE or NOI)
//vector<long> write_route_choice_result(Network *N, EdgeWeights *w, long destIdx, vector<pair<long, int>> orig_dep, node_info *node, int start_time, int policy, char *filename, int op_flag);
//vector<long> write_route_choice_result_new(Network *N, EdgeWeights *w, long destIdx, vector<pair<long, int>> orig_dep, node_info *node, int start_time, int policy, char *filename, int op_flag, int **realizationToEvent);
vector<long> write_route_choice_result_new(Network *N, EdgeWeights *w, long destIdx, vector<pair<long, int>> orig_dep, link_info *link, int start_time, int policy, char *filename, char *filename1, int op_flag, int **realizationToEvent,double runtime, double runtime1,double runtime3);
//vector<long> Policy2Path_NOI(Network *N, EdgeWeights *w, long destIdx, vector<pair<long, int>> orig_dep, node_info *node, int start_time, int policy, char *filename, int op_flag);
//vector<long> Policy2Path_CDPI(Network *N, EdgeWeights *w, long destIdx, vector<pair<long, int>> orig_dep, node_info *node, int start_time, int policy, char *filename, int op_flag);
#endif

