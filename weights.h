/*
  Edge Weight Classes
  Brian Dean
  1997
  Adapted by Song Gao for stochastic dependent link travel times
  April 29, 2004
  Song Gao, Added NOI related variables, Feburary 2007
*/
#pragma comment(lib,"ws2_32.lib")
#ifndef EDGEWEIGHTS
#define EDGEWEIGHTS

#define INT_INFINITY   (int)2147000000
#define FLOAT_INFINITY (float)1.0e38
#define DOUBLE_MAX		1.7e308

#include "network.h"

//use a piece-wise linear function to describe dynamic link travel times
typedef pair<long, float> bp;	   //breaking point
typedef vector<bp> bps; //breaking points 

class EdgeWeights          
{
private:
  	long 	_levels;   //number of time periods
  	long 	_numreal;  //number of joint realizations (R)
  	long	_agglength; //aggregation length of link travel times (i.e. length of a time interval) in seconds
	long	_numdefaults; //number of default breaking point
	long	_numtimeslices;  //number of time slices for dynamic link travel time step functions

public:
  	float 	**realization; //joint realizations of link travel times
  	float 	*prob; //probabilities of each joint realization
  	//float 	**expectedtt; //expected travel time for each time and each link
	//float	***marginaltt; //marginal travel time realizations for each time and link
	//float	***marginalprob; //marginal probabilities
	//int		**nummarginalprob; //number of marginal support points for each time and link
	long	*defaults; //the default breaking points entry time
	bps		**BP;   //breaking point pairs for each link and each support point (both deterministic and stochastic links)
	bps		*ExpectedBP; //breaking point pairs for expected travel times of each link
	//short	**det_flag; //flag for each link whether it is determinstic
	short	*det_flag;
 
	long	*slicestart; //entry time of each time slice (the same over all links): start with 12am and end with 11:59pm
	long	**numsupportpoints; //number of support points for each link and each time slice
	long	***supportpoints; //support points for each link and each time slice
	float	***probability;  //probability of each support point for each link at each time slice (#links x #time slices x #support points)
	float	*avgtime; //average travel times
	int ***blocked;

	void 	write_weights(char *filename, Network *N);
  	void 	calexpected(Network *N, long EV_start, long EV_end, long **event, long start_time);
	void	calexpectedBP(Network *N);
	//void	calexpected(Network *N);
	//void	calmarginal(Network *N);
  	float 	*pathtimes(long origIndex, long destIndex, long pathIndex, long t, Network *N, long EV_start, long EV_end, long flag, long **event=NULL);
  	float 	pathtime(long origIndex, long destIndex, long pathIndex, long t, Network *N, long r);

	EdgeWeights(Network *N, long __levels = 1, long __agglength = 1, long __numreal = 1) 
	{
		long i, numlinks = N->NumLinks();

		assert( __levels > 0 );
		_levels = __levels;

		assert( __agglength > 0);
		_agglength = __agglength;

		assert(__numreal > 0);
		_numreal = __numreal;
		//_numdefaults = -1;

		BP = new bps*[numlinks];
		ExpectedBP = new bps[numlinks];

		blocked = new int**[numlinks];
	}
  
  	~EdgeWeights() 
  		{
    	delete[] realization;
    	delete prob;
		delete ExpectedBP;
		delete[] det_flag;
		delete slicestart;
		delete[] numsupportpoints;
		delete[] probability; 
		delete[] BP;
		delete[] blocked;
  		}

	/*
    inline int Blocked(long linkIndex, int t, int r)
    {
		return blocked[linkIndex][t][r];
	}
	inline void SetBlocked(long linkIndex, int t, int r, int state)
	{
		blocked[linkIndex][t][r] = state;
	}
	*/

	inline long find_interval(long startslice, float time)
	{
		long	i;
		for (i = startslice; i < _numtimeslices - 1; i++)
			if (time >= slicestart[i] && time < slicestart[i + 1])  //left inclusive, right exclusive
				break;
		return i;
	}

	inline bp find_interval(bps *x, float time, int *i)
	//*i is the index of the right breaking point; prev is the left breaking point
		{
		bp	prev;

		for (*i = 0; *i < (*x).size() && (*x)[*i].first < time; (*i)++);

		if (*i == 0)
			prev.first = -2; //ray to the left
		else if (*i == (*x).size())
			{
			prev.first = -1; //ray to the right
			*i = *i - 1;
			}
		else
			prev = (*x)[*i - 1];
		return prev;
		}


  	inline float weight_link(long link, int time, Network *_N, bps *x) 
  		{
		float	ratio, starttime, endtime;
		int		i, start, end;
        bp      start_bp, end_bp;

    	assert(link >= 0 && link < _N->NumLinks());

		start_bp = find_interval(x, time, &i);
		start = start_bp.first;
		starttime = start_bp.second;

		end_bp = (*x)[i];
        end = end_bp.first;
        endtime = end_bp.second;
		
		if (start < 0) //static part
			return endtime;

        assert(time >= start && time <= end);

		ratio = (float)(time - start) / (end - start);
    	return starttime * (1 - ratio) + endtime * ratio;
  		}
    
  	inline float weight(long source, long dest, long time, Network *_N, bps *x) 
  		{
		long	link = _N->LinkIndex(source, dest);
    	assert(link != -1 );
    	return weight_link(link, time, _N, x);
  		}

	inline float weight_link(long link, int time, int r, Network *_N)
	{
		int	links = _N->NumLinks();
		assert(link >= 0 && links);
		return (realization[time * links + link][r]);
	}

	inline float weight(long source, long dest, long time, int r, Network *_N)
	{
		long link = _N->LinkIndex(source, dest);
		assert(link != -1);
		return weight_link(link, time, r, _N);
	}

  	inline float weight_in(long dest, long i, long time, Network *_N, bps *x) 
  		{
    	assert( i>=0 && i<_N->InDegree(dest) );
    	return weight_link( _N->LinkInIndex( dest, i ), time, _N, x);
  		}
  
  	inline float weight_out( long source, long i, long time, Network *_N, bps *x) 
  		{
    	assert( i>=0 && i<_N->OutDegree(source) );
    	return weight_link( _N->LinkOutIndex( source, i ), time, _N, x);
  		}
  
  	inline long timehorizon( void ) 
  		{
    	return _levels-1;
  		}
  
  	inline long levels( void ) 
  		{
    	return _levels;
  		}
  
  	inline long numreal(void) 
  		{
    	return _numreal;
  		}
  
  	inline long agglength(void) 
  		{
    	return _agglength;
  		}
  
	inline long numdefaults(void)
		{
		return _numdefaults;
		}

	inline long numtimeslices(void)
	{
		return _numtimeslices;
	}

  	
    inline void setweight(long source, long dest, long level, long real, float weight, Network *_N ) 
  		{
    	assert( level>=0 && level<_levels );
    	assert( _N->LinkIndex( source, dest ) != -1 );
    	assert(real >=0 && real<_numreal);
    	//_minweight = (weight < _minweight) ? weight : _minweight;
    	//_maxweight = (weight > _maxweight) ? weight : _maxweight;
    	realization[level * _N->NumLinks() + _N->LinkIndex( source, dest )][real] = weight;
  		}

  	inline void setweight_link( long linkIndex, long level, long real, float weight , Network *_N) 
  		{
    	assert( level>=0 && level<_levels );
    	assert( linkIndex>=0 && linkIndex<_N->NumLinks() );
    	assert(real >=0 && real<_numreal);
    	//_minweight = (weight < _minweight) ? weight : _minweight;
    	//_maxweight = (weight > _maxweight) ? weight : _maxweight;
    
    	realization[level * _N->NumLinks() + linkIndex][real] = weight;
  		}
     

  inline void setprob( long real, float probability) {
    assert(real >=0 && real<_numreal);
    assert( probability >0 && probability <= 1);
    prob[real] = probability;
  }

	inline void setlevels(long levels)
		{
		assert(levels > 0);
		_levels = levels;
		}

	inline void setagglength(long agglength)
		{
		assert(agglength > 0);
		_agglength = agglength;
		}

	inline void setnumreal(long numreal)
		{
	    assert(numreal > 0);
		_numreal = numreal;
		}

	inline void setnumdefaults(long numdefaults)
		{
		assert(numdefaults > 0);
		_numdefaults = numdefaults;
		}

	inline void setnumslices(long numtimeslices)
		{
		assert(numtimeslices > 0);
		_numtimeslices = numtimeslices;
		}

};
	
#endif
