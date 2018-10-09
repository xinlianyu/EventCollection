/*
  OD Demand Class
  Song Gao
  May 2004
*/

#pragma comment(lib,"ws2_32.lib")
#ifndef ODDEMAND
#define ODDEMAND

#include <utility>
#include <vector>
#include <map>

using namespace std;

typedef pair< int, int> OD;
typedef vector<OD> ODs;

typedef map<long, vector<pair<long, int>>> odpairs;   //for route choice (input only has orig, dest, departure time, no od flow)


class ODDemand { 
 private:
  int num_timeperiods;
  vector<int> listofdest;
  int start;
  int numveh;

 public:
  ODs *ods;
  vector< int > *demand;

  void find_list_of_dest();
  void printOD();

  ODDemand(int _timeperiods) {
    ods = new ODs[_timeperiods];
    demand = new vector< int >[_timeperiods];
    num_timeperiods = _timeperiods;
  }

  ~ODDemand() {  
  }

  inline void setOD(int orig, int dest, int time, int pos) {
    ods[time][pos].first = orig;
    ods[time][pos].second = dest;
  }

  inline void addOD(int origIndex, int destIndex, int time) { 
 
    OD od(origIndex, destIndex);
    ods[time].push_back(od);
 
  }

  inline void addDemand(int thisdemand, int time) {
    demand[time].push_back(thisdemand);
  }

  inline int getOrig(int time, int pos) {
    return ods[time][pos].first;
  }

  inline int getDest(int time, int pos) {
    return ods[time][pos].second;
  }

  
  inline int getNumDest() {
    return listofdest.size();
  }

  inline int getLevels() {
    return num_timeperiods;
  }

  inline vector<int> getDests() {
    return listofdest;
  }
  
  inline void setStart(int t){
    start = t;
  }

  inline int getStart(){
    return start;
  }

  inline void setNumVeh( int num) {
    numveh = num;
  }

  inline int getNumVeh() {
    return numveh;
  }

};

#endif
