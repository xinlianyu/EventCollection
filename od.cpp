/*
  OD Demand Class
  Song Gao
  May 2004
*/

#include "od.h"
#include <stdio.h>
#include <algorithm>

void ODDemand::find_list_of_dest() {
  //1. scan all the OD pairs in all time periods
  //2. sort them by destination index
  //3. get the distinct destinations. 

  vector <int> dests;
  int currentDest = -10000;

  for (int t=0;t<num_timeperiods;t++) 
    for (ODs::iterator it=ods[t].begin();it!=ods[t].end();it++) 
      dests.push_back((*it).second);
  
  sort(dests.begin(), dests.end());
  
  for (vector<int>::iterator it=dests.begin();it!=dests.end();it++) 
    if (*it != currentDest) {
      listofdest.push_back(*it);
      currentDest = *it;
    }
  
}

void ODDemand::printOD() {
  for (int t=0;t<num_timeperiods;t++) {
    printf("time %d: \n", t);
    for (ODs::iterator it=ods[t].begin();it!=ods[t].end();it++) 
      printf("%d %d \n", (*it).first, (*it).second);
  }
}
