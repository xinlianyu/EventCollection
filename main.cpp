/*
  Main program for policy based DTA
  Song Gao
  May 2004
  Revised by Jing Ding 2012
  Revised by Xinlian Yu 2014
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "main.h"
#include "algs.h"
#include "parser.h"
#include "weights.h"
#include "od.h"

#define paramfile "paralib.dat"

#ifndef lowerbound 
#define lowerbound 1e-5
#endif

#define NUMMODES 31
#define STARTMODE 0
typedef struct {
  //int alg;
  int id;
  int argsreq;
  char *name;
  char *syntax;
} Mode;

static Mode mode[NUMMODES] = {

  {-1, -1, "\nRouting Policy Based Stochastic Dynamic Traffic Assignment", "" },
  {-1, -1, "Song Gao, MIT, 2001-2004\n", ""},
  {-1, -1, "Usage:", ""},

  {-1, -1, "\n***Utilities***", ""},
  {0, 0, "clean", ""},
  {1, 1, "save", "#_DTA_iteration"},
  {2, 1, "resume", "#_DTA_iteration"},
  {3, 0, "avg", ""},
  {4, 2, "breakincsample", "breaking_size breaking_counter"},
  {5, 0, "formatpathshare", ""},

  {-1, -1, "\n***DTA Initialization and Postprocessing Components***", ""},
  {6, 0, "sampleinc", ""},
  {7, 3, "initloader", "#_DTA_iteration incident_counter demand_counter"},
  {8, 2, "initloader_demand", "#_DTA_iteration sample_counter"},
  {9, 0, "postDynaMIT", ""},
  {10, 3, "postloader", "#_DTA_iteration incident_counter demand_counter"},
  {11, 0, "sampleod", ""},
  {12, 1, "gap", "#_DTA_iteration"},
  {13, 2, "moe_policy", "#_DTA_iteration #_Policies"},
  {14, 1, "moe_path", "#_DTA_iteration"},
  {15, 1, "moe", "#_DTA_iteration"},

  {-1, -1, "\n***DTA Major Components***", ""},
  {16, 1, "routing_policy", "#_DTA_iteration"},
  {17, 1, "routing_CE", "#_DTA_iteration"},
  {18, 1, "routing_enroute_path", "#_DTA_iteration"},
  {19, 1, "routing", "#_DTA_iteration"},
  {20, 3, "choice_policy", "#_DTA_iteration MSA_Counter #_Policies"},
  {21, 2, "choice_path", "#_DTA_iteration MSA_Counter"},
  {22, 3, "choice", "#_DTA_iteration MSA_Counter #_Policies"},
  {23, 3, "routing_and_choice", "#_DTA_iteration MSA_Counter #_Policies"},
  {24, 3, "iterative_loader", "#_DTA_iteration loadingcounter #_Policies"},
};


void help( void ) {

  int x;
  for( x=STARTMODE; x<NUMMODES; x++ ) 
    printf( "%s%s %s\n", mode[x].argsreq==-1 ? "" : "dta  ", mode[x].name, mode[x].syntax );
  printf( "\n" );
}


int main( int argc, char *argv[] ) 
{

  	int x, firstarg=1;

  	if( argc <= 1 ) 
  		{
    	help();
    	exit(0);
  		}
  
	

  //read parameter files, some initializations here

	char netfile[256], demandfile[256], linktimefile[256], enroutepathpolicyfile[256], eventfile[256], policyfile[256], policyfile1[256], enroutepathpolicyflowfile[256], policyflowfile[256], pathfile[256], pathflowfile[256], basepathflowfile[256], pathpathflowfile[256], pathflownormfile[256], incsamplefile[256], dynamitincfile[256], trajfile[256], tmplinkfile[256], gapfile[256], moefile[256], basemoefile[256], pathmoefile[256], enroutepathmoefile[256], policymoefile[256], odsamplefile[256], odtimefile[256], ttdistfile[256], basettdistfile[256], pathttdistfile[256], enroutepathttdistfile[256], policyttdistfile[256], specialpolicyflowfile_enroutepath[256], specialpolicyflowfile_policy[256], pathsharefile[256];

  	int num_dest, num_timeperiods, lookbacklength, num_links, stepsize=60, samplesize, odsamplesize, choicesetsize=7, model=0, flag=0,num_inc, assign_start, assign_end, odstart, odend, odinterval, numods;
  	float alpha=1.0, beta=0.0;

  	int *segIndex, *start_a, *start_b, *duration_a, *duration_b;
  	float incrate, baseratio, pathratio, enroutepathratio, policyratio, *severity_a, *severity_b, *incscale, *meanod, *pho, *errorstd;

  	int int_temp;

  	Network *N;
  	EdgeWeights *w;
  	//node_info *node;
	link_info *link;
  	ODDemand *ods;
  
  	char s[256], s1[256], s2[256], s3[256];
  	char cmd[256];

	char *paramfile_new = argv[3];

  	FILE *fp;
  	fp = fopen(paramfile_new, "r");
  	if (!fp) 
  		{
    	printf("Error: Cannot open file %s!\n", paramfile_new);
    	exit(1);
  		}

  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Network File] : %s", netfile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Demand File] : %s", demandfile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Link Travel Time File] : %s", linktimefile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) )
    if( sscanf( s, " [Enroutepath Policy File] : %s", enroutepathpolicyfile ) == 1 ) 
      break;
  while (fgets(s, 255, fp))
	  if (sscanf(s, " [Event File] : %s", eventfile) == 1)
		  break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Policy Policy File] : %s", policyfile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Policy Policy File 1] : %s", policyfile1 ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Enroutepath Policy Flow File] : %s", enroutepathpolicyflowfile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Policy Policy Flow File] : %s", policyflowfile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Special Policy Flow File (Enroutepath)] : %s", specialpolicyflowfile_enroutepath ) == 1 ) 
      break; //for simple paths
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Special Policy Flow File (Policy)] : %s", specialpolicyflowfile_policy ) == 1 ) 
      break; //for simple paths
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Path File] : %s", pathfile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Working Path Flow File] : %s", pathflowfile ) == 1 ) 
      break; 
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Base Path Flow File] : %s", basepathflowfile ) == 1 ) 
      break; 
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Path Path Flow File] : %s", pathpathflowfile ) == 1 ) 
      break; 
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Path Flow Norm File] : %s", pathflownormfile ) == 1 ) 
      break; 
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Incident Sample File] : %s", incsamplefile ) == 1 ) 
      break;  
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [DynaMIT Incident File] : %s", dynamitincfile ) == 1 ) 
      break; 
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [OD Sample File] : %s", odsamplefile ) == 1 ) 
      break;  
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Trajectory File] : %s", trajfile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Temp Link Time File] : %s", tmplinkfile ) == 1 ) 
      break;  
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Gap File] : %s", gapfile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [MOE File] : %s", moefile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [OD Time File (by vehicle)] : %s", odtimefile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Travel Time Distribution File] : %s", ttdistfile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Travel Time Distribution File (Base)] : %s", basettdistfile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Travel Time Distribution File (Path)] : %s", pathttdistfile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Travel Time Distribution File (Enroutepath)] : %s", enroutepathttdistfile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Travel Time Distribution File (Policy)] : %s", policyttdistfile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Path Share File] : %s", pathsharefile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Base MOE File] : %s", basemoefile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Path MOE File] : %s", pathmoefile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [En Route Path MOE File] : %s", enroutepathmoefile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Policy MOE File] : %s", policymoefile ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Num Of Time Periods] : %d", &num_timeperiods ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Assignment Start Time] : %d", &assign_start ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Assignment End Time] : %d", &assign_end ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Step Size] : %d", &stepsize ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Look Back Length] : %d", &lookbacklength ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Num Links] : %d", &num_links ) == 1 ) 
      break;  
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Sample Size] : %d", &samplesize ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [OD Sample Size] : %d", &odsamplesize ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Choice Set Size] : %d", &choicesetsize ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Model (Deterministic=0, RUM=1)] : %d", &model ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [RUM for Deterministic (No=0, Yes=1)] : %d", &flag ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) )
    if( sscanf( s, " [RUM Alpha] : %f", &alpha ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [RUM Beta] : %f", &beta ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Base Penetration Ratio] : %f", &baseratio ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Path Penetration Ratio] : %f", &pathratio ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Enroutepath Penetration Ratio] : %f", &enroutepathratio ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Policy Penetration Ratio] : %f", &policyratio ) == 1 ) 
      break;

  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [OD Start Time (second)] : %d", &odstart ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [OD End Time (second)] : %d", &odend ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [OD Interval (second)] : %d", &odinterval ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Number of OD Pairs] : %d", &numods ) == 1 ) 
      break;

  meanod = new float[numods];
  pho = new float[numods];
  errorstd = new float[numods];

  for (int i=0;i<numods;i++) {
    while( fgets( s, 255, fp ) ) 
      if( sscanf( s, " [Mean Scale] : %f", &meanod[i] ) == 1 ) 
      break;
    while( fgets( s, 255, fp ) ) 
      if( sscanf( s, " [Pho for AR(1)] : %f", &pho[i] ) == 1 ) 
	break;
    while( fgets( s, 255, fp ) ) 
      if( sscanf( s, " [Standard Error for AR(1)] : %f", &errorstd[i] ) == 1 ) 
	break;
  }

  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Number of Incidents] : %d", &num_inc ) == 1 ) 
      break;
  while( fgets( s, 255, fp ) ) 
    if( sscanf( s, " [Incident Rate (/day/mile)] : %f", &incrate ) == 1 ) 
      break;

  incscale = new float[num_inc];
  segIndex = new int[num_inc];
  start_a = new int[num_inc];
  start_b = new int[num_inc];
  duration_a = new int[num_inc];
  duration_b = new int[num_inc];
  severity_a = new float[num_inc];
  severity_b = new float[num_inc];

  for (x=0;x<num_links;x++) {
    while( fgets( s, 255, fp ) )
      if( sscanf( s, " [Incident Segment Index] : %d", &segIndex[x]) == 1 ) 
      break;
    while( fgets( s, 255, fp ) )
      if( sscanf( s, " [Incident Scale] : %f", &incscale[x]) == 1 ) 
      break;
    while (fgets( s, 255, fp ) )
      if ( sscanf( s, " [Incident Start Time Range] : %d %d", &start_a[x], &start_b[x] ) == 2)      
	break; 
    while (fgets( s, 255, fp ) )
      if ( sscanf( s, " [Incident Duration Range] : %d %d", &duration_a[x], &duration_b[x] ) == 2) 
	break;
    while (fgets( s, 255, fp ) )
      if ( sscanf( s, " [Incident Severity Range] : %f %f", &severity_a[x], &severity_b[x] ) == 2) 
	break; 
  }
  
  fclose( fp );
 
	odpairs dynamic_ods;
   
	int policy = 0;
  /* Branch to appropriate subfunction */
  	for (x = STARTMODE; x < NUMMODES; x++) 
    	if( mode[x].argsreq != -1 && strcmp( argv[firstarg], mode[x].name ) == 0 ) 
    		{
      		if( argc-firstarg <= mode[x].argsreq ) 
      			{
				printf( "\nUsage: dta %s %s\n\n", mode[x].name, mode[x].syntax );
				exit(0);
      			}  

			if (mode[x].id > 15) //???
				{
				//N = readnetwork_mitsim(netfile);
                N = readnetwork_simple(netfile);  


				//N->DebugOut();
				//return 0;
				//N->buildpathtable(pathfile); //not needed for routing
				//strcat(linktimefile, argv[firstarg+1]);
				//w = readlinktimes(N, linktimefile, num_timeperiods, stepsize, samplesize);
				w = ReadPiecewiseLinkTimes(N, linktimefile, num_timeperiods, stepsize);


				//w = ReadStepLinkTimes(N, linktimefile, num_timeperiods, stepsize);

				dynamic_ods = ReadDynamicODPairs(N, demandfile);
				
				/*
				node = new node_info[N->NumNodes()];
				for (int n = 0; n < N->NumNodes(); n++) 
					node[n].init(w->levels(), w->numreal());

				//ods = readdemand(demandfile, w->levels(), stepsize, N);
		      	}
			    */

			    link = new link_info[N->NumLinks()];//????
				for (int n = 0; n < N->NumLinks(); n++) 
					link[n].init(w->levels(), w->numreal());

				//ods = readdemand(demandfile, w->levels(), stepsize, N);
		      	}
     
            switch (mode[x].id) 
      			{
				case 16:
					//NOI_StepTT(N, w, dynamic_ods, node, policyfile, assign_start, 0, &policy);
					for (odpairs::iterator it = dynamic_ods.begin(); it != dynamic_ods.end(); it++)
					{
						//CDPI(N, w, (*it).first, node, policyfile, assign_start,(*it).second, &policy);
						CDPI(N, w, (*it).first, link, eventfile, policyfile, policyfile1, assign_start,(*it).second, &policy);//(*it).second : destination?_Yu




					}
  		            
					break;
				default:
					break;
				}  	
			exit(0);
    		}
/*
      case 0:

	  strcpy(cmd, "rm ");
	  strcat(cmd, linktimefile);
	  strcat(cmd, "*");
	  system(cmd);

	if (fp=fopen(policyfile, "r")){
	  strcpy(cmd, "rm ");
	  strcat(cmd, policyfile);
	  strcat(cmd, "*");
	  system(cmd);
	}
	if (fp=fopen(enroutepathpolicyfile, "r")){
	  strcpy(cmd, "rm ");
	  strcat(cmd, enroutepathpolicyfile);
	  strcat(cmd, "*");
	  system(cmd);
	}
	if (fp=fopen(policyflowfile, "r")){
	  strcpy(cmd, "rm ");
	  strcat(cmd, policyflowfile);
	  strcat(cmd, "*");
	  system(cmd);
	}
	if (fp=fopen(enroutepathpolicyflowfile, "r")){
	  strcpy(cmd, "rm ");
	  strcat(cmd, enroutepathpolicyflowfile);
	  strcat(cmd, "*");
	  system(cmd);
	}
	if (fp=fopen(specialpolicyflowfile_policy, "r")){
	  strcpy(cmd, "rm ");
	  strcat(cmd, specialpolicyflowfile_policy);
	  strcat(cmd, "*");
	  system(cmd);
	}
	if (fp=fopen(specialpolicyflowfile_enroutepath, "r")){
	  strcpy(cmd, "rm ");
	  strcat(cmd, specialpolicyflowfile_enroutepath);
	  strcat(cmd, "*");
	  system(cmd);
	}
	if (fp=fopen(tmplinkfile, "r")){
	  strcpy(cmd, "rm ");
	  system(strcat(cmd, tmplinkfile));
	}
	if (fp=fopen(ttdistfile, "r")){
	  strcpy(cmd, "rm ");
	  strcat(cmd, ttdistfile);
	  //strcat(cmd, "*");
	  system(cmd);
	}
	if (fp=fopen(basettdistfile, "r")){
	  strcpy(cmd, "rm ");
	  system(strcat(cmd, basettdistfile));
	}
	if (fp=fopen(pathttdistfile, "r")){
	  strcpy(cmd, "rm ");
	  system(strcat(cmd, pathttdistfile));
	}
	if (fp=fopen(enroutepathttdistfile, "r")){
	  strcpy(cmd, "rm ");
	  system(strcat(cmd, enroutepathttdistfile));
	}
	if (fp=fopen(policyttdistfile, "r")){
	  strcpy(cmd, "rm ");
	  system(strcat(cmd, policyttdistfile));
	}
	if (fp=fopen(pathsharefile, "r")){
	  strcpy(cmd, "rm ");
	  system(strcat(cmd, pathsharefile));
	}
	
	if (fp=fopen(basepathflowfile, "r")){
	  strcpy(cmd, "rm ");
	  strcat(cmd, basepathflowfile);
	  strcat(cmd, "*");
	  system(cmd);
	}
	if (fp=fopen(pathpathflowfile, "r")){
	  strcpy(cmd, "rm ");
	  strcat(cmd, pathpathflowfile);
	  strcat(cmd, "*");
	  system(cmd);
	}
	if (fp=fopen(pathflownormfile, "r")){
	  strcpy(cmd, "rm ");
	  system(strcat(cmd, pathflownormfile));
	}
	if (fp=fopen(gapfile, "r")){
	  strcpy(cmd, "rm ");
	  strcat(cmd, gapfile);
	  //strcat(cmd, "*");
	  system(cmd);
	}
	if (fp=fopen(moefile, "r")){
	  strcpy(cmd, "rm ");
	  strcat(cmd, moefile);
	  //strcat(cmd, "*");
	  system(cmd);
	}
	if (fp=fopen(basemoefile, "r")){
	  strcpy(cmd, "rm ");
	  system(strcat(cmd, basemoefile));
	}
	if (fp=fopen(pathmoefile, "r")){
	  strcpy(cmd, "rm ");
	  system(strcat(cmd, pathmoefile));
	}
	if (fp=fopen(enroutepathmoefile, "r")){
	  strcpy(cmd, "rm ");
	  system(strcat(cmd, enroutepathmoefile));
	}
	if (fp=fopen(policymoefile, "r")){
	  strcpy(cmd, "rm ");
	  system(strcat(cmd, policymoefile));
	}
	
	break;
      case 1:
//4_10	 if (fp=fopen(moefile, "r")){
//4_10	  strcpy(cmd, "cp ");
//4_10	  strcpy(s, moefile);
//4_10	  strcpy(s1, moefile);
//4_10	  strcat(s1, argv[firstarg+1]);
//4_10	  strcat(s, " ");
//4_10	  strcat(s, s1);
//4_10	  system(strcat(cmd, s));
//4_10	  }
//4_10	if (fp=fopen(gapfile, "r")){
//4_10	  strcpy(cmd, "cp ");
//4_10	  strcpy(s, gapfile);
//4_10	  strcpy(s1, gapfile);
//4_10	  strcat(s1, argv[firstarg+1]);
//4_10	  strcat(s, " ");
//4_10	  strcat(s, s1);
//4_10	  system(strcat(cmd, s));
//4_10	  }
	if (policyratio>lowerbound) {
	if (fp=fopen(policyfile, "r")){
	  strcpy(cmd, "cp ");
	  strcpy(s, policyfile);
	  strcpy(s1, policyfile);
	  //strcat(s1, argv[firstarg+1]);
	  strcat(s1, "-old");
	  strcat(s, " ");
	  strcat(s, s1);
	  system(strcat(cmd, s));
	}
	
	if (fp=fopen(policyflowfile, "r")){
	  strcpy(cmd, "cp ");
	  strcpy(s, policyflowfile);
	  strcpy(s1, policyflowfile);
	  strcat(s1, argv[firstarg+1]);
	  strcat(s, " ");
	  strcat(s, s1);
	  system(strcat(cmd, s));
	}
	if (fp=fopen(specialpolicyflowfile_policy, "r")){
	  strcpy(cmd, "cp ");
	  strcpy(s, specialpolicyflowfile_policy);
	  strcpy(s1, specialpolicyflowfile_policy);
	  strcat(s1, argv[firstarg+1]);
	  strcat(s, " ");
	  strcat(s, s1);
	  system(strcat(cmd, s));
	}
	}
	if (enroutepathratio>lowerbound) {
	if (fp=fopen(enroutepathpolicyfile, "r")){
	  strcpy(cmd, "cp ");
	  strcpy(s, enroutepathpolicyfile);
	  strcpy(s1, enroutepathpolicyfile);
	  //strcat(s1, argv[firstarg+1]);
	  strcat(s1, "-old");
	  strcat(s, " ");
	  strcat(s, s1);
	  system(strcat(cmd, s));
	}
	if (fp=fopen(enroutepathpolicyflowfile, "r")){
	  strcpy(cmd, "cp ");
	  strcpy(s, enroutepathpolicyflowfile);
	  strcpy(s1, enroutepathpolicyflowfile);
	  strcat(s1, argv[firstarg+1]);
	  strcat(s, " ");
	  strcat(s, s1);
	  system(strcat(cmd, s));
	}
	if (fp=fopen(specialpolicyflowfile_enroutepath, "r")){
	  strcpy(cmd, "cp ");
	  strcpy(s, specialpolicyflowfile_enroutepath);
	  strcpy(s1, specialpolicyflowfile_enroutepath);
	  strcat(s1, argv[firstarg+1]);
	  strcat(s, " ");
	  strcat(s, s1);
	  system(strcat(cmd, s));
	}
	}
	if (pathratio>lowerbound) {
	if (fp=fopen(pathpathflowfile, "r")){
	  strcpy(cmd, "cp ");
	  strcpy(s, pathpathflowfile);
	  strcpy(s1, pathpathflowfile);
	  strcat(s1, argv[firstarg+1]);
	  strcat(s, " ");
	  strcat(s, s1);
	  system(strcat(cmd, s));
	}
	}
	break;
      case 2:
//4_10	 if (fp=fopen(moefile, "r")){
//4_10	  strcpy(cmd, "cp ");
//4_10	  strcpy(s, moefile);
//4_10	  strcat(s, argv[firstarg+1]);
//4_10	  strcat(s, " ");
//4_10	  strcat(s, moefile);
//4_10	  system(strcat(cmd, s));
//4_10	  }
//4_10	 if (fp=fopen(gapfile, "r")){
//4_10	  strcpy(cmd, "cp ");
//4_10	  strcpy(s, gapfile);
//4_10	  strcat(s, argv[firstarg+1]);
//4_10	  strcat(s, " ");
//4_10	  strcat(s, gapfile);
//4_10	  system(strcat(cmd, s));
//4_10	  }
	if (policyratio>lowerbound) {
	if (fp=fopen(policyfile, "r")){
	  strcpy(cmd, "cp ");
	  strcpy(s, policyfile);
	  //strcat(s, argv[firstarg+1]);
	  strcat(s, "-old");
	  strcat(s, " ");
	  strcat(s, policyfile);
	  system(strcat(cmd, s));
	}
	if (fp=fopen(policyflowfile, "r")){
	  strcpy(cmd, "cp ");
	  strcpy(s, policyflowfile);
	  strcat(s, argv[firstarg+1]);
	  strcat(s, " ");
	  strcat(s, policyflowfile);
	  system(strcat(cmd, s));
	}
	if (fp=fopen(specialpolicyflowfile_policy, "r")){
	  strcpy(cmd, "cp ");
	  strcpy(s, specialpolicyflowfile_policy);
	  strcat(s, argv[firstarg+1]);
	  strcat(s, " ");
	  strcat(s, specialpolicyflowfile_policy);
	  system(strcat(cmd, s));
	}
	}
	if (enroutepathratio>lowerbound) {
	if (fp=fopen(enroutepathpolicyfile, "r")){
	  strcpy(cmd, "cp ");
	  strcpy(s, enroutepathpolicyfile);
	  //strcat(s, argv[firstarg+1]);
	  strcat(s, "-old");
	  strcat(s, " ");
	  strcat(s, enroutepathpolicyfile);
	  system(strcat(cmd, s));
	}
	if (fp=fopen(enroutepathpolicyflowfile, "r")){
	  strcpy(cmd, "cp ");
	  strcpy(s, enroutepathpolicyflowfile);
	  strcat(s, argv[firstarg+1]);
	  strcat(s, " ");
	  strcat(s, enroutepathpolicyflowfile);
	  system(strcat(cmd, s));
	}
	if (fp=fopen(specialpolicyflowfile_enroutepath, "r")){
	  strcpy(cmd, "cp ");
	  strcpy(s, specialpolicyflowfile_enroutepath);
	  strcat(s, argv[firstarg+1]);
	  strcat(s, " ");
	  strcat(s, specialpolicyflowfile_enroutepath);
	  system(strcat(cmd, s));
	}
	}
	if (pathratio>lowerbound) {
	if (fp=fopen(pathpathflowfile, "r")){
	  strcpy(cmd, "cp ");
	  strcpy(s, pathpathflowfile);
	  strcat(s, argv[firstarg+1]);
	  strcat(s, " ");
	  strcat(s, pathpathflowfile);
	  system(strcat(cmd, s));
	}
	}
	break;
      case 3:
	avg(tmplinkfile, num_timeperiods, num_links);
	break;
      case 4:
	breakincsample(samplesize, stringToInteger(argv[firstarg+1]), incsamplefile, argv[firstarg+2]);
	break;
      case 5:
	N = readnetwork_mitsim(netfile);
	N->buildpathtable(pathfile);
	ods = readdemand(demandfile, num_timeperiods, stepsize, N);
	formatpathshare(pathflowfile, pathsharefile, num_timeperiods, num_links, N, ods);
	break;
      case 6:
	sampleincident(incsamplefile, samplesize, num_inc, incrate, incscale, segIndex, start_a, start_b, duration_a, duration_b, severity_a, severity_b, stepsize);
	break;
      case 7:
	N = readnetwork_mitsim(netfile);
	ods = readdemand(demandfile, num_timeperiods, stepsize, N);
	initloader(argv[firstarg+1], stringToInteger(argv[firstarg+2]), stringToInteger(argv[firstarg+3]), odsamplesize, linktimefile, tmplinkfile, num_timeperiods, num_links, incsamplefile, dynamitincfile, odsamplefile, demandfile, odstart, odend, odinterval, stepsize, ods, numods);
	break;
      case 8:
	N = readnetwork_mitsim(netfile);
	ods = readdemand(demandfile, num_timeperiods, stepsize, N);
	initloader_demand(argv[firstarg+1], stringToInteger(argv[firstarg+2]), linktimefile, tmplinkfile, num_timeperiods, num_links, odsamplefile, demandfile, odstart, odend, odinterval, stepsize, ods);
	break;
      case 9:
	N = readnetwork_mitsim(netfile);
	ods = readdemand(demandfile, num_timeperiods, stepsize, N);
	postDynaMIT(trajfile, tmplinkfile, odtimefile, ods->getNumVeh(), num_timeperiods, num_links, ods->getStart(), stepsize);
	break;
      case 10:
	N = readnetwork_mitsim(netfile);
	N->buildpathtable(pathfile);
	ods = readdemand(demandfile, num_timeperiods, stepsize, N);
	postloader(tmplinkfile, linktimefile, pathflowfile, ttdistfile, argv[firstarg+1], integerToString(stringToInteger(argv[firstarg+1])+1), stringToInteger(argv[firstarg+2]), stringToInteger(argv[firstarg+3]), odsamplesize, num_timeperiods, num_links, N, ods, assign_start, assign_end, stepsize);
	break;
      case 11:
	sampleARod(odsamplefile, odsamplesize, (odend-odstart)/odinterval, numods, meanod, pho, errorstd);
	break;
      case 12:
	gap(stringToInteger(argv[firstarg+1]),argv[firstarg+1], integerToString(stringToInteger(argv[firstarg+1])+1), linktimefile, samplesize, num_timeperiods, num_links, gapfile);
	break;
      case 13:
	N = readnetwork_mitsim(netfile);
	strcat(linktimefile, integerToString(stringToInteger(argv[firstarg+1])+1)); //work with the latest linktime file
	w = readlinktimes(N, linktimefile, num_timeperiods, stepsize, samplesize);
	ods = readdemand(demandfile, w->levels(), stepsize, N);
	N->buildpathtable(pathfile);
	moe(N, w, policyfile, policyflowfile, moefile, stringToInteger(argv[firstarg+2]), choicesetsize, ods->getNumDest(), ods, assign_start, assign_end, lookbacklength);
	break;
      case 14:
	N = readnetwork_mitsim(netfile);
	strcat(linktimefile, integerToString(stringToInteger(argv[firstarg+1])+1)); //work with the latest linktime file
	w = readlinktimes(N, linktimefile, num_timeperiods, stepsize, samplesize);
	ods = readdemand(demandfile, w->levels(), stepsize, N);
	N->buildpathtable(pathfile);
	moepath(N, w, pathpathflowfile, moefile, ods, assign_start, assign_end);
	break;
      case 15:
	//N = readnetwork_mitsim(netfile);
	//strcat(linktimefile, integerToString(stringToInteger(argv[firstarg+1])+1)); //work with the latest linktime file
	//w = readlinktimes(N, linktimefile, num_timeperiods, stepsize, samplesize);
	//ods = readdemand(demandfile, w->levels(), stepsize, N);	
	//N->buildpathtable(pathfile);
	//moe_all(N, w, policyfile, policyflowfile, enroutepathpolicyfile, enroutepathpolicyflowfile, specialpolicyflowfile_policy, specialpolicyflowfile_enroutepath, basepathflowfile, pathpathflowfile, moefile, stringToInteger(argv[firstarg+2]), choicesetsize, ods->getNumDest(), ods, assign_start, assign_end, lookbacklength, baseratio, pathratio, enroutepathratio, policyratio, basemoefile, pathmoefile, enroutepathmoefile, policymoefile, ttdistfile, basettdistfile, pathttdistfile, enroutepathttdistfile, policyttdistfile);
	moe_easy(samplesize, num_timeperiods, numods, ttdistfile, argv[firstarg+1], moefile);
	break;
      case 16:
	w->write_weights(policyfile, N);
	CDPI( N, w, ods->getDests(), node, policyfile);
	break;
      case 17:
	CE( N, w, ods->getDests(), node, policyfile, 0, 0);
	break;
      case 18:
	w->write_weights(enroutepathpolicyfile, N);
	OLF( N, w, ods->getDests(), node, enroutepathpolicyfile);
	break;
      case 19:
	if (policyratio>lowerbound) {
	  w->write_weights(policyfile, N);
	  CDPI( N, w, ods->getDests(), node, policyfile);
	}
	if (enroutepathratio>lowerbound) {
	  w->write_weights(enroutepathpolicyfile, N);
	  N->buildpathtable(pathfile);
	  OLF( N, w, ods->getDests(), node, enroutepathpolicyfile);
	}
	break;
      case 20:
	choicemodel(policyfile, policyflowfile, specialpolicyflowfile_policy, stringToInteger(argv[firstarg+2]), stringToInteger(argv[firstarg+3]), ods->getNumDest(), N, w, choicesetsize, model, flag, alpha, beta, ods, 0, num_timeperiods, lookbacklength);
	break;
      case 21:
	pathchoice(pathpathflowfile, N, w, ods, stringToInteger(argv[firstarg+2]), alpha, beta, 0, num_timeperiods);
	break;
      case 22:
	if (enroutepathratio>lowerbound)
	  choicemodel(enroutepathpolicyfile, enroutepathpolicyflowfile,  specialpolicyflowfile_enroutepath, stringToInteger(argv[firstarg+2]), stringToInteger(argv[firstarg+3]), ods->getNumDest(), N, w, choicesetsize, model, flag, alpha, beta, ods, 0, num_timeperiods, lookbacklength);
	if (policyratio>lowerbound)
	  choicemodel(policyfile, policyflowfile,  specialpolicyflowfile_policy, stringToInteger(argv[firstarg+2]), stringToInteger(argv[firstarg+3]), ods->getNumDest(), N, w, choicesetsize, model, flag, alpha, beta, ods, 0, num_timeperiods, lookbacklength);
	if (pathratio>lowerbound) {
	  pathchoice(pathpathflowfile, N, w, ods, stringToInteger(argv[firstarg+2]), alpha, beta, 0, num_timeperiods);
	}
	break;
      case 23:
	if (policyratio>lowerbound) {
	  w->write_weights(policyfile, N);
	  CDPI( N, w, ods->getDests(), node, policyfile);
	}
	if (enroutepathratio>lowerbound) {
	  w->write_weights(enroutepathpolicyfile, N);
	  OLF( N, w, ods->getDests(), node, enroutepathpolicyfile);
	}
	if (enroutepathratio>lowerbound)
	  choicemodel(enroutepathpolicyfile, enroutepathpolicyflowfile, specialpolicyflowfile_enroutepath, stringToInteger(argv[firstarg+2]),stringToInteger(argv[firstarg+3]), ods->getNumDest(), N, w, choicesetsize, model, flag, alpha, beta, ods, 0, num_timeperiods, lookbacklength);
	if (policyratio>lowerbound)
	  choicemodel(policyfile, policyflowfile, specialpolicyflowfile_policy, stringToInteger(argv[firstarg+2]), stringToInteger(argv[firstarg+3]), ods->getNumDest(), N, w, choicesetsize, model, flag, alpha, beta, ods, 0, num_timeperiods, lookbacklength);
	if (pathratio>lowerbound) {
	  pathchoice(pathpathflowfile, N, w, ods, stringToInteger(argv[firstarg+2]), alpha, beta, 0, num_timeperiods);
	}
	break;
      case 24:
	loader(stringToInteger(argv[firstarg+2]), policyfile, enroutepathpolicyfile, tmplinkfile, policyflowfile, enroutepathpolicyflowfile, specialpolicyflowfile_policy, specialpolicyflowfile_enroutepath, basepathflowfile, pathpathflowfile, ods, pathflowfile, pathflownormfile, ods->getNumDest(), stringToInteger(argv[firstarg+3]), choicesetsize, N, w, 0, num_timeperiods, lookbacklength, baseratio, pathratio, enroutepathratio, policyratio);
	break;
      }
*/
  
  	help();

  	return 0;
}
