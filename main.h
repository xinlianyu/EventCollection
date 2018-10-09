/*
  Main program for policy based DTA
  Song Gao
  May 2004
  Revised by Jing Ding 2012
*/
#pragma comment(lib,"ws2_32.lib")
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int stringToInteger(char* string) {

  int result=0;
  int i=0;

  do {
    result=10*result+(int)(string[i++])-(int)('0');
  }while(string[i]!=0);

  return result;

}

char *integerToString(int i) {
  //printf("the input interger: %d\n", i);
  char *result = new char[256];
  char *temp = new char[256]; temp[0]='\0';
  char *current = new char[2]; current[1]='\0';
  int a, counter=0;

  int original = i;

  while (i!=0) {
    a = i - 10*(i/10);
    //printf("a: %d\n", a);
    current[0] = (char)(a+(int)('0')); 
    //printf("current: %s\n", current);
    strcat(temp, current);
    //printf("temp: %s\n", temp);
    i/=10;
    counter++;
  }
  
  for (int x=0;x<counter;x++)
    result[x] = temp[counter-1-x];

  result[counter]='\0';

  return result;
}

void avg(char *tmplinkfile, int num_timeperiods, int num_links) {
  
    FILE *fp;
    fp=fopen(tmplinkfile, "rt");
    if  (!fp)
    {
        printf("Error: Cannot open file %s!\n", tmplinkfile);
        exit(1);
    }

    float **avglinktt;
    avglinktt = new float*[num_timeperiods];
    for (int t=0;t<num_timeperiods;t++)
        avglinktt[t] = new float[num_links];

  for (int t=0;t<num_timeperiods;t++)
    for (int j=0;j<num_links;j++)
      avglinktt[t][j] = 0;

  float tmp;

  int counter = 0;
  while(!feof(fp)) {
    counter++; 
    for (int t=0;t<num_timeperiods;t++) {
      for (int j=0;j<num_links;j++) {
	fscanf(fp, "%f ", &tmp);
	//cout << tmp << " ";
	avglinktt[t][j] = avglinktt[t][j]*(counter-1)/counter + tmp/counter;
      }
      //cout << endl;
    }
  }

  fclose(fp);

  //output file
  fp=fopen(tmplinkfile, "wt");
  if(!fp) {
    printf("Error: Cannot open file %s!\n", tmplinkfile);
    exit(1);
  }

  for (int t=0;t<num_timeperiods;t++){
      for (int j=0;j<num_links;j++) 
	fprintf(fp, "%f ", avglinktt[t][j]);
      fprintf(fp, "\n");
  }
  
  fclose(fp);

}

void breakincsample(int numreal, int breaksize, char *incsamplefile, char *d) {
  
  char filename[256];

  FILE *fp1, *fp2;
  fp1=fopen(incsamplefile, "rt");
  if (!fp1) {
    printf("Error: Cannot open file %s!\n", incsamplefile);
    exit(0);
  }

  strcpy(filename, incsamplefile);
  strcat(filename, d);
  
  fp2=fopen(filename, "wt");
  if (!fp2) {
    printf("Error: Cannot open file %s!\n", filename);
    exit(0);
  }

  int startindex=stringToInteger(d)*breaksize;
  int segIndex, start, end;
  float severity;

  char s[256];

  int i=0, int_tmp=0;
  while( fgets( s, 255, fp1 ) && i<breaksize+1) {
    if ( sscanf( s, "[ r = %d ]", &int_tmp ) == 1 ) {
      if (int_tmp >= startindex+1 && int_tmp < startindex+breaksize+1) {
	fprintf(fp2, "\n[ r = %d ]\n", i+1);
	i++;
      }
    }
    else if( sscanf( s, "%d %d %d %f", &segIndex, &start, &end, &severity) == 4 ) {
      if (int_tmp >= startindex+1 && int_tmp < startindex+breaksize+1) 
	fprintf(fp2, "%d %d %d %f\n", segIndex, start, end, severity);
    }
  }
    
  fprintf(fp2, "\n[ r = %d ]\n", breaksize+1);

  fclose(fp1);
  fclose(fp2);
}

