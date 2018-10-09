// heapSort Algorithm

#include "heapSort.h"

inline void swap(int b[], int i, int j){ //swap two elements of an array
  int t = b[i];
  b[i] = b[j];
  b[j] = t;
}

//void heapSort(float a[], int b[], int N) { //a stores the keys, b stores the sorted indices
void heapSort(double a[], int b[], int N) {
  makeHeap(a, b, N);

  for (int i=N-1; i>0; i--) { 
    swap (b, i, 0);
    sift (a, b, i-1, 0);
  }
}

//void makeHeap (float a[], int b[], int N) {
void makeHeap (double a[], int b[], int N) {

 
  for (int i=(N-1)/2;i>=0;i--) sift (a, b, N-1, i);

}

//void sift (float a[], int b[], int m, int i) {
void sift (double a[], int b[], int m, int i) {

  int k, j, l, index;
  float val;
  
  k = i; val = a[b[i]]; index = b[i];

  do {
    
    j = k; //barometer operation

    if (2*j+1 <= m ) {
      
      l = 1;
      if ( 2*j+2 <= m && a[b[2*j+2]] > a[b[2*j+1]] )  l = 2;
      if ( a[b[2*j+l]] > val ) {
	b[k] = b[2*j+l];
	k = 2*j+l;
      }
    }

  }while (j != k);

  if (k != i ) b[k] = index;

} 


      
