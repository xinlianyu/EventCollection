/*
Copyright (c) 2000                      RIPE NCC


All Rights Reserved

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that copyright notice and this permission notice appear in
supporting documentation, and that the name of the author not be
used in advertising or publicity pertaining to distribution of the
software without specific, written prior permission.

THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING
ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS; IN NO EVENT SHALL
AUTHOR BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY
DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN
AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

/*
-------------------------------------------------------------------------------
Module Header
Filename          : Percentiles.C
Author            : Rene Wilhelm
Date              : 17-OCT-2000
Revised for Linux : 15-JUN-2001
Description       : Implementation of Percentiles class
Language Version  : C++
OSs Tested        : Solaris 2.6 , Debian Linux 2.2
Todo		  : turn this into a template class, so same code can be used
for distributions of any numeric type (int, float etc.)

Percentiles objects store collections of (now floating point) numbers,
such that in the end various 'percentile' levels of the distribution
can be retrieved.

$Id: Percentiles.C,v 1.3 2001/06/18 17:13:01 wilhelm Exp $
-------------------------------------------------------------------------------
*/

#include "Percentiles.h"

static char const rcsid[] = "$Id: ";

void Percentiles::Add(float value) {

	if (entries == size) {
		//have to resize
		if (size > 4000) {
			size += 2048;
		}
		else {
			size += size; // double the size
		}
		array = (float *)realloc(array, sizeof(float)*size);
	}
	array[entries] = value;
	entries++;
};

Percentiles::Percentiles() {

	size = 256;
	entries = 0;

	// can't use "new" operator, since we want to be able to realloc 
	array = (float *)malloc(sizeof(float)*size);
};

float Percentiles::GetLevel(float level) {
	// get specified percent level

	int element = (int)(level / 100 * entries);
	if (element == entries) {
		element--;	   // array boundaries 0..entries-1
	}
	return select(element);
};

/**************************************************************/
//
// Select the kth smallest value in the elements of array.
// The input array will be rearranged to have the return value in
// arr[k], with all smaller elements moved to arr[1..k-1] (in
// arbitrary order) and all larger elements in arr[k+1..n] (also
// in arbitrary order). The lowest value is 0!
// 
// From Numerical Recipes in C, page 342
// 
// Modifications: only one arguments k 
// others are private data members of the Percentiles Class
// 
/**************************************************************/

/* note #undef's at end of file */
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

float Percentiles::select(int k) {
	int i, ir, j, l, mid;
	float a, temp;

	if (entries == 0) return(-1);

	l = 0;
	ir = entries - 1;

	for (;;) {
		if (ir <= l + 1) {
			if (ir == l + 1 && array[ir] < array[l]) {
				SWAP(array[l], array[ir])
			}
			return array[k];
		}
		else {
			mid = (l + ir) >> 1;
			SWAP(array[mid], array[l + 1])
				if (array[l] > array[ir]) {
					SWAP(array[l], array[ir])
				}
			if (array[l + 1] > array[ir]) {
				SWAP(array[l + 1], array[ir])
			}
			if (array[l] > array[l + 1]) {
				SWAP(array[l], array[l + 1])
			}
			i = l + 1;
			j = ir;
			a = array[l + 1];
			for (;;) {
				do i++; while (array[i] < a);
				do j--; while (array[j] > a);
				if (j < i) break;
				SWAP(array[i], array[j])
			}
			array[l + 1] = array[j];
			array[j] = a;
			if (j >= k) ir = j - 1;
			if (j <= k) l = i;
		}
	}
}
#undef SWAP

