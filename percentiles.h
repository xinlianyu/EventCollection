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
Filename          : Percentiles.h
Author            : Rene Wilhelm
Date              : 17-OCT-2000
Description       : Percentiles class definition
Language Version  : C++
OSs Tested        : Solaris 2.6

Percentiles objects store collections of (now floating point) numbers,
such that in the end various 'percentile' levels of the distribution
can be retrieved (via the GetLevel() method).

$Id: Percentiles.h,v 1.1 2002/08/09 23:59:37 wilhelm Exp $
-------------------------------------------------------------------------------
*/

#ifndef  PERCENTILES_H
#define  PERCENTILES_H

#include <stdlib.h>

class Percentiles {
private:
	float *array;
	int   size;
	int   entries;
	float select(int k);
public:
	void Add(float value);	       // Add a number to the collection
	float GetLevel(float level);   // return item at specified level (0-100)
	Percentiles();
	void Reset()   { entries = 0; };   // Reset; discard old information    
	~Percentiles() { free(array); };
};

#endif

