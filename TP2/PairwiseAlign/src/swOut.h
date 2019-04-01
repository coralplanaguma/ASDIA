#ifndef _SWOUT_H_
#define _SWOUT_H_

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "swCalc.h"
#include "swCost.h"


/* Find highest scoring local alignment(s) in mat, and print to stdout
   the corresponding best alignments.
   mat must have been filled with scores and prevs.
   cost is provided so mismatches with negative scores can be lowercased.
*/
void printAlis(struct matrix *mat, char *s1, char *s2);
void maxIndex(struct matrix *mat, unsigned int* index);
void printAlisAlt(struct matrix *mat_d, struct matrix *mat_v, struct matrix *mat_h, char *s1, char *s2);
#endif
