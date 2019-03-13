#ifndef _SWOUT_H_
#define _SWOUT_H_

#include <stdlib.h>
#include <ctype.h>
#include "swCalc.h"
#include "swCost.h"

/* Find highest scoring local alignment(s) in mat, and print to stdout
   the corresponding best alignments.
   mat must have been filled with scores and prevs.
   cost is provided so mismatches with negative scores can be lowercased.
*/
void printBestAlis(struct matrix *mat, char *s1, char *s2) ;
void printAlis(struct matrix *mat, char *s1, char *s2, int* index);
int **maxIndex(struct matrix *mat);
#endif
