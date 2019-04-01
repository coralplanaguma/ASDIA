#include <stdio.h>
#include <stdlib.h>

#include "swCost.h"
#include "swCalc.h"
#include "swGetSeq.h"
#include "swOut.h"

int main(void)
{
	char *s1 ;
	while((s1 = getSeq(1)) == NULL) {
		// nothing to do
	}
	char *s2 ;
	while((s2 = getSeq(1)) == NULL) {
	}
	printf("Sequences read:\ns1\t%s\ns2\t%s\n\n", s1, s2) ;

	/* BLOSUM62 prot subst cost with affine cost for short indels */
	struct cost *cost = costProt(-10,-0.5);
	struct matrix *mat_d = swInitMat(s1,s2);
	struct matrix *mat_v = swInitMat(s1,s2);
	struct matrix *mat_h = swInitMat(s1,s2);
	swFillMat(mat_d, cost, s1, s2);
	/* for debugging you can uncomment:
	swPrintMat(mat, s1, s2);*/
	printAlisAlt(mat_d, mat_v, mat_h, s1, s2);

	swFreeMat(mat_d);
	swFreeMat(mat_v);
	swFreeMat(mat_h);
	free(cost);
	free(s1);
	free(s2);
	return(0);
}
