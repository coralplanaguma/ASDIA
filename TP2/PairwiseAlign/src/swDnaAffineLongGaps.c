#include <stdio.h>
#include <stdlib.h>

#include "swCost.h"
#include "swCalc.h"
#include "swGetSeq.h"
#include "swOut.h"

int main(void)
{
	char *s1 ;
	while((s1 = getSeq(0)) == NULL) {
		// nothing to do
	}
	char *s2 ;
	while((s2 = getSeq(0)) == NULL) {
	}
	printf("Sequences read:\ns1\t%s\ns2\t%s\n\n", s1, s2) ;

	/* affine cost for long gaps eg spliced RNA on genome */
	struct cost *cost = costDna(-100,-0.05);
	struct matrix *mat_d = swInitMat(s1,s2);
	struct matrix *mat_v = swInitMat(s1,s2);
	struct matrix *mat_h = swInitMat(s1,s2);

	swFillMatAlt(mat_d, mat_v, mat_h, cost, s1, s2);
	/* for debugging you can uncomment:
	   swPrintMat(mat); */
	printBestAlisAlt(mat_d, mat_v, mat_h, s1, s2);

	swFreeMat(mat_d);
	swFreeMat(mat_v);
	swFreeMat(mat_h);
	free(cost);
	free(s1);
	free(s2);
	return(0);
}
