#include <stdio.h>
#include <stdlib.h>

/* Demonstration du remplissage *partiel* et *à l'envers* de strings */

int main(void)
{
	/* allocation de deux strings == tableaux de char */
	char out1[10];
	char out2[10];

/* pointeurs vers les caractères courants de out1 et out2, on
   commence a la fin des chaines et on va remonter en remplissant */
	char *out1_cur = out1 + 9;
	char *out2_cur = out2 + 9;

// une string finit par \0
	*out1_cur = '\0' ;
	*out2_cur = '\0' ;

/* on remplit quelques caracs dans out1 et out2 depuis la fin,
   en remontant d'un char chaque pointeur courant avant de remplir */

	*(--out1_cur) = 'A' ;
	*(--out2_cur) = 'A' ;

	*(--out1_cur) = 'c' ;
	*(--out2_cur) = 't' ;

	*(--out1_cur) = 'A' ;
	*(--out2_cur) = 'A' ;

	*(--out1_cur) = 'T' ;
	*(--out2_cur) = 'T' ;

	*(--out1_cur) = 'G' ;
	*(--out2_cur) = 'G' ;

/* On a rempli moins de 10 caracteres dans chaque string mais ce n'est
   pas un probleme, si on a fini de remplir (debut d'un alignement optimal?)
   on peut printer ce qui a été rempli */
	printf("out1 %s\nout2 %s\n\n", out1_cur, out2_cur);
	return(0) ;
}
