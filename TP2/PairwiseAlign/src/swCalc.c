#include "swCalc.h"
#include "mem.h"

/*******************************************
  This file fills, frees and prints matrices
 *******************************************/



/* allocate and initialize (first row and col) a matrix for SW
    alignment of strings s1 and s2
*/
struct matrix *swInitMat(char *s1, char *s2){
  /*
  on allocate en premier temps la matrice
  */
  struct matrix *matrice = NULL;
  matrice = malloc(sizeof(struct matrix));
  matrice->w = strlen(s1) + 1;
  matrice->h = strlen(s2) + 1;
  matrice->cells = malloc(matrice->w*matrice->h*sizeof(struct cell));

  /*
    on initialize a 0 la première ligne
  */

  for(unsigned int  col = 0; col < matrice->w; col++){
    matrice->cells[col].score = 0;
    matrice->cells[col].prevs = 0;
  }

  /*
    on initialize a 0 la première colonne
  */

  for(unsigned int  row = 0; row < matrice->h; row++){
    matrice->cells[row*matrice->w].score = 0;
    matrice->cells[row*matrice->w].prevs = 0;
  }
  return matrice;
}

/* free all allocated memory in mat */
void swFreeMat(struct matrix *mat){
  free(mat->cells);
  free(mat);
}

/* print contents of matrix, for debugging */
void swPrintMat(struct matrix *mat, char *s1, char *s2){
  for(unsigned int  i = 0; i < strlen(s1); i++){
    printf("              %c", s1[i]);
  }
  printf("\n");
  for(unsigned int  h = 0; h < (mat->w); h++){
    printf("--------------");
  }
  for(unsigned int  i = 0; i < (mat->h)-1; i++){
    printf("\n");
    printf(" |");
    for(unsigned int  j = 0; j < (mat->w); j++){
      printf(" %f, %i |", (mat->cells)[i*mat->w+j].score, (mat->cells)[i*mat->w+j].prevs);
    }
    printf("\n");
    printf("%c ", s2[i]);
    for(unsigned int  h = 0; h < (mat->w); h++){
      printf(" -------------");
    }
  }
  printf("\n");
  printf(" |");
  for(unsigned int  j = 0; j < (mat->w); j++){
    printf(" %f, %i |", (mat->cells)[(mat->h-1)*mat->w+j].score, (mat->cells)[(mat->h-1)*mat->w+j].prevs);
  }
  printf("\n");
}

/* Fill the mat matrix, using Smith-Waterman with a linear indel model
   using cost->indelOpen, or Gotoh with an affine indel model using
   cost->indelOpen and cost->indelExtend.
   Preconditions:
   - mat is correctly allocated and initialized (by swInitMat)
   - cost->subst is defined for each pair of letters in s1 and s2
*/
void swFillMat(struct matrix *mat, struct cost *cost, char *s1, char *s2){
  unsigned int w = mat->w;
  double* score = malloc(3*sizeof(double));
  uint8_t prevs = 0;
  double max;
  for(unsigned int  i = 1; i < mat->h; i++){
    for(unsigned int  j = 1; j < mat->w; j++){
        max = 0;
        prevs = 0;
        // diag
        score[0] = mat->cells[(i-1)*w+ j-1].score + cost->subst(s1[j-1], s2[i-1]);
        // left
        score[1] = mat->cells[i*w + j-1].score + cost->indelOpen;
        // top
        score[2] = mat->cells[(i-1)*w+j].score + cost->indelOpen;
        }
        for(unsigned int i=0; i<3; i++){
          if (score[i] > max){
            max = score[i];
            prevs = 1<<i;
          } else {
            if (score[i] == max){
              prevs += 1<<i;
            }
          }
        }
        mat->cells[w*i+j].score = max;
        mat->cells[w*i+j].prevs = prevs;
    }
  }
}
