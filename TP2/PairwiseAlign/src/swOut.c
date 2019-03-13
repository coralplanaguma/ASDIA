#include "swOut.h"


/* Find highest scoring local alignment(s) in mat, and print to stdout
   the corresponding best alignments.
   mat must have been filled with scores and prevs.
   cost is provided so mismatches with negative scores can be lowercased.
*/
void printBestAlis(struct matrix *mat, char *s1, char *s2){
  int **indexListe = maxIndex(mat);
  int i =0;
  while(indexListe[i] != NULL){
    printf("le meilleur chemin ayant pour score %f pour la case (%d, %d) est :\n",mat->cells[indexListe[i][0]*mat->w+indexListe[i][1]].score, indexListe[i][0],indexListe[i][1]);
    printAlis(mat, s1, s2, indexListe[i]);
    i++;
  }
}



void printAlis(struct matrix *mat, char *s1, char *s2, int* index){
  char c1[mat->w+mat->h+1];
  char c2[mat->w+mat->h+1];
  char * chaine1 = c1 + mat->w+mat->h;
  char * chaine2 = c2 + mat->w+mat->h;
  *chaine1 = '\0';
  *chaine2 = '\0';
  struct cell current_cell = mat->cells[index[0]*mat->w+index[1]];
  while(current_cell.score != 0){
    if((current_cell.prevs & 1) == 1){
      if(s1[index[1] - 1] == s2[index[0] - 1]){
        *(--chaine1) = s1[index[1] - 1];
        *(--chaine2) = s2[index[0] - 1];
      }
      else{
        *(--chaine1) = tolower(s1[index[1] - 1]);
        *(--chaine2) = tolower(s2[index[0] - 1]);
      }
      index[0]--;
      index[1]--;
      current_cell = mat->cells[index[0]*mat->w+index[1]];
    }
    else if((current_cell.prevs & 2) == 2){
        *(--chaine1) = tolower(s1[index[1] - 1]);
        *(--chaine2) = '-';
      index[1]--;
      current_cell = mat->cells[index[0]*mat->w+index[1]];
    }
    else{
      *(--chaine1) = '-';
      *(--chaine2) = tolower(s2[index[0] - 1]);
      index[0]--;
      current_cell = mat->cells[index[0]*mat->w+index[1]];
    }

  }
  printf("chaine 1 : %s\nchaine 2 : %s\n\n", chaine1, chaine2);
}



int **maxIndex(struct matrix *mat){
  int max = 0;
  int maxNumber = 0;
  for(unsigned int i = 0; i < mat->h; i++){
    for(unsigned int j = 0; j < mat->w; j++){
      if(mat->cells[i*mat->w+j].score > max){
        maxNumber = 1;
        max = mat->cells[i*mat->w+j].score;
      }
      else if(mat->cells[i*mat->w+j].score == max){
        maxNumber += 1;
      }
    }
  }
  int maxIndex = 0;
  int **index = NULL;
  index = malloc(maxNumber*sizeof(int*));
  for(int i = 0; i <maxNumber; i++){
    index[i] = malloc(2*sizeof(int));
  }
  for(unsigned int i = 0; i < mat->h; i++){
    for(unsigned int j = 0; j < mat->w; j++){
      if(mat->cells[i*mat->w+j].score == max){
        index[maxIndex][0] = i;
        index[maxIndex][1] = j;
        //printf("le max num %d est: (%d,%d)\n", maxIndex, i,j);
        maxIndex++;
      }
    }
  }
  return index;
}
