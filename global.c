/* An example of global sequence alignment by dynamic programming.
 *
 * The program is ANSI C and should compile on any machine that has
 * a C compiler, with a command line like:
 *     gcc -o global global.c
 * To run the program:
 *     ./global 
 *
 * To change the sequences or the scoring system, you have to edit the
 * appropriate #define's at the top of the program, and recompile.
 * This is meant to be a simple working example of the mechanics of
 * dynamic programming sequence alignment. The details of reading
 * FASTA files and such are left as an exercise to the reader...
 *
 * SRE, Wed Jun  2 08:54:40 2004
 */



/* Define the problem here: two sequences SEQX and SEQY,
 * and the simple scoring system (MATCH, MISMATCH, and 
 * the linear gap score INDEL).
 */
#define SEQX "TTCATA"
#define SEQY "TGCTCGTA"

#define MATCH      5
#define MISMATCH  -2
#define INDEL     -6



#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int
main(void)
{
  char      *x,*y;		/* the two sequences, x_1..M and y_1..N */
  int        M,N;               /* lengths of sequences x,y             */
  int        i,j;               /* residue coords in x and y            */
  int      **S;                 /* dynamic programming matrix           */
  char      *ax,*ay;		/* result: the two aligned sequences    */
  int        K;			/* length of pairwise alignment         */
  int        k;			/* index in pairwise alignment, 0..K-1  */
  int        sc;                /* tmp variable used to find best path  */
  char       tmp;		/* tmp variable used to swap chars      */


  /*****************************************************************
   * Some initial bookkeeping.
   *****************************************************************/
  /* Initialize the sequences x,y, using SEQX and SEQY.
   *
   * In C, arrays are indexed 0..N-1, but we want our sequences to be
   * indexed 1..N, and we need row/column 0 in the DP matrix for
   * initialization conditions. That's why we do some +1 fiddling
   * here, so our sequences are x_1..x_M and y_1..y_N.
   */
  M = strlen(SEQX);
  N = strlen(SEQY);
  x = malloc(sizeof(char) * (M+2)); /* +2: leading blank, and trailing \0 */
  y = malloc(sizeof(char) * (N+2));
  strcpy(x+1, SEQX);		    /* x_1..x_M now defined; x_0 undefined */
  strcpy(y+1, SEQY);		    /* y_1..y_N now defined; y_0 undefined */

  /* Allocate the dynamic programming matrix, 0..M rows by 0..N columns.
   */
  S = malloc(sizeof(int *) * (M+1));
  for (i = 0; i <= M; i++)
    S[i] = malloc(sizeof(int) * (N+1));




  /*****************************************************************
   * Recursion: the heart of the DP algorithm.
   *****************************************************************/
  /* Initialization.
   */
  S[0][0] = 0;
  for (i = 1; i <= M; i++)  S[i][0] = i * INDEL;
  for (j = 1; j <= N; j++)  S[0][j] = j * INDEL;

  /* The dynamic programming, global alignment (Needleman/Wunsch) recursion.
   */
  for (i = 1; i <= M; i++)
    for (j = 1; j <= N; j++)
      {
                                              /* case #1: i,j are aligned */
        if (x[i] == y[j]) S[i][j] = S[i-1][j-1] + MATCH;
        else              S[i][j] = S[i-1][j-1] + MISMATCH;

	sc = S[i-1][j] + INDEL;               /* case #2: i aligned to -  */
	if (sc > S[i][j]) S[i][j] = sc;

        sc = S[i][j-1] + INDEL;               /* case #3: j aligned to -  */
        if (sc > S[i][j]) S[i][j] = sc;
      }
  /* The result (optimal alignment score) is now in S[M][N]. */






  /*****************************************************************
   * Traceback
   *****************************************************************
   * 
   * Traceback is a little fiddly, maybe a little obvious, and
   * there's more than one way to do it, so it often gets glossed
   * over in intros to DP.
   * 
   * There are two general strategies that people use for 
   * tracebacks. You can keep a "shadow matrix" of traceback
   * pointers in the recursion, so that you remember which
   * case won in every cell S(i,j). Alternatively, you can
   * recapitulate the calculations to figure out which case
   * won.
   * 
   * Recapitulation is usually a little simpler, for small examples;
   * but shadow matrices are probably more common in real-world
   * implementations. Recapitulation is used here.
   * 
   * Because we trace backwards, and we don't know the length of 
   * the alignment 'til we're done, we have a bookkeeping problem
   * with creating the aligned sequences. There's more than one
   * way to solve it. A simple (but brutish) way is to allocate
   * for the maximum length alignment (M+N; happens if x,y are
   * completely unaligned to each other), make the alignment backwards,
   * then reverse the aligned strings.
   */

  /* Allocate for our aligned strings (which will be 0..K-1)
   */
  ax = malloc(sizeof(char) * (M+N+1));
  ay = malloc(sizeof(char) * (M+N+1));

  /* Start at M,N; work backwards */
  i = M;
  j = N;
  k = 0;
  while (i > 0 && j > 0)
    {
      /* Case 1: best was x_i aligned to x_j?
       */
      if (i > 0 && j > 0) { /* watch out for boundaries */
	sc = S[i-1][j-1];
	if (x[i] == y[j]) sc += MATCH; else sc += MISMATCH;
	if (sc == S[i][j]) {              /* residue i aligns to j:     */
	  ax[k] = x[i]; ay[k] = y[j]; k++; /* record the alignment       */
	  i--; j--; 		           /* move to next cell in trace */
	  continue; 
	} 
      }

      /* Case 2: best was x_i aligned to -?
       */
      if (i > 0) {  /* watch out for boundaries */
	if (S[i-1][j] + INDEL == S[i][j]) {
	  ax[k] = x[i]; ay[k] = '-'; k++; /* record the alignment       */
	  i--; 			          /* move to next cell in trace */
	  continue; 
	}    
      }

      /* Case 3: best was - aligned to y_j?
       */
      if (j > 0) { /* watch out for boundaries */
	if (S[i][j-1] + INDEL == S[i][j]) { 
	  ax[k] = '-'; ay[k] = y[j]; k++; /* record the alignment       */
	  j--; 			          /* move to next cell in trace */
	  continue;
	}
      }
    }
  /* Done with the traceback.
   * Now ax[0..k-1] and ay[0..k-1] are aligned strings, but backwards;
   * and k is our alignment length.
   * Reverse the strings and we're done. 
   */
  K = k;			/* K = remember the alignment length */
  for (k = 0; k < K/2; k++)
    {	/* a sly, efficient in-place reversal by swapping pairs of chars: */
      tmp = ax[K-k-1]; ax[K-k-1] = ax[k]; ax[k] = tmp; 
      tmp = ay[K-k-1]; ay[K-k-1] = ay[k]; ay[k] = tmp; 
    }


  /*****************************************************************
   * Output
   *****************************************************************
   */
  printf("Sequence X: %s\n", x+1);
  printf("Sequence Y: %s\n", y+1);
  printf("Scoring system: %d for match; %d for mismatch; %d for gap\n",
	 MATCH, MISMATCH, INDEL);
  printf("\n");
  
  printf("Dynamic programming matrix:\n");
  for (j = -1; j <= N; j++)
    printf("   %c ", j > 0 ? y[j] : ' ');
  printf("\n");
  for (i = 0; i <= M; i++)
    {
      printf("   %c ", i > 0 ? x[i] : ' ');
      for (j = 0; j <= N; j++)
	printf("%4d ", S[i][j]);
      printf("\n");
    }
  printf("\n");

  printf("Alignment:\n");
  printf("X: %s\n", ax);
  printf("Y: %s\n", ay);

  printf("Optimum alignment score: %d\n\n", S[M][N]);

  /* Free the memory we allocated; then exit.
   */
  free(x);  free(y);
  free(ax); free(ay);
  for (i = 0; i <= M; i++) free(S[i]);
  free(S);
  exit(0);
}

