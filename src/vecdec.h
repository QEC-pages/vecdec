#ifndef VECDEC_H
#define VECDEC_H
/**
 * @file vecdec.h
 *
 * @brief vecdec - a simple vectorized decoder
 *
 * @author Leonid Pryadko (University of California, Riverside)
 *
 * Copyright (C) 2022 Leonid Pryadko
 * University of California, Riverside
 * All rights reserved.
 *
 */
#ifdef __cplusplus
extern "C"{
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/**< minimum probability for LLR calculations */
#define MINPROB (1e-7)  

/** @brief structure to read in one line of data */
typedef struct ONE_PROB_T {
  double p; /**< probability */
  int n1;   /**< number of entries in H column */
  int n2;   /**< number of entries in L column */
  int idx[]; /**< flexible array to store `n1+n2` entries */
} one_prob_t;

/** @bried helper structure to sort by probabilities */
typedef struct IPPAIR_T {int index; double prob; } ippair_t; 

/** structure to hold global variables */

typedef struct PARAMS_T {
  int nrows; /** rows in `H` or `r` (set in the `input file`) */
  int ncws;  /** how many codewords `k` (set in the `input file`) */
  int n;     /** columns in H `n` (set in the `input file`) */
  int colw;  /** max column weight (default: 10, `deprecated`) */ 
  int steps; /** number of random window decoding steps, default: `1` */
  int nvec;  /** max number of syndromes to process in one bunch (default: `16`) */
  int ntot;  /** total number of syndromes to generate (default: `1`) */
  int nfail; /** when non-zero, num of fails to terminate the run (default: `0`, do not terminate) */
  int swait; /** gauss decoding steps with no vectors changed to stop (default: `0`, do not stop) */
  int lerr;  /** local search after gauss up to this weight (default: `0`) */
  int mode;  /** operation mode, see help */
  int debug; /** `debug` information */ 
  char *fdem; /** `input file` name for detector error model (`DEM`) */
  char *fdet; /** `input file` name for detector events */
  char *fobs; /** `input file` name for observables */
  int seed;  /** rng `seed`, set=0 for automatic */
  double *vP; /** probability vector (total of `n`) */
  double *vLLR; /** vector of LLRs (total of `n`) */
  int nzH, nzL; /** count of non-zero entries in `H` and `L` */
  csr_t *mH; /** sparse version of H (by rows) */
  csr_t *mHt; /** sparse version of H (by columns) */
  csr_t *mL; /** sparse version of L (by rows) */
  csr_t *mLt; /** sparse version of L (by columns) */
  int maxJ;  /** memory to initially allocate for local storage */
  double LLRmin;
  double LLRmax;
} params_t;

extern params_t prm;

/** 
 * @brief The help message. 
 * 
 * @todo: This has to be checked and updated, especially the `debug` options.
 */
#define USAGE                                                           \
  "%s:  vecdec - vectorized decoder and LER estimator\n"                \
  "  usage: %s param=value [[param=value] ... ]\n"			\
  "\t Command line arguments are processed in the order given.\n"	\
  "\t Supported parameters:\n"						\
  "\t --help\t: give this help (also '-h' or just 'help')\n"            \
  "\t fdem=[string]\t: name of the input file with detector error model\n" \
  "\t fdet=[string]\t: input file with detector events (01 format)\n"   \
  "\t fobs=[string]\t: file with observables (01 matching lines in fdet)\n" \
  "\t\t\t (space is OK in front of file name to enable shell completion)\n" \
  "\t steps=[integer]\t: num of random window decoding steps (default: 1)\n" \
  "\t lerr =[integer]\t: local search level after gauss (0, no search)\n" \
  "\t swait=[integer]\t: steps w/o new errors to stop (0, do not stop)\n" \
  "\t nvec =[integer]\t: max vector size for decoding (default: 16)\n"  \
  "\t\t\t (list size in distance or energy calculations)\n"             \
  "\t ntot =[integer]\t: total syndromes to generate (default: 1)\n"     \
  "\t nfail=[integer]\t: total fails to terminate (0, do not terminate)\n" \
  "\t seed= [integer]\t: RNG seed or use time(NULL) if 0 (default)\n"	\
  "\t mode= [integer]\t: operation mode (default: 0)\n"                 \
  "\t\t* 0: use basic decoder; read events from file if spec'd\n"       \
  "\t\t* 1: (reserved for BP)\n"                                        \
  "\t\t* 2: generate most likely fault vectors, estimate Prob(Fail)\n"  \
  "\t debug=[integer]\t: bitmap for aux information to output (default: 1)\n" \
  "\t\t*   0: clear the entire debug bitmap to 0.\n"                    \
  "\t\t*   1: output misc general info (on by default)\n"		\
  "\t\t*   2: output matrices for verification\n"                       \
  "\t See program documentation for input file syntax.\n"               \
  "\t Multiple `debug` parameters are XOR combined except for 0.\n"	\
  "\t Use debug=0 as the 1st argument to suppress all debug messages.\n"

#ifdef __cplusplus
}
#endif

#endif /* VECDEC_H */
