#ifndef VECDEC_H
#define VECDEC_H
/**
 * @file utils.h
 *
 * @brief Collection of utility functions
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

/** @brief structure to read in one line of data */
typedef struct ONE_PROB_T {
  double p; /**< probability */
  int n1;   /**< number of entries in H column */
  int n2;   /**< number of entries in L column */
  int idx[]; /**< flexible array to store `n1+n2` entries */
} one_prob_t;

/** global variables */

typedef struct PARAMS_T {
  int nrows; /** rows in `H` or `r` (set in the input file) */
  int ncws;  /** how many codewords `k` (set in the input file) */
  int n;     /** columns in H `n` (set in the input file) */
  int colw;  /** max column weight (default: `10`, increase at the command line if needed) */ 
  int steps; /** how many random decoding steps, default: `1` */
  int nvec;  /** max number of syndromes to process in one bunch */
  int debug; /** `debug` information */ 
  char *fin; /**< `input file` name for error model */
  int seed;  /**< rng `seed`, set=0 for automatic */
} params_t;




// extern params_t prm;
//params_t prm={1,3,1,NULL,NULL,0,0,0,0,0,0,0,0};

/** 
 * @brief The help message. 
 * 
 * This has to be checked and updated, especially the `debug` options.
 */
#define USAGE                                                           \
  "%s: a simple vectorized random information set decoder.\n"           \
  "  usage: %s param=value [[param=value] ... ]\n"			\
  "\t Command line arguments are processed in the order given.\n"	\
  "\t Supported parameters:\n"						\
  "\t f=[string]\t: name of the input file with the error model\n"      \
  "\t steps=[integer]\t: how many random window decoding cycles (default: 1)\n" \
  "\t nvec =[integer]\t: suggested max vector size for decoding\n"      \
  "\t colw =[integer]\t: maximum column weight (default: 0)\n"          \
  "\t seed=[integer]\t: RNG seed or use time(NULL) if 0 (default)\n"	\
  "\t debug=[integer]\t: bitmap for aux information to output\n"	\
  "\t\t*   0: clear the entire debug bitmap to 0.\n"                    \
  "\t\t*   1: output misc general info (on by default)\n"		\
  "\t See program documentation for input file syntax.\n"               \
  "\t Multiple `debug` parameters are XOR combined except for 0.\n"	\
  "\t Use debug=0 as the 1st argument to suppress all debug messages.\n"


#ifdef __cplusplus
}
#endif

#endif /* VECDEC_H */
