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

/** global variables */

typedef struct PARAMS_T {
  int nrows; /** rows in `H` */
  int ncols; /** `suggested` maximum number of columns */ 
  int n; /** columns in H */
  int ncws; /** how many codewords */
  int steps; /** how many random decoding steps */
  int debug; /** `debug` information */ 
  char *fin; /**< `input file` name for error model */
  int seed;/* rng seed, set=0 for automatic */
} params_t;




// extern params_t prm;
//params_t prm={1,3,1,NULL,NULL,0,0,0,0,0,0,0,0};

/** 
 * @brief The help message. 
 * 
 * This has to be checked and updated, especially the `debug` options.
 */
#define USAGE                                                           \
  "%s: a simple vectorized decoder.\n"                                  \
  "  usage: %s param=value [[param=value] ... ]\n"			\
  "\t Command line arguments are processed in the order given.\n"	\
  "\t Supported parameters:\n"						\
  "\t steps=[integer]\t: how many random window decoding cycles (default: 1)\n" \
  "\t nrows=[integer]\t: how many rows in the check matrix\n"           \
  "\t ncws =[integer]\t: how many codewords\n"                          \
  "\t ncols =[integer]\t: suggested max number of columns when decoding\n" \
  "\t seed=[integer]\t: RNG seed or use time(NULL) if 0 (default)\n"	\
  "\t f=[string]\t: name of an input file with the model to process\n"  \
  "\t debug=[integer]\t: bitmap for aux information to output\n"	\
  "\t\t*   0: clear the entire debug bitmap to 0.\n"                    \
  "\t\t*   1: output misc general info (on by default)\n"		\
  "\t See program documentation for input file/string syntax.\n"	\
  "\t Multiple `debug` parameters are XOR combined except for 0.\n"	\
  "\t Use debug=0 as the 1st argument to suppress all debug messages.\n"


#ifdef __cplusplus
}
#endif

#endif /* VECDEC_H */
