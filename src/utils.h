#ifndef UTILS_H
#define UTILS_H
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

#include "tinymt64.h"

extern tinymt64_t tinymt;

#define ERROR(fmt,...)							\
  do{									\
    fprintf (stderr, "%s:%d: *** ERROR in function '%s()' ***\n", __FILE__, __LINE__, __FUNCTION__); \
    printf("[31;1m " fmt " [0m \n",##__VA_ARGS__);			\
    exit(-1);								\
  }									\
  while(0)

  

/** 
 * @brief uniformly distributed pseudo random double `x`: `0 <= x < 1`.
 * 
 * Actually, just a call to `tinymt64_generate_double()`.
 * @return a uniformly distributed pseudo random double in `[0,1)`.
 */
//double drandom(void);
static inline double drandom(void){
  return tinymt64_generate_double(&tinymt);
};

/**
 * @brief generate an exponentially-distributed random number with average =1.
 *
 * Distribution probability density is exp(-x).  Used to get next time interval
 * in a Poisson process.
 */

static inline double rnd_exponential(void){
  return -log(tinymt64_generate_doubleOO(&tinymt));
  /**< no need to check for zero or one values */
};

/** @brief Hagenauer boxplus operator 
 * 
 * WARNING: this requires full `LLR` values `x` and `y`.  Resulting
 * LLR is calculated using the equivalent form (`**verified**`):
 * `arctahn(tanh(x)*tanh(y))=log[cosh[(x+y)/2]/cosh[(x-y)/2]]`, which
 * is further simplified using the fact that `cosh` is an even
 * function, by introducing `a=abs(x+y); b=abs(x-y);` and
 * `logexp(x)=log(1-exp(-x))`, to give, finally
 * `0.5*(a-b)+logexp(a)-logexp(b);`

 *  @return LLR of the probability to have only one non-zero
 */

static inline double boxplus(const double x, const double y){
  const double a = fabs(x+y);
  const double b = fabs(x-y);
  return 0.5*(a-b)+log((1+exp(-a))/(1+exp(-b)));
  //  return 0.5*(a-b); /** min-sum version */
}

  
void read_dem_file(char *fnam, void * ptrs[3], int debug);
  
double * dbl_mm_read(const char * const fin, int *nrows, int *ncols, int *siz, double *  arr);
  
void dbl_mm_write( char * const fout, const char fext[],
		   const int rows, const int cols, const double buf[],
		   const char comment[]);

rci_t read_01(mzd_t *M, FILE *fin, rci_t *lineno, const char* fnam,
	      const int debug);  
  
#ifdef __MINGW32__ /** windows compiler */
/**
 * @brief home-grown version of unix function with the same name.
 * 
 * Read up to (and including) a newline from STREAM into `*line` (and
 * null-terminate it). On input, `*line` is a pointer returned from
 * malloc (or NULL), pointing to `*n` characters of space.  It is
 * realloc'd as necessary.  
 * 
 * @param line pointer to an array of `n` characters or NULL.
 * @param *n pointer to the size of allocated buffer or NULL.
 * @param fp a FILE open for reading.
 *
 * @return number of characters read (not including the '\0'
 *         terminator), or -1 on error or EOF.
 */
int getline(char **line, size_t *n, FILE *fp);
#endif /* __MINGW32__ */

/** @brief just a minimum of two integers */
inline int min(const int a, const int b){ return a<b? a:b; }
/** @brief a maximum of two integers */
inline int max(const int a, const int b){ return a>b? a:b; }

#define _maybe_unused __attribute__((unused)) 

#ifdef __cplusplus
}
#endif

#endif /* UTILS_H */
