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
#include "uthash.h" /** hashing storage macros */


#define ERROR(fmt,...)							\
  do{									\
    fprintf (stderr, "%s:%d: *** ERROR in function '%s()' ***\n", __FILE__, __LINE__, __FUNCTION__); \
    printf("[31;1m " fmt " [0m \n",##__VA_ARGS__);			\
    exit(-1);								\
  }									\
  while(0)

  
/** pseudo-random number functions **********************************************/
  extern tinymt64_t tinymt;
  
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

  /** hash storage helper functions *** use `uthash.h` *************************************/


  /**< @brief structure to hold sparse vectors in a hash */
  typedef struct ONE_VEC_T {
    UT_hash_handle hh;
    double energ; /**< sum of LLRs */
    int weight; /**< number of integers in the list */
    int cnt; /** how many times this vector was encountered */
    //  size_t len; /** `weight*sizeof(int)` (is this really needed?) */
    int arr[0]; /** array of `weight` integers, the actual key  */
  } one_vec_t;

  /** @brief structure to read in one line of data */
  typedef struct ONE_PROB_T {
    double p; /**< probability */
    int n1;   /**< number of entries in H column */
    int n2;   /**< number of entries in L column */
    int idx[]; /**< flexible array to store `n1+n2` entries */
  } one_prob_t;

  /** @bried helper structure to sort by probabilities */
  typedef struct IPPAIR_T {int index; double prob; } ippair_t; 

  /** @brief print entire `one_vec_t` structure by pointer */
  void print_one_vec(const one_vec_t * const pvec);

  /** @brief compare two `one_vec_t` structures by energy */
static inline int by_energy(void *a, void *b){
  const one_vec_t *pa = (one_vec_t *) a;
  const one_vec_t *pb = (one_vec_t *) b;
  if (pa->energ < pb->energ)
    return -1;
  else if (pa->energ > pb->energ)
    return +1;
  else /** Ea == Eb */
    return 0;
}

  /** @brief helper function to sort `int`
   *  use `qsort(array, len, sizeof(rci_t), cmp_rci_t);`
   */
  static inline int cmp_rci_t(const void *a, const void *b){
    const int va= *((int *) a);
    const int vb= *((int *) b);
    return va-vb;
  }
/** @brief open a new `NZLIST` file for writing */
  FILE * nzlist_w_new(const char fnam[], const char comment[]);

/** @brief write a `one_vec_t` to an open `NZLIST` file */
  int nzlist_w_append(FILE *f, const one_vec_t * const vec);

/** @brief prepare to read from an `NZLIST` file */
  FILE * nzlist_r_open(const char fnam[], long int *lineno);

/** @brief read one item from an `NZLIST` file.
 * The structure in `vec` will be reallocated if necessary.
 * @param f an open file to read from.
 * @param[in] vec structure to hold the vector or NULL.            
 * @param[in,out] lineno pointer to current line number in the file
 * @return the pointer to the structure containing the data or NULL. */
  one_vec_t * nzlist_r_one(FILE *f, one_vec_t * vec, const char fnam[], long int *lineno);  

  /** extra io functions ******************************************************************/
void read_dem_file(char *fnam, void * ptrs[3], int debug);
  
double * dbl_mm_read(const char * const fin, int *nrows, int *ncols, int *siz, double *  arr);
  
void dbl_mm_write( char * const fout, const char fext[],
		   const int rows, const int cols, const double buf[],
		   const char comment[]);

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
