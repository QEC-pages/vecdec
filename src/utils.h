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
