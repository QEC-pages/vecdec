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

#ifndef linux 
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
#endif /* linux */

/** @brief just a minimum of two integers */
inline int min(const int a, const int b){ return a<b? a:b; }
/** @brief a maximum of two integers */
inline int max(const int a, const int b){ return a>b? a:b; }

#define _maybe_unused __attribute__((unused)) 

#ifdef __cplusplus
}
#endif

#endif /* UTILS_H */
