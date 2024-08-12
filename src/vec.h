#ifndef VEC_H
#define VEC_H
/**
 * @file vec.h
 *
 * @brief vecdec - a simple vectorized decoder
 *
 * @author Leonid Pryadko (University of California, Riverside)
 *
 * Copyright (C) 2024 Leonid Pryadko
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
#include  "qllr.h" 



typedef struct VEC_T{
  int wei; /** current weight */
  int max; /** allocated */
  int vec[0];
} vec_t;

  /** @brief return empty `vec` with `max` positions reserved */
  static inline vec_t *vec_init(const int max){
    if(max<0)
      ERROR("invalid max=%d",max);
    vec_t * ans = calloc(1, sizeof(vec_t)+max*sizeof(int));
    if(!ans)
      ERROR("memory allocation");
    ans->max = max;
    //    ans->wei = 0;
    return ans;
  }
  
/** @brief print entire `vec_t` structure by pointer */
static inline void vec_print(const vec_t * const pvec){
  printf(" w=%d [ ",pvec->wei);
  for(int i=0; i < pvec->wei; i++)
    printf("%d ",pvec->vec[i]);
  printf("]\n");
}


// v1[:] = v0[:] + mat[row,:]
// @ return weight of the resulting vector   
static inline int vec_csr_row_combine(vec_t * const v1, const vec_t * const v0,
				  const csr_t * const mat, const int row){
#ifndef NDEBUG
  if ((!v1) || (!v0) || (!mat))
    ERROR("all arguments must be allocated: v1=%p v0=%p mat=%p\n",v1,v0,mat);
  if(v1 == v0)
    ERROR("the two vectors should not be the same !");
  if((row<0) || (row >= mat->rows) ||
     (v1->max < mat->cols) ||
     (v0->max < mat->cols))
    ERROR("this should not happen\n");
#endif
  int iM, i0=0, i1=0; /** iterators */
  for (iM = mat->p[row]; iM < mat->p[row+1]; iM++){
    const int ic = mat->i[iM];
    while((i0 < v0->wei) && (v0->vec[i0] < ic))
      v1->vec[i1++] = v0->vec[i0++];
    if(i0 >= v0->wei)
      break;
    if(v0->vec[i0]==ic)
      i0++; /** `1+1=0` just skip this position */
    else
      v1->vec[i1++] = ic;
  }
  if(i0 >= v0->wei) /** remaining `mat[row]` entries */
    for (                      ; iM < mat->p[row+1]; iM++){
      const int ic = mat->i[iM];
      v1->vec[i1++] = ic;
    }
  else /** remaining `v0` entries */
    while(i0 < v0->wei)
      v1->vec[i1++] = v0->vec[i0++];
  
  v1->wei = i1;
  return i1; /** weith of the out vector */
}

/** @brief insert `j` (originally absent) into ordered array, return position */
static inline int vec_ordered_ins(vec_t * const err, const int j){
  int pos=err->wei-1;
  while(j < err->vec[pos]){
    err->vec[pos+1] = err->vec[pos];
    pos--;
  }
#ifndef NDEBUG  
  if (j == err->vec[pos]) 
    ERROR("Unexpected! vec[%d]=%d is already present!",pos,j);
#endif   
  err->vec[pos+1]=j;
#ifndef NDEBUG
  if(err->wei >= err->max)
    ERROR("unexpected!"); /** before increment */
  for(int i=0; i < err->wei; i++)
    if(err->vec[i] >= err->vec[i+1]){
      printf("check ordering at i=%d! ",i); /** before increment */
      vec_print(err);
      ERROR("unexpected");
    }
#endif   
  err->wei ++;
  return pos+1;
}


/** @brief find `val` in ordered array 
 *  @return position of `val` if found, -1 otherwise 
*/
static inline int vec_ordered_search(vec_t * const err, const int val){
  /** binary search for pos of `val` */
  int bot=0, top=err->wei , mid=0;
#ifndef NDEBUG  
  if (!top)
    return -1;
#endif   
  while(top - bot > 1){
    mid = (top+bot) >> 1;
#ifndef NDEBUG  
  if (mid>=err->wei)
    ERROR("this should not happen");
#endif   
    if (err->vec[mid] <= val)
      bot = mid;
    else
      top = mid;
    //    printf("bot=%d mid=%d top=%d err->vec[mid]=%d val=%d\n",bot,mid,top,err->vec[mid],val);
  }
  if ( err->vec[bot] == val)
    return bot;
  else
    return -1;
} 

/** @brief delete `val` (if originally present) from ordered array 
 *  @return 1 of `val` was found, 0 otherwise 
*/
static inline int vec_ordered_find_del(vec_t * const err, const int val){
  /** binary search for pos of `j` */
  int bot=0, top=err->wei, mid=0;
#ifndef NDEBUG  
  if (!top)
    return 0;
#endif   
  while(top - bot > 1){
    mid = (top+bot) >> 1;
    if (err->vec[mid] <= val)
      bot = mid;
    else
      top = mid;
  }
  if ( err->vec[mid] != val)
    return 0;
    //    ERROR("unexpected! value val=%d not found",val);
  
  for(int i=mid; i < err->wei; i++)
      err->vec[i] = err->vec[i+1];
  err->wei --;
  return 1;
}


/** @brief delete `val` in known position `pos` from ordered array */
static inline void vec_ordered_pos_del(vec_t * const err, _maybe_unused const int val, const int pos){
#ifndef NDEBUG
  if ((pos<0) || (pos >= err->wei) || (err->wei == 0) || (err->vec[pos] != val))
    ERROR("this should not happen!");
#endif
  err->wei --; 
  for(int i=pos; i < err->wei; i++)
      err->vec[i] = err->vec[i+1];
}

#ifdef __cplusplus
}
#endif
    
#endif /* VEC_H */
