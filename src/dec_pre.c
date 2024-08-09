#include <inttypes.h>
#include <strings.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <m4ri/m4ri.h>
#include <m4ri/mzd.h>
#include "utils.h"
#include "util_m4ri.h"
#include "vecdec.h"
#include "qllr.h"
#include "vec.h"

/** *** UFL : UNION-FIND LOOKUP DECODER *********************/

typedef struct POINT_T {
  int index; /** index of the `v` or `c` node */
  struct POINT_T *next;  /** pointer to next node or NULL */
} point_t; 
point_t *v_nodes; /** [`nvar`] non-zero value is the number of the next
		 `v` in cluster */
point_t *c_nodes; /** [`nchk`] non-zero value is the number of the next
		 `c` in cluster (only non-zero syndrome bits). */
int *label;  /** [`nchk+1`] */

typedef struct CLUSTER_T {
  int label; /** (by the position of non-zero syndrome bits).
		 Positive value: label of this cluster.  Negative
		 value: reference label.*/
  int num_v;
  point_t *first_v; 
  point_t *last_v;
  int num_c;
  point_t *first_c;
  point_t *last_c;
} cluster_t;

cluster_t *clus; /** [`nchk+1`] list of clusters and associated `v` and `c` lists. */
int max_clus; /** maximum cluster index */
int num_clus; 

/** functions to define:
- [ ] `merge(i,j)`: merge `j` to `i` (reference `label[j]->i`, update linked lists and the counters)  
- [ ] add vertex `v` to cluster `j` (and merge if needed)
- [ ] init clusters (from given syndrome vector)
  - init `v_next`, `c_next` to zero
  - init `label` to zero (?)
  - init `clus` to zero.
- [ ] grow clusters
- [ ] reset clusters 
- [ ] look up decoding (given the set of clusters, try to look up the syndrome vectors in the hash).
- [ ] RIS decoding: 
 - form local_HT (set of rows in the cluster)
 - add local syndrome vector 
 - see if these are linearly independent 
 - if not, transpose the matrices and do the decoding.

*/


csr_t * do_vv_graph(const csr_t * const mH, const csr_t * const mHT, const params_t *const p){
  int max_buf=10, nz=0;
  //! WARNING: cannot reuse this matrix!!!!
  int *buf = malloc(max_buf*sizeof(int));
  int *idx = malloc((mHT->rows + 1)*sizeof(int)); // row indices
  if((!idx)||(!buf))
    ERROR("memory allocation");
  for(int v0 = 0; v0< mHT->rows; v0++){
    idx[v0]=nz;
    int cnt=0;
    for(int ic0 = mHT->p[v0]; ic0 < mHT->p[v0+1]; ic0++){
      const int c0=mHT->i[ic0];
      assert(c0 < mH->rows);
      for(int iv1=mH->p[c0]; iv1 < mH->p[c0+1]; iv1++){
	const int v1=mH->i[iv1];
	if(nz >= max_buf){
	  max_buf*=2;
	  buf=realloc(buf,max_buf*sizeof(buf[0]));
	  if(!buf)
	    ERROR("memory allocation");
	}	  
	buf[nz++] = v1;
	cnt++;
      }
    }    
    qsort(buf+idx[v0], cnt, sizeof(rci_t), cmp_rci_t);
    int i, j;
    for(i = j = idx[v0]; j<nz; ){
      if(buf[j] == v0){ /** skip it */
	cnt--;
	j++;
      }
      else{
	if(i<j)
	  buf[i]=buf[j];
	j++;
	while( (j < nz) && (buf[i] == buf[j]) ){
	  j++; /** skip equal values */
	  cnt--;
	}
	i++;
      }
    }
    nz = idx[v0]+cnt;
    assert(nz == i);
  }
  idx[mHT->rows]=nz;
  csr_t *ans=malloc(sizeof(csr_t));
  if(!ans)
    ERROR("memory allocation");
  ans->rows = ans->cols = mHT->rows;
  ans->nz = -1;
  ans->nzmax = max_buf;
  ans->p = idx;
  ans->i = buf;

  if((p->debug&32)&&((ans->cols < 40)||(p->debug&2048))){
    printf("v_v_graph:\n");
    csr_out(ans);
  }

  return ans;
}

two_vec_t * two_vec_init(const int wgt, const int err[], params_t * const p){
  vec_t *v0 = vec_init(p->mHt->rows), *v1 = vec_init(p->mHt->rows), *v2;
  int w=0;
  for(int i=0; i < wgt; i++){
    w=vec_csr_row_combine(v1,v0,p->mHt,err[i]);
    v2=v1;v1=v0;v0=v2; /** `swap the pointers */    
  }
  //  printf("before sort: ");  vec_print(v0);
  qsort(v0->vec, w, sizeof(rci_t), cmp_rci_t);
  //  printf("after  sort: ");  vec_print(v0);
  two_vec_t *ans = malloc(sizeof(two_vec_t)+(w+wgt)*sizeof(int));
  if(!ans)
    ERROR("memory allocation");
  int i;
  for(i=0; i < w; i++)
    ans->arr[i]=v0->vec[i];
  ans->cnt = 0;
  ans->w_e = wgt;
  ans->w_s = w;
  ans->w_tot=wgt+w;
  ans->vec = &( ans->arr[i] );
  for(int j=0 ; j < wgt; j++)
    ans->vec[j] = err[j];

  free(v0);
  free(v1);
  
  return ans;  
}

/** @brief alphabetically next vector of weight `w` 
 *  0<= err[0] < err[1] < ... < err[w] < max
 * 
 * use similar approach to non-recursively construct connected error
 * clusters starting from a given check or variable node 
 * 
 * 
 */ 
int next_vec(const int w, int arr[], const int max){
  if(w>=max)
    return 0;
  int i, top=max;
  for(i=w-1; i>=0; i--){
    top--;
    if(arr[i]<top)
      break;
  }
  if(i<0)
    return 0; /** no more vectors */
  arr[i]++;
  for(int j=i;j+1<w;j++)
    arr[j+1]=arr[j]+1;
  return 1;
}

void do_clusters(params_t * const p){
  int wmax=p->uW;
  int *err = malloc((wmax)*sizeof(int));
  two_vec_t *entry, *pvec;
  int max = p->mHt->rows; 
  //  csr_out(p->mHt);

  // w=0 vector 
  entry = two_vec_init(0,err,p);
  two_vec_print(entry);
  const size_t keylen = entry->w_e * sizeof(rci_t);
  HASH_FIND(hh, p->hashU_error, entry->vec, keylen, pvec);
  if(!pvec) /** vector not found, inserting */	     
    HASH_ADD(hh, p->hashU_error, vec, keylen, entry); /** store in the `hash` */	

  for(int w=1; w<=wmax; w++){
    printf("w=%d max=%d \n",w,max);
    int cnt=0; // initial invalid vector 
    for(int j=0; j<w; j++)
      err[j]=j;
    err[w-1]=w-2;
      
    while(next_vec(w,err,max)){
      cnt++;
      entry = two_vec_init(w,err,p);
      two_vec_print(entry);
      const size_t keylen = entry->w_e * sizeof(rci_t);
      HASH_FIND(hh, p->hashU_error, entry->vec, keylen, pvec);
      if(!pvec) /** vector not found, inserting */	     
	HASH_ADD(hh, p->hashU_error, vec, keylen, entry); /** store in the `hash` */
      else
	ERROR("unexpected");
    }
    int num=1, den=1; /** expected {max \choose w} */
    for(int i=0; i<w; i++){
      num *= (max-i);
      den *= (i+1);
    }
    printf("w=%d max=%d cnt=%d expected=%d=%d/%d\n\n",w,max,cnt,num/den,num,den);
  }

  printf("# here here2 #######################\n\n");
  HASH_SORT(p->hashU_error, by_syndrome);

  printf("move entries syndrome ordering:\n");
  two_vec_t *tmp=NULL, *good=NULL;
  HASH_ITER(hh, p->hashU_error, entry, pvec) {
    if(good){
      if(by_syndrome(entry,good)==0){
	int cmp = by_error(good,entry);
	switch(cmp){
	case -1: /** good is better */
	  printf(" -1, deleting entry:\n");
	  two_vec_print(entry);
	  HASH_DEL(p->hashU_error, entry);
	  break;
	case +1: /** entry is better */
	  printf(" +1, deleting good:\n");
	  two_vec_print(good);
	  HASH_DEL(p->hashU_error,good);
	  good=entry;
	  break;
	case 0:
	  ERROR(" 0, this should not happen!");
	}
      }
      else{
	printf("moving 'good' to new hash\n");
	two_vec_print(good);
	HASH_DEL(p->hashU_error,good);
	HASH_FIND(hh, p->hashU_syndr, good->arr, good->w_s*sizeof(int), tmp);
	if(!tmp) /** vector not found, inserting */	     
	  HASH_ADD(hh, p->hashU_syndr, arr, good->w_s*sizeof(int), good);
	else
	  ERROR("unexpected");
	good=NULL;
      }
    }
    else
      good=entry;
  }

  /** deal with final `good` entry */
  printf("moving final 'good' to new hash\n");
  two_vec_print(good);
  HASH_DEL(p->hashU_error,good);
  HASH_FIND(hh, p->hashU_syndr, good->arr, good->w_s*sizeof(int), tmp);
  if(!tmp) /** vector not found, inserting */	     
    HASH_ADD(hh, p->hashU_syndr, arr, good->w_s*sizeof(int), good);
  else
    ERROR("unexpected");
  good=NULL;

  printf("removing all items by syndrome ##################\n");
  HASH_ITER(hh, p->hashU_syndr, entry, pvec) {
    printf("removing\n");
    two_vec_print(entry);
    HASH_DEL(p->hashU_syndr,entry);    
  }


}
      
    
    
