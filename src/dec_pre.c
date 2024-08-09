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

/**
- *** UFL : union-find lookup decoder 
- [x] prepare_v_v_graph (sparse form of vv connectivity graph).
- [x] given the error vector, init two_vec_t structure
- [ ] check and optionally insert vector in hash (by error vector).  
- [ ] Sort by syndrome vectors and (if multiple `e` per `s`) pick the most likely `e`; 
- [ ] Check and insert vector in hash (by syndrome).
- [ ] clean up the hash 
- [ ] recursively construct connected clusters from the list of non-zero syndrome bits, label each with smallest `v`; 
- [ ] maintain (or construct) the list of non-zero syndrome bits associated with each cluster (???)
- [ ] is it a valid cluster (decoding is possible)?
- [ ] grow clusters by 1 (or by energy incremement).  

  int v_in_clus[] list of clusters by initial position `-1`: not in a
     cluster `int`: position of the cluster head (if `in_clus[x]==x`,
     this is the head).  Several levels of references are allowed.
     The actual algorithm: after increment, go over neighbors of each position in a cluster, 
     if allowed and not yet in this cluster, add (or cluster merge).
  int c_in_clus[] 

- [ ] given the position, find the corresponding cluster, optionally update pointers 
      along the path to the head.  
      Return: -1 (not in a cluster), X>=0: cluster starting at X. 
- [ ] RW decoding in a cluster

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
  int *err = malloc((wmax+2)*sizeof(int));
  two_vec_t *entry, *pvec;
  csr_out(p->mHt);
  if((0<=wmax)&&(1)){
    entry = two_vec_init(0,err,p);
    two_vec_print(entry);
    const size_t keylen = entry->w_s * sizeof(rci_t);
    HASH_FIND(hh, p->hashU, entry->arr, keylen, pvec);
    if(!pvec){ /** vector not found, inserting */	     
      qllr_t energ=0;
      for(int i=0; i < entry->w_e; i++) 
	energ += p->vLLR[entry -> vec[i]];
      entry->energ=energ;	
      HASH_ADD(hh, p->hashU, arr, keylen, entry); /** store in the `hash` */	
    }
  }
  printf("# here here\n");

  int w=1;
  if(w<=wmax){
    for(int i=0; i < p->mHt->rows; i++){
      err[0]=i;
      printf("adding i=%d: [ ",i);
      for(int i=0;i<w;i++)
	printf("%d%s",err[i],i+1<w?" ":" ]\n");
      entry = two_vec_init(w,err,p);
      two_vec_print(entry);
      const size_t keylen = entry->w_s * sizeof(rci_t);
      HASH_FIND(hh, p->hashU, entry->arr, keylen, pvec);
      if(!pvec){ /** vector not found, inserting */	     
	//	qllr_t energ=0;
	//	for(int i=0; i < entry->w_e; i++) 
	//	  energ += p->vLLR[entry -> vec[i]];
	//   entry->energ=energ;	
	HASH_ADD(hh, p->hashU, arr, keylen, entry); /** store in the `hash` */	
      }
    }
  }
  printf("# here here2 #######################\n\n");

  if(1){
    HASH_ITER(hh, p->hashU, entry, pvec) {
      printf("deleting:\n");
      two_vec_print(entry);
      HASH_DEL(p->hashU, entry);
      printf("done.\n");
    }
  }
  
  if(0){ /** exercise `next_vec()` */ 
    int err[6], max=6;
    for(int w=1; w<=6; w++){
      printf("w=%d max=%d \n",w,max);
      int cnt=1; // initial vector 
      for(int j=0; j<w; j++)
	err[j]=j;
      for(int i=0; i<w; i++)
	printf("%d%s",err[i],i+1==w?"\n":"");
      while(next_vec(w,err,max)){
	cnt++;
	printf("w=%d cnt=%d: ",w,cnt);
	for(int i=0; i<w; i++)
	  printf("%d%s",err[i],i+1==w?"\n":"");
      }
      int num=1, den=1; /** expected {max \choose w} */
      for(int i=0; i<w; i++){
	num *= (max-i);
	den *= (i+1);
      }
      printf("w=%d max=%d cnt=%d expected=%d=%d/%d\n\n",w,max,cnt,num/den,num,den);
    }
  }
}
      
    
    
/** @brief prepare syndrome vectors for small clusters */
long long int xdo_clusters(const csr_t * const mH, const csr_t * const mHT, const params_t * const p){
  
  return 0;
}
