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

typedef struct VNODE_T {
  UT_hash_handle hh;
  int v; /** key: variable node */
  int clus; /** cluster reference */
} vnode_t;

typedef struct CLUSTER_T {
  int label; /** (by the position of non-zero syndrome bits).
		 Positive value: label of this cluster.  Negative
		 value: reference label.*/
  int num_poi_v;
  point_t *first_v; /** linked list for associated v-nodes */
  point_t *last_v;
  int num_poi_c;
  point_t *first_c; /** linked list for associated c-nodes */
  point_t *last_c;
} cluster_t;

/** @brief structure with full information on cluster structure for UFL decoder */
typedef struct UFL_T {
  const int nvar;
  const int nchk;
  vnode_t * nodes; /** hash storage for occupied nodes */
  /** TODO: compare performance for taking an int array of nodes */
  int num_v; /** number of used `v_nodes` */
  int num_c; /** number of used `c_nodes` */
  int num_clus; /** number of defined clusters */
  vec_t *error; /** [`nvar`] current error vector */
  vec_t *syndr; /** [`nchk`] remaining syndrome bits */
  point_t *v_nodes; /** [`nvar`] non-zero value is the number of the next
		 `v` in cluster */
  point_t *c_nodes; /** [`nchk`] non-zero value is the number of the next
		 `c` in cluster (only non-zero syndrome bits). */
  vnode_t * spare;  //  int *label;  /** [`nchk`] */
  cluster_t clus[0];/** [`nchk`] list of clusters and associated `v` and `c` lists. */
} ufl_t;

/** @brief construct an empty ufl structure */
ufl_t *ufl_init(const params_t * const p){
  ufl_t *ans = calloc(1,sizeof(ufl_t) + sizeof(cluster_t)*p->nchk);
  if(!ans) ERROR("memory allocation");
  /** assign constant members */
  int *dummy_ =(void *) &(ans->nvar); *dummy_ = p->nvar;
  dummy_ =(void *) &(ans->nchk); *dummy_ = p->nchk;
  ans->error = vec_init(ans->nvar);
  ans->syndr = vec_init(ans->nchk);
  ans->spare = malloc(ans->nvar * sizeof(vnode_t));
  ans->v_nodes = malloc(sizeof(point_t) * ans->nvar);
  ans->c_nodes = malloc(sizeof(point_t) * ans->nchk);
  if ((ans->error==NULL) || (ans->syndr == NULL) ||
      (ans->spare == NULL) ||
      (ans->v_nodes == NULL) ||
      (ans->c_nodes == NULL))
    ERROR("memory allocation");
  return ans;
}

  /** @brief destruct a ufl structure */
  ufl_t *ufl_free( ufl_t *s){

    vnode_t *nod, *tmp;
    HASH_ITER(hh, s->nodes, nod, tmp) {
      HASH_DEL(s->nodes,nod);
    }    
    free(s->error);
    free(s->syndr);
    free(s->v_nodes);
    free(s->c_nodes);
    free(s->spare);
    free(s);
    return NULL;
  }


static inline void cluster_print(const int cl, const cluster_t * const u){
  printf("# cluster %d label=%d num_v=%d num_c=%d (%s)\n",cl, u->label, u->num_poi_v, u->num_poi_c,
	 cl == u->label ? "original" : "reference");
  if(u->num_poi_v){
    printf("# nodes_v: ");
    int cnt=0;
    for(point_t * tmp = u->first_v; tmp != NULL; tmp = tmp ->next, cnt++)      
      printf("%d%s", tmp->index, tmp->next == NULL ?
	     tmp == u->last_v ? "\n" : " invalid termination!\n" : " ");
    if(cnt!=u->num_poi_v)
      ERROR("num_v=%d cnt=%d mismatch",u->num_poi_v, cnt);
  }
  if(u->num_poi_c){
    printf("# nodes_c: ");
    int cnt=0;
    for(point_t * tmp = u->first_c; tmp != NULL; tmp = tmp ->next, cnt++){
      printf("%d%s", tmp->index, tmp->next == NULL ?
	     tmp == u->last_c ? "\n" : "invalid termination!\n" : " ");
    }
    if(cnt!=u->num_poi_c)
      ERROR("num_c=%d cnt=%d mismatch",u->num_poi_c, cnt);
  }
}

static inline int cluster_verify(const int cl, const cluster_t * const clus){
  int num_err = 0; /** error counter */
  if(clus->label != cl){ /* not native */
    if((clus->num_poi_v)||(clus->first_v)||(clus->last_v)||
       (clus->num_poi_c)||(clus->first_c)||(clus->last_c))
      printf("err %d: cluster %d label %d : non-empty\n", ++num_err, cl, clus->label);
  }
  else{ /** native */
    if(clus->num_poi_v){  /** `v` nodes present */
      int cnt=0;
      for(point_t * tmp = clus->first_v; tmp != NULL; tmp = tmp ->next, cnt++)      
	if((tmp->next == NULL)&&(tmp != clus->last_v))
	  printf("err %d: invalid termination v!\n", ++num_err);	
      if(cnt!=clus->num_poi_v)
	printf("err %d: num_v=%d cnt=%d mismatch\n",++num_err, clus->num_poi_v, cnt);
    }
    else{  /** no `v` nodes */
      if((clus->first_v)||(clus->last_v))
	printf("err %d: non-null v pointers\n", ++num_err);
    }

    if(clus->num_poi_c){  /** `v` nodes present */
      int cnt=0;
      for(point_t * tmp = clus->first_c; tmp != NULL; tmp = tmp ->next, cnt++)      
	if((tmp->next == NULL)&&(tmp != clus->last_c))
	  printf("err %d: invalid termination c!\n", ++num_err);	
      if(cnt!=clus->num_poi_c)
	printf("err %d: num_c=%d cnt=%d mismatch\n",++num_err, clus->num_poi_c, cnt);
    }
    else{  /** no `v` nodes */
      if((clus->first_c)||(clus->last_c))
	printf("err %d: non-null c pointers\n", ++num_err);
    }    
  }
  if(num_err)
    ERROR("errors encountered\n");
  return 0;
}

static inline int ufl_verify(const ufl_t * const u){
  int num_err = 0;
  int num_c = 0, num_v = 0;
  for (int ic=0; ic< u->num_clus; ic++){
    cluster_verify(ic, & u->clus[ic]);
    num_c += u->clus[ic].num_poi_c;
    num_v += u->clus[ic].num_poi_v;
  }
  if ((num_c!=u->num_c) || (num_c!=u->num_c))
    ERROR("unexpected: num_c=%d cnt=%d num_v=%d cnt=%d",u->num_c, num_c, u->num_v, num_v);
  vnode_t *nod, *tmp;
  int cnt=0;
  HASH_ITER(hh, u->nodes, nod, tmp) {
    if(cnt > u->nvar){
      num_err++; break;
    }
    cnt++;    
  }
  if((cnt != u->num_v) || (num_err))
    ERROR("counted %d expected num_v=%d num_err=%d",cnt, u->num_v, num_err);
  return 0;
}

static inline void ufl_print(const ufl_t *const u){
  printf("# ufl_strucure num_v=%d num_c=%d num_clus=%d\n",u->num_v, u->num_c, u->num_clus);
  for(int i=0 ; i< u->num_clus; i++)
    cluster_print(i,& u->clus[i]);
  vnode_t *tmp, *nod;
  printf("vnodes in hash: ");
  int cnt=0;
  HASH_ITER(hh, u->nodes, nod, tmp) {
    if(cnt > u->nvar)
      break;      
    printf("(%d %d -> %d) ",nod->v, nod->clus, u->clus[nod->clus].label);
    cnt++;
  }
  if(cnt != u->num_v)
    ERROR("counted %d expected num_v=%d",cnt, u->num_v);
}

static inline int merge_clus(const int c1, const int c2, ufl_t * const u){
  assert(c1 != c2); /** sanity check */

  u->clus[c1].num_poi_v += u->clus[c2].num_poi_v;
  u->clus[c2].num_poi_v = 0;
  if( u->clus[c2].last_v ){
    if( u->clus[c1].last_v )
      u->clus[c1].last_v -> next = u->clus[c2].first_v;
    else
      u->clus[c1].first_v = u->clus[c2].first_v;
    u->clus[c1].last_v = u->clus[c2].last_v;
    u->clus[c2].last_v = u->clus[c2].first_v = NULL;
  }

  if( u->clus[c2].last_c){
    u->clus[c1].num_poi_c += u->clus[c2].num_poi_c;
    u->clus[c2].num_poi_c = 0;
    if( u->clus[c1].last_c )
      u->clus[c1].last_c -> next = u->clus[c2].first_c;
    else
      u->clus[c1].first_c = u->clus[c2].first_c;
    u->clus[c1].last_c = u->clus[c2].last_c;
  }
  u->clus[c2].last_c = u->clus[c2].first_c = NULL;
    
  u->clus[c2].label = c1;
  return 0;
}

/** @brief add variable node `v` to cluster `cl` 
 * @return 0 if nothing was done, 1 if just added, 2 if clusters merged 
 */
static inline int add_v(const int v, int cl, ufl_t * const u){
  while(u->clus[cl].label != cl)  /** dereference */
    cl = u->clus[cl].label;
  vnode_t * nod;
  HASH_FIND_INT(u->nodes, &v, nod);
  if(nod){ /** vertex already in a cluster */
    int lbl = u->clus[nod->clus].label;
    while(u->clus[lbl].label != lbl)  /** dereference */
      lbl = u->clus[lbl].label;
    if(lbl == cl){  /** same `cl`, nothing needs to be done */
      //      printf("found %d, already in %d\n",v,cl);      
      return 0;
    }
    /** otherwise merge `cl` and `lbl` clusters */
    //    printf("found %d, already in %d !=  %d\n",v,lbl,cl);      
    int c1, c2; /** merge `c2` into `c1` */
    if (cl < lbl){
      c1=cl; c2=lbl;
    }
    else {
      c1=lbl; c2=cl;
    }
    merge_clus(c1,c2,u);
#ifndef NDEBUG    
    //    printf("verify c1=%d\n",c1);
    cluster_verify(c1, & u->clus[c1]);
    //    printf("verify c2=%d\n",c2);
    cluster_verify(c2, & u->clus[c2]);
#endif     
    return 2; /** merged clusters */
  }
  /** otherwise add `v` to `cl` */
  //  printf("not found, adding %d to %d\n",v,cl);
  point_t *tmp = & (u->v_nodes[u->num_v]);
  tmp->index = v;
  tmp->next = NULL;
  if(u->clus[cl].last_v){
    if((u->clus[cl].first_v == NULL) ||
       (u->clus[cl].num_poi_v == 0)){
      printf("invalid cluster: cl=%d last_v!=NULL num_poi_v=%d\n",cl,u->clus[cl].num_poi_v);
      ufl_print(u);
      ERROR("here");
    }
    u->clus[cl].last_v -> next = tmp;
  }
  else{ /** last_v == NULL */
    if((u->clus[cl].num_poi_v != 0) || 
       (u->clus[cl].first_v != NULL)){
      printf("invalid cluster: cl=%d last_v==NULL num_poi_v=%d\n",cl,u->clus[cl].num_poi_v);
      ufl_print(u);
      ERROR("here");
    }
    u->clus[cl].first_v = tmp;
  }
  u->clus[cl].last_v = tmp;
  u->clus[cl].num_poi_v ++; 
  vnode_t *x = &(u->spare[u->num_v ++]);
  x->v = v;
  x->clus = cl;
  HASH_ADD_INT(u->nodes,v,x);

  return 1; /** added vertex */
}

/** @brief start ufl decoding */
/** assume `u` has been initialized but may need cleaning */
  int dec_ufl_start(const vec_t * const s, ufl_t * const u, const params_t * const p){
    /** clean up */
    vnode_t *nod, *tmp;
    HASH_ITER(hh, u->nodes, nod, tmp) {
      HASH_DEL(u->nodes,nod);
    }
     for(int i=0; i< u->num_clus; i++){       
       u->clus[i].first_v = u->clus[i].last_v = NULL;
       u->clus[i].num_poi_v=0; 
       u->clus[i].first_c = u->clus[i].last_c = NULL;
       u->clus[i].num_poi_c=0; 
     } 
    u->num_v = u->num_c = u->num_clus =0;

    //    csr_out(p->mH);
    /** start a cluster for each position in `s` */
    for (int ic = 0; ic < s->wei; ic++){
      u->num_clus ++ ;
      u->clus[ic].label = ic;
      u->num_c ++; /** c nodes in use total */
      point_t *tmp = & u->c_nodes[ic];
      int c = tmp->index = s->vec[ic]; /** check node index = row of `H` */
      tmp->next = NULL;
      u->clus[ic].first_c = u->clus[ic].last_c = tmp;
      u->clus[ic].num_poi_c = 1;
      u->clus[ic].first_v = u->clus[ic].last_v = NULL;
      u->clus[ic].num_poi_v = 0;
      /** attach all neighboring vertices, merging clusters if needed */
      assert(c < p->mH->rows);
      for(int iv = p->mH->p[c]; iv < p->mH->p[c+1]; iv++){
	int v = p->mH->i[iv]; /** variable node index */
	//	printf("adding v=%d to cluster %d\n",v,ic);
	add_v(v,ic,u);
      }
    }    
#ifndef NDEBUG
    ufl_verify(u);
    if(p->debug & 32){
      printf("# after dec_ufl_start()\n# s: ");      
      vec_print(s);
      ufl_print(u);
    }
#endif

    return 0;    
  }


/** functions to define:
- [x] `merge(i,j)`: merge `j` to `i` (reference `label[j]->i`, update linked lists and the counters)  
- [x] add vertex `v` to cluster `j` (and merge if needed)
- [x] init clusters (from given syndrome vector)
  - init `v_next`, `c_next` to zero
  - init `label` to zero (?)
  - init `clus` to zero.
- [ ] grow clusters
- [x] reset clusters 
- [x] verify clusters 
- [ ] look up decoding (given the set of clusters, try to look up the syndrome vectors in the hash).
- [ ] RIS decoding: 
 - form local_HT (set of rows in the cluster)
 - add local syndrome vector 
 - see if these are linearly independent 
 - if not, transpose the matrices and do the decoding.
- [ ] items TODO: 
  - make sure the `U` hash has near-ML vectors; 
  - also, use connected clusters to increase `uW`
  - read/write the constructed clusters into / from nz list
    (non-trivial if ML)
  - if the resulting clusters are sufficiently small, verify
    decodability and use `RIS` decoding 


- [ ] simplified `pre` decoding 
  - given `s` vector, try to find it as a whole.  If success, return `e` vector
  - construct clusters
  - for each cluster, try to find the syndrome -> `e`; add `v`-nodes
    to `e`-vector.  If success in every cluster, sort and return the result.
    vector.  
- [ ] implementation
  - [ ] non-zero `uW` or `maxU` triggers `pre` (also, `do_clusters()`
        and `kill_clusters()` at the end.)
  - [ ] add LLR calculation to clusters and comparison
  - [ ] adjust counters (add counters for `pre`)
  - [ ] Read (or generate) `e` or `He` by rows (each row = syndrome vector)
  - [ ] pointer vector `poivec` of size `nvec` to store rows positions to be
        processed by BP or RIS
  - [ ] each vector, try with `pre`, if success, write `e` to
        `out-matrix`; otherwise add row number to `poivec` and increment
        the counter
  - [ ] allocate small `HeT` matrix; go over remaining rows (in
        `poivec`) and assign its elements 
  - [ ] copy `output` to proper rows of `out-matrix`
  - [ ] if internally comparing `obs`, come up with the protocol to
        update success counters.

## Additional notes:

The small parameter is roughtly `x= z_r z_c p` where `z_r` and `z_c` are row
and column weights.  This parameter is expected to remain more or less
the same for phenomenological error model and for the corresponding
circuit error model.  With clusters of size up to `m`, the
residual error rate should scale as `n*p*x^m /m!`.

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
  vec_t *v0 = p->v0, *v1 = p->v1, *tmp;
  if(!v0){
    v0 = p->v0 = vec_init(p->mHt->rows);
    v1 = p->v1 = vec_init(p->mHt->rows);
  }

  int w=0;
  for(int i=0; i < wgt; i++){
    w=vec_csr_row_combine(v1,v0,p->mHt,err[i]);
    tmp=v1;v1=v0;v0=tmp; /** `swap the pointers */    
  }
  qsort(v0->vec, w, sizeof(rci_t), cmp_rci_t);
  two_vec_t *ans = malloc(sizeof(two_vec_t)+(w+wgt)*sizeof(int));
  if(!ans)
    ERROR("memory allocation");
  for(int i=0; i < w; i++)
    ans->arr[i]=v0->vec[i];
  ans->cnt = 0;
  ans->w_e = wgt;
  ans->w_s = w;
  ans->w_tot=wgt+w;
  ans->vec = &( ans->arr[w] );
  for(int j=0 ; j < wgt; j++)
    ans->vec[j] = err[j];

 
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
  //  two_vec_print(entry);
  const size_t keylen = entry->w_e * sizeof(rci_t);
  HASH_FIND(hh, p->hashU_error, entry->vec, keylen, pvec);
  if(!pvec) /** vector not found, inserting */	     
    HASH_ADD(hh, p->hashU_error, vec, keylen, entry); 
  else
    ERROR("unexpected");

  for(int w=1; w<=wmax; w++){
    //    printf("w=%d max=%d \n",w,max);
    int cnt=0; // initial invalid vector 
    for(int j=0; j<w; j++)
      err[j]=j;
    err[w-1]=w-2;/** this is a kludge to ensure all vectors show up below */
      
    while(next_vec(w,err,max)){
      cnt++;
      entry = two_vec_init(w,err,p);
      //      two_vec_print(entry);
      const size_t keylen = entry->w_e * sizeof(rci_t);
      HASH_FIND(hh, p->hashU_error, entry->vec, keylen, pvec);
      if(!pvec) /** vector not found, inserting */	     
	HASH_ADD(hh, p->hashU_error, vec, keylen, entry); /** store in the `hash` */
      else
	ERROR("unexpected");
    }
#if 0    
    int num=1, den=1; /** expected {max \choose w} */
    for(int i=0; i<w; i++){
      num *= (max-i);
      den *= (i+1);
    }
    printf("w=%d max=%d cnt=%d expected=%d=%d/%d\n\n",w,max,cnt,num/den,num,den);
#endif     
  }
  /** TODO: use hash by syndrome (`det`) with secondary hash by `obs`
   *  to enable near-ML decoding for these low-weight clusters 
   */
  HASH_SORT(p->hashU_error, by_syndrome);

  //  printf("move entries syndrome ordering:\n");
  two_vec_t *tmp=NULL, *good=NULL;
  HASH_ITER(hh, p->hashU_error, entry, pvec) {
    if(good){
      if(by_syndrome(entry,good)==0){
	int cmp = by_error(good,entry);
	switch(cmp){
	case -1: /** good is better */
	  //	  printf(" -1, deleting entry:\n");
	  //	  two_vec_print(entry);
	  HASH_DEL(p->hashU_error, entry);
	  free(entry);
	  entry=NULL;
	  break;
	case +1: /** entry is better */
	  //	  printf(" +1, deleting good:\n");
	  //	  two_vec_print(good);
	  HASH_DEL(p->hashU_error,good);
	  free(good);
	  good=entry;
	  entry=NULL;
	  break;
	case 0:
	  ERROR(" 0, this should not happen!");
	}
      }
      else{ /** `entry` with new syndrome */
	//	printf("moving 'good' to new hash\n");
	//	two_vec_print(good);
	HASH_DEL(p->hashU_error,good);
	HASH_FIND(hh, p->hashU_syndr, good->arr, good->w_s*sizeof(int), tmp);
	if(!tmp) /** vector not found, inserting */	     
	  HASH_ADD(hh, p->hashU_syndr, arr, good->w_s*sizeof(int), good);
	else
	  ERROR("unexpected");
	good=entry;
	entry=NULL;
      }
    }
    else
      good=entry;
  }

  /** deal with final `good` entry */
  if(good){
    //    printf("moving final 'good' to new hash\n");
    //    two_vec_print(good);
    HASH_DEL(p->hashU_error,good);
    HASH_FIND(hh, p->hashU_syndr, good->arr, good->w_s*sizeof(int), tmp);
    if(!tmp) /** vector not found, inserting */	     
      HASH_ADD(hh, p->hashU_syndr, arr, good->w_s*sizeof(int), good);
    else
      ERROR("unexpected");
    good=NULL;
  }
  assert(p->hashU_error == NULL);
  free(err);
}

void kill_clusters(params_t * const p){
  two_vec_t *entry, *tmp;
  printf("removing all errors by syndrome ##################\n");
  HASH_ITER(hh, p->hashU_syndr, entry, tmp) {
    //    printf("removing\n");
    //    two_vec_print(entry);
    HASH_DEL(p->hashU_syndr,entry);
    free(entry);
  }
  assert(p->hashU_syndr == NULL);
}

void do_cl_dec(params_t * const p){
  ufl_t *ufl = ufl_init(p);
  two_vec_t *entry, *tmp;
  printf("############### start do_cl_dec()\n");
  HASH_ITER(hh, p->hashU_syndr, entry, tmp) {
    printf("decoding w_s=%d\n",entry->w_s);
    //    two_vec_print(entry);
    vec_t *svec = vec_init(entry->w_s);
    for(int i=0; i< entry->w_s; i++)
      svec->vec[i] = entry->arr[i];
    svec->wei = entry->w_s;
    dec_ufl_start(svec, ufl, p);
#ifndef NDEBUG    
    ufl_verify(ufl);
    //    printf("\n\n");
#endif    
    free(svec);
  }

  ufl_free(ufl);
}
      
    
    
