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

/** @brief print the cluster (its number `cl` is for information purpose) */
static inline void cluster_print(const int cl, const cluster_t * const u){
  printf("# cluster %d label=%d num_v=%d num_c=%d (%s)\n",cl, u->label, u->num_poi_v, u->num_poi_c,
	 cl == u->label ? "proper" : u->label < 0? "deleted" : "reference");
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

/** @brief verify structural integrity of a cluster */
static inline int cluster_verify(const int cl, const cluster_t * const clus){
  int num_err = 0; /** error counter */
  if(clus->label != cl){ /* not proper */
    if(clus->label >= 0){ /** not erazed */
      if((clus->num_poi_v)||(clus->first_v)||(clus->last_v)||
	 (clus->num_poi_c)||(clus->first_c)||(clus->last_c))
	printf("err %d: cluster %d label %d : non-empty\n", ++num_err, cl, clus->label);
    }
    else{
      printf("erazed cluster %d with nv=%d nc=%d\n",cl, clus->num_poi_v, clus->num_poi_c);
    }
  }
  else{ /** proper */
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

/** @brief verify structural integrity of a `ufl` and its clusters */
static inline int ufl_verify(const ufl_t * const u){
  int num_err = 0;
  int num_c = 0, num_v = 0, num_prop=0;
  for (int ic=0; ic< u->num_clus; ic++){
    cluster_verify(ic, & u->clus[ic]);
    num_c += u->clus[ic].num_poi_c;
    num_v += u->clus[ic].num_poi_v;
    if(u->clus[ic].label == ic)
      num_prop++;
  }
  if ((num_c!=u->num_c) || (num_c!=u->num_c))
    ERROR("unexpected: num_c=%d cnt=%d num_v=%d cnt=%d",u->num_c, num_c, u->num_v, num_v);
  if(num_prop != u->num_prop)
    ERROR("number of proper clusters %d mismatch, expected %d",num_prop, u->num_prop);
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

/** @brief print out the `ufl` and its clusters */
static inline void ufl_print(const ufl_t *const u){
  printf("# ufl_strucure num_v=%d num_c=%d num_clus=%d\n",u->num_v, u->num_c, u->num_clus);
  for(int i=0 ; i< u->num_clus; i++){
    cluster_print(i,& u->clus[i]);
  }
  vnode_t *tmp, *nod;
  printf("# vnodes in hash: ");
  int cnt=0;
  HASH_ITER(hh, u->nodes, nod, tmp) {
    if(cnt > u->nvar)
      break;
    int cls = nod->clus;
    printf("(%d %d",nod->v, cls);
    while((cls>0) && (cls != u->clus[cls].label)){
      cls = u->clus[cls].label;
      printf(" -> %d", cls);
    }
    printf(") ");
    cnt++;
  }
  if(cnt != u->num_v)
    ERROR("counted %d expected num_v=%d",cnt, u->num_v);
  printf("\n");
}

/** @brief merge `proper` ufl clusters `c2` into `c1` */
static inline int merge_clus(const int c1, const int c2, ufl_t * const u){
  //  printf("start of merge_clus(%d, %d, ...):\n", c1, c2);
  //  ufl_print(u);
  assert((c1 != c2)&&(u->num_prop>=2)&&
	 (c1>=0)&&(c1 < u->num_clus)&&(u->clus[c1].label==c1)&&
	 (c2>=0)&&(c2 < u->num_clus)&&(u->clus[c2].label==c2)); /** sanity check */

  u->num_prop--;
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
static inline int add_v(const int v, const int cls, ufl_t * const u){
  int cl = cls;
  if(cl<0)
    ERROR("deleted cluster %d encountered\n",cl);
  while(u->clus[cl].label != cl){  /** dereference */   
    cl = u->clus[cl].label;
    if(cl<0)
      ERROR("deleted cluster %d encountered\n",cl);
  } /** WARNING: this will not work with multiple  cluster removal / growth */
    
  vnode_t * nod;
  HASH_FIND_INT(u->nodes, &v, nod);
  if(nod){ /** vertex already in a cluster */
    assert(nod->clus >= 0); /** just in case */
    int lbl = u->clus[nod->clus].label;
    while((lbl >= 0) &&(u->clus[lbl].label != lbl))  /** dereference */
      lbl = u->clus[lbl].label;
    if(lbl == cl){  /** same `cl`, nothing needs to be done */
      //      printf("found %d, already in %d\n",v,cl);      
      return 0;
    }
    else if (lbl < 0){
      ufl_print(u);
      ERROR("deleted cluster %d encountered\n",nod->clus);
    }
    /** otherwise merge `cl` and `lbl` clusters */
    //    printf("found %d, already in %d !=  %d\n",v,lbl,cl);      
    int c1=lbl, c2=cl; /** merge `c2` into `c1` */
    if (cl < lbl){ c1=cl; c2=lbl; } /** `c1` has the smaller index */
    merge_clus(c1,c2,u);
#ifndef NDEBUG    
    cluster_verify(c1, & u->clus[c1]);
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
  else{ /** `last_v == NULL`, empty cluster*/
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

/** @brief start ufl decoding 
 *
 *  given the syndrome vector `s`, construct its cluster decomposition 
*/

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
  u->num_v = u->num_c = u->num_clus = u->num_prop = 0;
  u->error->wei = u->syndr->wei = 0;
  //  csr_out(p->mH);
  /** start a cluster for each position in `s` */
  for (int ic = 0; ic < s->wei; ic++){
    u->num_clus ++ ;
    u->num_prop ++ ; 
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
      //      printf("adding v=%d to cluster %d\n",v,ic);
      add_v(v,ic,u);
    }
  }    
#ifndef NDEBUG
  ufl_verify(u);
  if(p->debug & 32){
    printf("# after dec_ufl_start()\n# s: ");      
    vec_print(s);
    ufl_print(u);
    printf("\n");
  }
#endif

  return 0;    
}

/** @brief try look-up decoding of the syndrome vectors 
 * 
 * make sure that `p->hashU_syndr` has been constructed, 
 * and the clusters in `u` have been defined.
 *
 * @return 1 if success (i.e., found syndrome vector in each cluster)
 * still need to verify the `obs`
 */
int dec_ufl_lookup(ufl_t * const u, const params_t * const p){
  int num_prop=0;
  u->syndr->wei = u->error->wei = 0;
  for (int ic=0; ic< u->num_clus; ic++){
    if(u->clus[ic].label == ic){
      int beg_c = u->syndr->wei;
      cluster_t *clus = & u->clus[ic];
      int cnt=0;
      for(point_t * tmp = clus->first_c; tmp != NULL; tmp = tmp ->next, cnt++)
	vec_push(tmp->index, u->syndr);
      /** sort check node indices for this cluster */
      int * beg = u->syndr->vec + beg_c; 
      qsort(beg, cnt, sizeof(rci_t), cmp_rci_t);
      /** actual lookup */
      two_vec_t *tmp;
      HASH_FIND(hh, p->hashU_syndr, beg, cnt*sizeof(int), tmp);
      if(tmp){ /** vector found, add to `error` */
	for(int i=0; i<tmp->w_e; i++)
	  vec_push(tmp->err[i], u->error);
	u->clus[ic].label = -1 - ic; /**< eraze this cluster */
	u->num_prop --; 
      }
      else{ /** we dealt with this cluster but failed to match */
	num_prop++;
      }
    }
    if (num_prop == u->num_prop)
      break; /** all done, we only need to check all proper clusters */
  }
  if(u->num_prop==0){ /** success: matched syndrome in every cluster */    
    qsort(u->error->vec, u->error->wei, sizeof(rci_t), cmp_rci_t);
    qsort(u->syndr->vec, u->syndr->wei, sizeof(rci_t), cmp_rci_t);    
    assert(u->syndr->wei == u->num_c);
    return 1; /** `success`; decoded */    
  }
  /** TODO: insert some sort of residual decoding here.  Cluster growth? */
  return 0; /** some clusters remain undecoded */
}


/** functions to define:
- [ ] grow clusters
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


- [ ] implementation
  - [x] non-zero `uW` or `maxU` triggers `pre` (also, `do_clusters()`
        and `kill_clusters()` at the end.)
  - [x] add LLR calculation to clusters and comparison
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

/** @brief for `err`or vector of weight `wgt`, compute syndrome vector and init `two_vec_t` 
 */
two_vec_t * two_vec_init(const vec_t * const err, params_t * const p){
  vec_t *v0 = p->v0, *v1 = p->v1;
  if(!v0){ /** first time use => init temp storage. */
    v0 = p->v0 = vec_init(p->nchk);
    v1 = p->v1 = vec_init(p->nchk);
  }
  qllr_t energ=0;
  if(p->vLLR){
    for(int i=0; i < err->wei; i++) 
      energ += p->vLLR[err->vec[i]];
  }
  else
    energ = err->wei;
  
  v0 = csr_vec_mul(v1,v0,p->mHt,err,1);
  two_vec_t *ans = malloc(sizeof(two_vec_t)+(v0->wei + err->wei)*sizeof(int));
  if(!ans)
    ERROR("memory allocation");
  for(int i=0; i < v0->wei; i++)
    ans->arr[i]=v0->vec[i];
  ans->cnt = 0;
  ans->w_e = err->wei;
  ans->w_s = v0->wei;
  ans->w_tot = err->wei + v0->wei;
  ans->err = &( ans->arr[v0->wei] );
  for(int j=0 ; j < err->wei; j++)
    ans->err[j] = err->vec[j];
  ans->energ = energ;
  
  return ans;  
}

/** @brief alphabetically next vector of weight `w` 
 *  0<= err[0] < err[1] < ... < err[w] < max
 * 
 * TODO: use similar approach to non-recursively construct connected error
 * clusters starting from a given check or variable node 
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
  vec_t *err = vec_init(wmax);
  two_vec_t *entry, *pvec;
  int max = p->mHt->rows; 

  long long int cnt_err=0, cnt_synd=0;
  for(int w=1; w<=wmax; w++){
    //        printf("w=%d max=%d \n",w,max);
    int cnt=0; 
    for(int j=0; j<w; j++) // initial invalid vector 
      err->vec[j]=j;
    err->vec[w-1]=w-2;/** this is a kludge to ensure all vectors show up below */
    err->wei=w;
    
    while(next_vec(w,err->vec,max)){
      cnt++;
      entry = two_vec_init(err,p);
      //            two_vec_print(entry);
      const size_t keylen = entry->w_e * sizeof(rci_t);
      HASH_FIND(hh, p->hashU_error, entry->err, keylen, pvec);
      if(!pvec){ /** vector not found, inserting */	     
	HASH_ADD(hh, p->hashU_error, err, keylen, entry); /** store in the `hash` */
	if((p->maxU > 0) && (cnt_err + cnt >= p->maxU)){
	  if(p->debug&1)
	    printf("# reached limit of maxU=%lld errors\n",p->maxU);	  
	  break;
	}
      }
      else
	ERROR("unexpected");
    }
#if 0 /** generate stats */
    int num=1, den=1; /** expected {max \choose w} */
    for(int i=0; i<w; i++){
      num *= (max-i);
      den *= (i+1);
    }
    if(p->debug&2)
      printf("w=%d max=%d cnt=%d expected=%d=%d/%d\n\n",w,max,cnt,num/den,num,den);
#endif
    cnt_err += cnt;
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
	//	int cmp = by_error(good,entry);
	int cmp = by_energy(good,entry);
	switch(cmp){
	case 0: case -1: /** `good` is same or better, keep it */
	  //	  printf(" -1, deleting entry:\n");
	  //	  two_vec_print(entry);
	  HASH_DEL(p->hashU_error, entry);
	  free(entry);
	  entry=NULL;
	  break;
	case +1: /** `entry` is better, replace `good` */
	  //	  printf(" +1, deleting good:\n");
	  //	  two_vec_print(good);
	  HASH_DEL(p->hashU_error,good);
	  free(good);
	  good=entry;
	  entry=NULL;
	  break;
	default:
	  ERROR("this should not happen!");
	  break;
	}
      }
      else{ /** `entry` with new syndrome */
	cnt_synd++;
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
  cnt_synd++;
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
  if(p->debug&1)
    printf("# generated %lld errors / %lld distinct syndrome vectors\n",cnt_err, cnt_synd);
#ifndef NDEBUG
  if(p->debug&64){
    int cnt=0;
    printf("############################# error vectors in hash:\n");
    HASH_ITER(hh, p->hashU_syndr, entry, pvec) {
      cnt++;
      two_vec_print(entry);
    }
    printf("############################# total of %d\n\n",cnt);
  }    
#endif   
}

void kill_clusters(params_t * const p){
  two_vec_t *entry, *tmp;
  //  printf("removing all errors by syndrome ##################\n");
  HASH_ITER(hh, p->hashU_syndr, entry, tmp) {
    //    printf("removing\n");
    //    two_vec_print(entry);
    HASH_DEL(p->hashU_syndr,entry);
    free(entry);
  }
  assert(p->hashU_syndr == NULL);
}


int dec_ufl_one(const mzd_t * const srow, params_t * const p){
  int res=0, clus=0, ans=0;
  vec_t *err = p->err, *svec = p->svec, *v1 = p->v1;
  ufl_t *ufl = p->ufl;
  if(!p->ufl)
    ufl = p->ufl = ufl_init(p);
  if(!err)
    err = p->err = vec_init(p->nvar);
  if(!svec)
    svec = p->svec = vec_init(p->nchk);  
  if(!v1){
    v1 = p->v1 = vec_init(p->nchk);
    p->v0 = vec_init(p->nchk);
  }
  
  /** convert mzd_t *srow -> vec_t *svec */
  int idx=0, w=0;
  const word * const rawrow = srow->rows[0];
  while(((idx=nextelement(rawrow,srow->width,idx))!=-1)&&(idx<srow->ncols)){
    svec->vec[w++]=idx++;
  }
  svec->wei = w;
  
#ifndef NDEBUG    
  if(p->debug&2){
    printf("# decoding w_s=%d\n",err->wei);
    if(p->debug&64)
      vec_print(svec);
  }
#endif     
  two_vec_t *tmp;
  if(svec->wei == 0){
    cnt[CONV_TRIVIAL]++;
    ans = 1; /** trivial syndrome */
    //    ufl->syndr->wei =0;
    ufl->error->wei =0;      
  } 
  else{  /** `look up the entire syndrome -- a bit faster` */
    HASH_FIND(hh, p->hashU_syndr, svec->vec, svec->wei*sizeof(int), tmp);
    if(tmp){
      ans = 2; /** low weight match */
      cnt[CONV_LOWW]++;
      for(int i=0; i<tmp->w_e; i++)
	ufl->error->vec[i] = tmp->err[i];
      ufl->error->wei = tmp->w_e;
      for(int i=0; i<tmp->w_s; i++)
	ufl->syndr->vec[i] = tmp->arr[i];
      ufl->syndr->wei = tmp->w_s;
    }    
    else{ /** try to match as several clusters */
      dec_ufl_start(svec, ufl, p);
#ifndef NDEBUG    
      ufl_verify(ufl);
#endif
      clus = ufl->num_prop;     
      res=dec_ufl_lookup(ufl, p);
      if(res){
	ans = 3; /** cluster lookup matched syndrome */
	cnt[CONV_CLUS]++;
      }
    }
  }

  if(p->debug&16){
    switch(ans){
    case 0:
      printf("# pre-decoder match failed\n"); break;
    case 1:
      printf("# trivial syndrome vector w_s=0\n"); break;
    case 2:
      printf("# low-weight error matched w_s=%d\n",svec->wei); break;
    case 3:
      printf("# %d syndrome clusters matched\n",clus); break;
    default:
      ERROR("this should not happen!");
    }
    if(p->debug&64)
      ufl_print(ufl);
  }
  return ans;
}
		
void dec_ufl_exercise(params_t * const p){
  ufl_t *ufl = ufl_init(p);
  vec_t *err;
  if(p->useP == 0)
    ERROR("please set useP=... value\n");
  
  err=vec_init(p->nvar);
  int num_conv=0, num_conv1=0, num_corr1=0, num_corr=0, num_fail=0; /** count of converged, correct */
  for(int i=0; i< p->ntot;i++){
    vec_rnd(p->useP, err);
    

  }
  printf("# ufl lookup decoding: tot=%lld conv1=%d corr1=%d conv=%d corr=%d fail=%d p_f=%g\n",
	 p->ntot, num_conv1, num_corr1, num_conv, num_corr, num_fail, (double) num_fail/p->ntot);
  
  free(err);
  ufl_free(ufl);
}
      
    
    
