/************************************************************************ 
 * @brief functions for exact star-polygon transformations 
 * 
 * code modification functions using partial summation over certain qubits, 
 * see L.P. Pryadko "On maximum-likelihood decoding..." arXiv:1909.06732
 *
 * author: Leonid Pryadko <leonid.pryadko@ucr.edu> 
 * some code borrowed from various sources 
 ************************************************************************/
#include <m4ri/m4ri.h>
#include <stdio.h>
// #include <"copy_m4ri.h"
#include "mmio.h"
#include <ctype.h>

#include "utils.h"
#include "util_m4ri.h"
#include "qllr.h"

/** @brief structure to store one column of DEM (using LLRs) */
typedef struct ONE_PROB_T {
  qllr_t llr; /**< LLR corresponding to the probability value */
  int wH;   /**< number of entries in H column */
  int wL;   /**< number of entries in L column */
  int idx[]; /**< flexible array to store `n1+n2` entries */
} one_prob_t;




/** @brief transform LLR coefficients `[A,B;C]:=0.5*(boxplus(A+B,C)+boxplus(A-B,C));`
 * @param A the input LLR 
 * @param B the input LLR 
 * @param C the input LLR
 * @return the transformed LLR
 */
qllr_t transform3(const qllr_t A, const qllr_t B, const qllr_t C){
  const qllr_t C1=boxplus(A+B,C);
  const qllr_t C2=boxplus(A-B,C);
  const double CC = 0.5 * (dbl_from_llr(C1) + dbl_from_llr(C2));
  return llr_from_dbl(CC);
}


/** @brief replace the DEM with star-triangle transformed 
 *
 * 
 * `zero-column removal`:
 * any column which is zero both in `Hx` and  `Lx` (weight-1 row in `Hz`) is just dropped
 * 
 * `inverse decoration transformation`: 
 * from a pair of identical columns in `Hx` and `Lx` (weight-2 row in `Hz`) keep only one, 
 * with `K=boxplus(K1,K2)`
 * 
 * `star-triangle`: 
 * for every triplet of columns `(i1,i2,i3)` which sum to zero both in
 * `Hx` and `Lx` [weight-3 row in `Hz`], in effect, replace these columns 
 *  of `Hz` with `(b2+b3,b1+b2,b1+b3)`.
 * 1. set to zero the column `i1` in these matrices, 
 * 2. add weight-3 row supported on `(i1,i2,i3)` to Hx,
 * 3. modify LLRs: `K1<-[K2,K3;K1]`, `K2<-[K1,K3;K2]`, `K3<-[K1,K2;K3]`, 
 *    where `[A,B;C]:=0.5*(boxplus(A+B,C)+boxplus(A-B,C));`
 * The results depend on the order of columns.  
 * Run several times to ensure no 3-cycles remain.

 * @param[in,out] Ht transposed matrix `Hx` -> the result
 * @param[in,out] Lt transposed matrix `Lx` -> the result
 * @param[in,out] initial LLR vector -> the modified LLR vector
 * @param debug bitmap for debugging info
 * @return the new number of columns (???) 
 */
int star_triangle(csr_t * Ht, csr_t * Lt, qllr_t *LLR, const one_vec_t * const codewords,
		  _maybe_unused const long int debug){
    /** sanity check */
  assert(Ht->nz == -1); 
  assert(Lt->nz == -1); /** CSR, not LOP form */
  assert(Ht->rows == Lt->rows);
  const int n=Ht->rows;
  mzd_t *used = mzd_init(1,n);
  one_prob_t **cols = calloc(n, sizeof(one_prob_t *)); /* storage for processed columns */
  int n1=0, r1=Ht->cols, nzH1=0, nzL1=0; /** new `n`, rows of new `H`, and non-zero elements */
  if((!used) || (!cols))
    ERROR("memory allocation failed!");

  for(one_vec_t const * pvec = codewords; pvec != NULL; pvec=(one_vec_t *)(pvec->hh.next)){

    if((pvec->weight ==1) && (mzd_read_bit(used, 0, pvec->arr[0])==0)){
      /** `zero-column removal`: just mark off this column */
      mzd_flip_bit(used, 0, pvec->arr[0]); 
    }
    else if((pvec->weight ==2) &&
	    (mzd_read_bit(used, 0, pvec->arr[0])==0) &&
	    (mzd_read_bit(used, 0, pvec->arr[1])==0)){
      /** `inverse decoration transformation`: combine two identical columns into one */
      for(int i=0; i<2; i++)
	mzd_flip_bit(used, 0, pvec->arr[i]);
      int idx = pvec->arr[0]; /** working on this column */      
      int wH =  Ht->p[idx+1] - Ht->p[idx] ;
      nzH1 += wH;
      int wL = Lt->p[idx+1] - Lt->p[idx];
      nzL1 += wL;
      if((cols[n1] = malloc(sizeof(one_prob_t)+(wH+wL)*sizeof(int)))==NULL)
	ERROR("memory allocation failed!");
      cols[n1]->llr = boxplus(LLR[pvec->arr[0]],LLR[pvec->arr[0]]);
      cols[n1]->wH = wH;
      cols[n1]->wL = wL;
      int nz=0;
      for(int j=Ht->p[idx]; j < Ht->p[idx+1]; j++)
	cols[n1]->idx[nz++] = Ht->i[j];
      for(int j=Lt->p[idx]; j < Lt->p[idx+1]; j++)
	cols[n1]->idx[nz++] = Lt->i[j];
      n1++;
    }
    else if((pvec->weight ==3) &&
	    (mzd_read_bit(used, 0, pvec->arr[0])==0) &&
	    (mzd_read_bit(used, 0, pvec->arr[1])==0) &&
	    (mzd_read_bit(used, 0, pvec->arr[2])==0)){
      for(int i=0; i<3; i++)
	mzd_flip_bit(used, 0, pvec->arr[i]); 

      for(int i=0; i<3; i++){
	const int i0=(i+1)%3, i1=(i+2)%3, i2=(i+0)%3; 
	int idx = pvec->arr[i2]; /** working on this column */
	qllr_t llr = transform3(LLR[pvec->arr[i0]], LLR[pvec->arr[i1]], LLR[pvec->arr[i2]]);
	int wH = i2==0 ? 1 : Ht->p[idx+1] - Ht->p[idx] + 1 ;
	nzH1 += wH;
	int wL = i2==0 ? 0 : Lt->p[idx+1] - Lt->p[idx];
	nzL1 += wL;
	if((cols[n1] = malloc(sizeof(one_prob_t)+(wH+wL)*sizeof(int)))==NULL)
	  ERROR("memory allocation failed!");
	cols[n1]->llr = llr;
	cols[n1]->wH = wH;
	cols[n1]->wL = wL;
	int nz=0;
	if(i2!=0)
	  for(int j=Ht->p[idx]; j < Ht->p[idx+1]; j++)
	    cols[n1]->idx[nz++] = Ht->i[j];
	cols[n1]->idx[nz++] = r1; /** new (1,1,1) row */
	assert(nz==wH);
	if(i2!=0)
	  for(int j=Lt->p[idx]; j < Lt->p[idx+1]; j++)
	    cols[n1]->idx[nz++] = Lt->i[j];
	assert(nz==wH+wL);
	n1++;
      }
      r1++;
    }     
  }
  /** go over the remaining columns */
  for (int idx=0; idx<n; idx++){
    if(mzd_read_bit(used, 0, idx)==0){
      mzd_flip_bit(used, 0, idx);
      int wH = Ht->p[idx+1] - Ht->p[idx];
      nzH1 += wH;
      int wL = Lt->p[idx+1] - Lt->p[idx];
      nzL1 += wL;
      if((cols[n1] = malloc(sizeof(one_prob_t)+(wH+wL)*sizeof(int)))==NULL)
	ERROR("memory allocation failed!");
      cols[n1]->llr = LLR[idx];
      cols[n1]->wH = wH;
      cols[n1]->wL = wL;
      int nz=0;
      for(int j=Ht->p[idx]; j < Ht->p[idx+1]; j++)
	cols[n1]->idx[nz++] = Ht->i[j];
      for(int j=Lt->p[idx]; j < Lt->p[idx+1]; j++)
	cols[n1]->idx[nz++] = Lt->i[j];
      n1++;
    }
  }

  assert(n1<=n);//!    `ERROR`("unexpected n1=%d > n=%d",n1,n);
  /** move the resulting data to `Ht`, `Lt`, and `LLR` */
  csr_init(Ht,n1,r1,      nzH1); /** only reallocate `p` and `i` arrays */
  csr_init(Lt,n1,Lt->cols,nzL1);
  int jH=0, jL=0;
  for(int ir = 0; ir < n1; ir++){
    int j=0;
    Ht->p[ir]=jH;
    for(/*none*/ ; j < cols[ir]->wH; j++)
      Ht->i[jH++] = cols[ir]->idx[j];
    Lt->p[ir]=jL;
    const int max = cols[ir]->wL + cols[ir]->wH;
    for(/*none*/ ; j < max; j++)
      Lt->i[jL++] = cols[ir]->idx[j];
    LLR[ir] = cols[ir]->llr;
  }
  Ht->p[n1]=jH;
  Ht->nz = -1;
  Lt->p[n1]=jL;
  Lt->nz = -1;  

  printf("################ here!\n"); fflush(stdout);
  
  /** memory clean-up */
  mzd_free(used);
  for(int i=0; i<n1; i++)
    if(cols[i])
      free(cols[i]);
  free(cols);
  
  return 0;
}

/** @brief create K=Lz matrix with minimum weight rows from a list of codewords in hash */
csr_t * do_K_from_C(const csr_t * const mLt, const one_vec_t * const codewords,
		    const int k, const int n, const int minW, int maxW,
		    _maybe_unused const int debug){
  int classical = (mLt == NULL) ? 1 : 0;
  //printf("k=%d minW=%d maxW=%d\n",k,minW,maxW);
  mzd_t *vLt = NULL;
  if (!classical)
    vLt = mzd_init(k, mLt->cols);
  else 
    vLt = mzd_init(k, n);
  one_vec_t const **list = malloc(k*sizeof(one_vec_t *)); /** chosen vectors (pointers) */
  if(!list) ERROR("memory allocation");
  one_vec_t const * pvec;
  int rank=0, nz=0; /** how many rows already found; `nz` = non-zero bits total */
  int skip_checkW = minW > maxW ? 1 : 0;
  if (maxW < minW)
    maxW=minW; 
  for(int iw = minW; iw <= maxW && rank < k; iw++){
    if(debug&2)
      printf("rank=%d of %d; checking codewords of weight %d\n", rank, k, iw);
    for(pvec = codewords; pvec != NULL; pvec=(one_vec_t *)(pvec->hh.next)){
      if(skip_checkW || (pvec->weight == iw)){
	//	printf("iw=%d rank=%d ",iw,rank); print_one_vec(pvec);
	mzd_row_clear_offset(vLt, rank, 0); 
	if(!classical){
	  for(int i=0; i < pvec->weight; i++){
	    const int pos = pvec->arr[i];
	    for(int j=mLt->p[pos]; j < mLt->p[pos+1] ; j++)
	      mzd_flip_bit(vLt,rank,mLt->i[j]);
	  }
	}
	else{ /* classical */
	  for(int i=0; i < pvec->weight; i++){
	    const int pos = pvec->arr[i];
	    mzd_flip_bit(vLt,rank,pos);
	  }    
	}
	if(!mzd_row_is_zero(vLt, rank)){	 
	  //	  mzd_print_row(vLt,rank);
	  int rk=mzd_gauss_delayed(vLt,0,1);
	  if(rk > rank){
	    list[rank++]=pvec;
	    nz += pvec->weight;
	    if(rank==k)
	      break;
	  }
	  else{
	    //	    mzd_print(vLt);
	  }
	}
      }
    }
  }
  if(rank<k){
    printf("we have minW=%d maxW=%d (check dW parameter)\n",minW,maxW);
    ERROR("Number of codewords is not sufficient to construct K=Lz matrix");
  }
  
  mzd_free(vLt);

  /** create CSR matrix from the list of vectors */
  csr_t *ans = csr_init(NULL,k,n,nz); 
  for(int j=0, pos=0; j<k; j++){ /** row index */
    for(int i=0; i < list[j]->weight; i++){
      int idx=list[j]->arr[i]; /** col index */
      ans->i[pos++] = idx;
    }
    ans->p[j+1] = pos;
  }
  ans->nz = -1;/** csr compressed form */
  free(list);
  return ans;
}

/** @brief create G=Hz matrix with minimum weight rows from a list of codewords in hash */
csr_t * do_G_from_C(const csr_t * const mLt, const one_vec_t * const codewords,
		    const int num_need, const int minW, int maxW,
		    _maybe_unused const int debug){
  assert(mLt);
  assert(codewords);
  //printf("k=%d minW=%d maxW=%d\n",k,minW,maxW);
  mzd_t *mat = mzd_init(num_need, mLt->rows); /** needed to ensure the sufficient rank */
  mzd_t *vLt = mzd_init(1, mLt->cols); 
  one_vec_t const * pvec;
  int rank=0, num=0, nz=0; /** `rank` how many rows already in mat; `num` vectors, `nz` = non-zero bits total */
  int skip_checkW = minW > maxW ? 1 : 0; /** this may happen if `dW=-1`, i.e., include any weight */
  if (maxW < minW)
    maxW=minW;
  printf("minW=%d maxW=%d \n",minW, maxW);
  /** first round: count codewords, non-zero entries, and ensure sufficient rank */
  for(int iw = minW; iw <= maxW ; iw++){
    for(pvec = codewords; pvec != NULL; pvec=(one_vec_t *)(pvec->hh.next)){
      mzd_row_clear_offset(vLt, 0, 0); 
      if(rank < num_need)
	mzd_row_clear_offset(mat, rank, 0); 
      if(skip_checkW || (pvec->weight == iw)){
	for(int i=0; i < pvec->weight; i++){
	  const int idx = pvec->arr[i];
	  for(int j=mLt->p[idx]; j < mLt->p[idx+1] ; j++)
	    mzd_flip_bit(vLt,0,mLt->i[j]);
	}
	if(mzd_row_is_zero(vLt, 0)){
	  num++;
	  nz+=pvec->weight;
	  if(rank < num_need){
	    for(int i=0; i < pvec->weight; i++)
	      mzd_flip_bit(mat, rank, pvec->arr[i]);
	    rank=mzd_gauss_delayed(mat,0,1);
	  }
	  //	  printf("num=%d nz=%d rank=%d\n",num,nz,rank);
	}	
      }
    }
  }
  if(rank<num_need)
    ERROR("Number of codewords is not sufficient to construct G=Hz matrix: rk=%d < rankG=%d",rank,num_need);
  
  mzd_free(mat);

  /** second round: actually create the CSR matrix */
  csr_t *ans = csr_init(NULL,num,mLt->rows,nz);
  int jvec=0; /** vector index */
  int pos=0; /** position index */
  for(int iw = minW; iw <= maxW ; iw++){ /** same cycle as before  */
    for(pvec = codewords; pvec != NULL; pvec=(one_vec_t *)(pvec->hh.next)){
      if(skip_checkW || (pvec->weight == iw)){
	for(int i=0; i < pvec->weight; i++){
	  const int idx = pvec->arr[i];
	  for(int j=mLt->p[idx]; j < mLt->p[idx+1] ; j++)
	    mzd_flip_bit(vLt,0,mLt->i[j]);
	}
	if(mzd_row_is_zero(vLt, 0)){ /** this vector is trivial */
	  for(int i=0; i < pvec->weight; i++){
	    const int idx=pvec->arr[i]; /** col index again */
	    ans->i[pos++] = idx;
	  }
	  ans->p[++jvec] = pos;	
	}
      }
    }	
  }
  ans->nz = -1;/** csr compressed form */

  mzd_free(vLt);
  return ans;
}

		    
/** @brief create `generator` matrix orthogonal to rows of `mH` and
 *  `mL`.
 *
 *  WARNING:  Only finds cycles of length 3 (this is OK to
 *  construct `G` matrix for a DEM constructed by `stim`.  
 *
 *  If fails, use `do_G_from_C()` instead 
 *
 * @param rankG required rank of the matrix 
 */
csr_t * do_G_matrix(const csr_t * const mHt, const csr_t * const mLt, const qllr_t LLR[], const int rankG, 
		    const int debug){
  /** sanity check */
  assert(mHt->nz == -1); 
  assert(mLt->nz == -1); /** CSR, not LOP form */
  assert(mHt->rows == mLt->rows);
  
  /** init hash for combined columns of `mH` and `mL` */
  one_vec_t *hash=NULL;
  int buflen = mLt->rows * sizeof(one_vec_t) +
    mHt->p[mHt->rows]*sizeof(int) + mLt->p[mLt->rows]*sizeof(int);
  // printf("buflen=%zu nzH=%d nzL=%d\n",buflen,mHt->p[mHt->rows], mLt->p[mLt->rows]);
  char *pos, *buf;
  pos = buf = calloc(buflen, sizeof(char));
  if(!buf)
    ERROR("memory allocation failed");
  unsigned int max=0; /** max column weight */
  one_vec_t **HL=malloc(sizeof(one_vec_t *) * mLt->rows); /**array for combined cols */
  if (!HL)
    ERROR("memory allocation failed");
  for(int i=0; i< mLt->rows; i++){
    one_vec_t *pvec = HL[i] = (one_vec_t *) pos;
    int *vec = pvec->arr;     
    unsigned int len=0, need_sorting=0;
    for(int ir=mHt->p[i]; ir< mHt->p[i+1]; ir++){ /** row `i` of `Ht` */
      vec[len]=mHt->i[ir];
      if((len>0)&&(vec[len-1]>vec[len]))
	need_sorting=1;
      len++;
    }
    for(int ir=mLt->p[i]; ir< mLt->p[i+1]; ir++){/** row `i` of `Lt` */
      vec[len]=mLt->i[ir] + mHt->cols;      
      if((len>0)&&(vec[len-1]>vec[len]))
	need_sorting=1;
      len++;
    }
    if(need_sorting)
      qsort(vec, len, sizeof(rci_t), cmp_rci_t);
    assert(sizeof(rci_t) == sizeof(int));
      
    //for(int ir=0; ir< len; ir++)  printf("%d%s",vec[ir],ir+1<len ? " " : "\n");
    pos += sizeof(one_vec_t) + len * sizeof(int);
    //    printf("i=%d len=%d pos=%ld of %zu\n",i,len,pos-buf,buflen);
    assert(pos-buf <=buflen);
    pvec->weight = len;
    pvec->energ = LLR[i];
    pvec->cnt = i; /* use it to record the column number */
    HASH_ADD(hh, hash, arr, (len*sizeof(int)), pvec);
    if(max<len)
      max=len;
  }

  /** init the pairs to construct the matrix */
  int nz=0, nzmax= 3*(mLt->rows - mLt->cols), rowsG=0;
  assert(nzmax>0);
  int_pair *prs = malloc(sizeof(int_pair)*nzmax);
  if(!prs)
    ERROR("memory allocation fail");

  /** combine columns `i` and `j` and search in hash */
  int *cvec = malloc(2*max*sizeof(int)); /* temp storage */
  if(!cvec)
    ERROR("memory allocation fail!");
  
  for(int i=0; i < mLt->rows; i++){
    const one_vec_t *cw1 = HL[i];
    const int *vec1=cw1->arr;
    assert(cw1->weight>0); /** Just in case. An all-zero columns
			       should be removed from DEM (and are
			       never there in practice for DEM files
			       produced by Stim). */
    
    for(int j=i+1; j<mLt->rows; j++){     
      const one_vec_t *cw2 = HL[j];
      const int *vec2=cw2->arr;
      
      int ir2=0;
      size_t cnum=0;
      for(int ir1=0; ir1 < cw1->weight; ir1++){
	while((ir2 < cw2->weight) && (vec2[ir2] < vec1[ir1]))
	  cvec[cnum++]=vec2[ir2++]; /** insert coords of the second vector */
	if ((ir2 < cw2->weight) && (vec2[ir2] == vec1[ir1])){        
	  ir2++; /** just skip equal bits */
	  continue;
	}
	else 
	  cvec[cnum++]=vec1[ir1]; /** insert coord of the 1st vector */
      }
      while((ir2 < cw2->weight))
	cvec[cnum++]=vec2[ir2++]; /** insert coords of the second vector */
      assert(cnum <= 2*max);
      assert(cnum>0); /** Just in case.  Would have `cnum=0` if two
			  identical columns were present in DEM file
			  ---does not happen with Stim */
      one_vec_t *cwn = NULL;
      HASH_FIND(hh, hash, cvec, (cnum*sizeof(int)), cwn);
      if((cwn)&&(j < cwn->cnt)){ /** new triplet */
	if(nz>=nzmax){
	  nzmax=2*nzmax;
	  prs=realloc(prs,nzmax*sizeof(int_pair));
	  assert(prs!=NULL);
	}
	prs[nz++]=(int_pair ){ rowsG, i };
	prs[nz++]=(int_pair ){ rowsG, j };
	prs[nz++]=(int_pair ){ rowsG, cwn->cnt };
	rowsG++;
      if(debug & 32)
	printf("found %d (%d %d %d): ", rowsG, i, j, cwn->cnt);
      //      for(int ir=0; ir < cwn->weight; ir++)
      //	printf("%d%s",cwn->arr[ir],ir+1 < cwn->weight ? " " : "\n");
      }

    }    
  }
  free(cvec);
  
  /** init `G` matrix */
  csr_t *ans=csr_from_pairs(NULL, nz, prs, rowsG, mHt->rows);
  if(debug&64)
    csr_out(ans);
  free(prs); /** not needed anymore */

  /** free the hash */  
  one_vec_t *cw, *tmp;
  HASH_ITER(hh, hash, cw, tmp) {
#if 0    
    int *vec=cw->arr;
    printf("LLR=%g w=%d ",cw->energ,cw->weight);
    for(int ir=0; ir < cw->weight; ir++)
      printf("%d%s",vec[ir],ir+1 < cw->weight ? " " : "\n");
#endif 
    HASH_DEL(hash, cw);
  }
  free(buf);
  free(HL);
  
  /** verify the rank (should be `n`-`k`) ???? *********************** */
  int got_rankG=rank_csr(ans);
  if(debug&1){
    printf("# n=%d k=%d rankG=%d\n"
	   "# Created matrix G of size %d x %d (all weight-3 rows)\n",
	   mLt->rows, mLt->cols, rankG, ans->rows, ans->cols);
  }

  if( rankG != got_rankG )    
    ERROR("Some longer cycles are missing from G\n"
	  "expect rankG=%d != got rankG=%d", rankG, got_rankG );
  /** This would require some extensive changes to the code.  DEMs
      from Stim with depolarizing noise enabled have the shortest
      cycles of length 3, and these are sufficient to generate the
      full G matrix (may not be the case with some error models). */ 

  
  /** construct the actual matrix and clean-up */
  return ans;
}

