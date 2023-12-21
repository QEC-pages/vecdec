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

/** @brief replace the DEM with star-triangle transformed 
 *
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
 * @return the number of transformation steps done 
 */
int star_triangle(csr_t * Ht, csr_t * Lt, qllr_t *LLR, const long int debug){
    /** sanity check */
  assert(Ht->nz == -1); 
  assert(Lt->nz == -1); /** CSR, not LOP form */
  assert(Ht->rows == Lt->rows);

}



/** @brief create `generator` matrix orthogonal to rows of `mH` and
 *  `mL`.
 *
 *  WARNING: Currently, only finds cycles of length 3 (this is OK to
 *  construct `G` matrix for a DEM constructed by `stim`. 
 */
csr_t * do_G_matrix(const csr_t * const mHt, const csr_t * const mLt, const qllr_t LLR[], 
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
  int rankG=rank_csr(ans);
  int rankL=rank_csr(mLt); /** `k` of the code */
  int rankH=rank_csr(mHt); 
  if(debug&1){
    printf("# n=%d k=%d rankG=%d rankH=%d\n"
	   "# Created matrix G of size %d x %d (all weight-3 rows)\n",
	   mLt->rows, rankL, rankG, rankH, ans->rows, ans->cols);
  }

  if(rankH+rankL + rankG != mLt->rows )    
    ERROR("FIXME: some longer cycles are missing from G\n"
	  "n=%d != (k=%d) + (rankG=%d) + (rankH=%d)",
	  mLt->rows, rankL, rankG, rankH );
  /** This would require some extensive changes to the code.  DEMs
      from Stim with depolarizing noise enabled have the shortest
      cycles of length 3, and these are sufficient to generate the
      full G matrix (may not be the case with some error models). */ 

  
  /** construct the actual matrix and clean-up */
  return ans;
}

