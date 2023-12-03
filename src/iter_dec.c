/**
 *  @file iter_dec.c
 *
 * @brief iter_dec - iterative decoder for quantum and classical codes
 *
 * @author Leonid Pryadko (University of California, Riverside)
 *
 * Copyright (C) 2023 Leonid Pryadko
 * University of California, Riverside
 * All rights reserved.
 *
 *
 */

#include <inttypes.h>
#include <strings.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>
#include <time.h>
#include <unistd.h>
#include <m4ri/m4ri.h>
#include "utils.h"
#include "util_m4ri.h"
#include "iter_dec.h"

params_t prm={ .nchk=0, .nvar=0, .ncws=0, .steps=50,
  .lerr=0, .useP=0.0, //.swait=0,
  .nvec=1024,
  .ntot=1, .nfail=0, .seed=0, 
  .debug=1, .fdem=NULL, .fdet=NULL, .fobs=NULL, .finH=NULL, .finP=NULL,
  .finG=NULL, .finL=NULL, .internal=0, 
  .mode=0, 
  .LLRmin=0, .LLRmax=0, 
  .vP=NULL, .vLLR=NULL, .mH=NULL, .mHt=NULL,
  .mL=NULL, .mLt=NULL,
  .file_det=NULL, .file_obs=NULL, .line_det=0, .line_obs=0 };

typedef enum EXTR_T { TOTAL, CONV_TRIVIAL, CONV_BP, CONV_BP_AVG,
  SUCC_TRIVIAL, SUCC_BP, SUCC_OSD0, SUCC_OSD1, SUCC_OSD2,
  SUCC_TOT, EXTR_MAX } extr_t;
long long int cnt[EXTR_MAX];
long long int iter1[EXTR_MAX]; /** sums of BP iteration numbers */
long long int iter2[EXTR_MAX]; /** sums of BP iteration numbers squared */


void cnt_out(int print_banner){
  if(print_banner)
    printf("# TOTAL CONV_TRIVIAL CONV_BP CONV_BP_AVG SUCC_TRIVIAL SUCC_BP SUCC_TOT\n");
  printf(" %lld %lld %lld %lld %lld %lld %lld\n",
	 cnt[TOTAL], cnt[CONV_TRIVIAL],cnt[CONV_BP], cnt[CONV_BP_AVG],
	 cnt[SUCC_TRIVIAL], cnt[SUCC_BP], cnt[SUCC_TOT]);
}

void cnt_update(extr_t which, int iteration){
  const long long int max = (LLONG_MAX >> 3);
  cnt[which]++;
  iter1[which]+=iteration;
  iter2[which]+=iteration*iteration;
  if((cnt[which] > max) ||(iter1[which] > max) || (iter2[which]>max)){
    cnt_out(1);
    ERROR("too many iterations\n"); 
  }
}

int var_init(int argc, char **argv, params_t *p){

  int dbg=0;
  double val;
  if(argc<=1)
    ERROR("try \"%s -h\" for help",argv[0]);

  for(int i=1; i<argc; i++){  /** `debug` */
    if(sscanf(argv[i],"debug=%d",& dbg)==1){
      if(dbg==0)
	p->debug = 0;
      else{
        if(i==1)
          p->debug = dbg; /** just assign if in the `1st position` */
        else
          p->debug ^= dbg; /** otherwise `XOR` */
        if(p->debug &1)
	  printf("# read %s, debug=%d octal=%o\n",argv[i],p->debug,p->debug);
      }
    }
    else if(sscanf(argv[i],"mode=%d",& dbg)==1){
      if(dbg==0)
	p->mode = 0;
      else{
	p->mode ^= dbg;
        if(p->debug&1)
          printf("# read %s, mode=%d octal=%o\n",argv[i],p->mode,p->mode);
      }
    }
    else if (sscanf(argv[i],"nvec=%d",&dbg)==1){ /** `nvec` */
      p -> nvec = dbg;
      if (p->debug&1)
	printf("# read %s, nvec=%d\n",argv[i],p-> nvec);
    }
    else if (sscanf(argv[i],"steps=%d",&dbg)==1){ /** `steps` */
      p -> steps = dbg;
      if (p->debug&1)
	printf("# read %s, steps=%d\n",argv[i],p-> steps);
    }
    else if (sscanf(argv[i],"lerr=%d",&dbg)==1){ /** `lerr` */
      p -> lerr = dbg;
      if (p->debug&1)
	printf("# read %s, lerr=%d\n",argv[i],p-> lerr);
    }
    else if (sscanf(argv[i],"ntot=%d",&dbg)==1){ /** `ntot` */
      p -> ntot = dbg;
      if (p->debug&1)
	printf("# read %s, ntot=%d\n",argv[i],p-> ntot);
    }
    else if (sscanf(argv[i],"nfail=%d",&dbg)==1){ /** `nfail` */
      p -> nfail = dbg;
      if (p->debug&1)
	printf("# read %s, nfail=%d\n",argv[i],p-> nfail);
    }
    else if (sscanf(argv[i],"seed=%d",&dbg)==1){ /** `seed` */
      p->seed=dbg;
      if (p->debug&1)
	printf("# read %s, seed=%d\n",argv[i],p->seed);
    }
    else if (sscanf(argv[i],"useP=%lg",&val)==1){ /** `useP` */
      p -> useP = val;
      if (p->debug&1)
	printf("# read %s, useP=%g\n",argv[i],p-> useP);
    }
    else if (0==strncmp(argv[i],"f=",2)){/** back compatibility */
      if(strlen(argv[i])>2)
        p->fdem = argv[i]+2;
      else
        p->fdem = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, (fdem) f=%s\n",argv[i],p->fdem);
    }
    else if (0==strncmp(argv[i],"fdem=",5)){
      if(strlen(argv[i])>5)
        p->fdem = argv[i]+5;
      else
        p->fdem = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, fdem=%s\n",argv[i],p->fdem);
    }
    else if (0==strncmp(argv[i],"finH=",5)){
      if(strlen(argv[i])>5)
        p->finH = argv[i]+5;
      else
        p->finH = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, finH=%s\n",argv[i],p->finH);
    }
    else if (0==strncmp(argv[i],"finP=",5)){
      if(strlen(argv[i])>5)
        p->finP = argv[i]+5;
      else
        p->finP = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, finP=%s\n",argv[i],p->finP);
    }
    else if (0==strncmp(argv[i],"finL=",5)){
      if(strlen(argv[i])>5)
        p->finL = argv[i]+5;
      else
        p->finL = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, finL=%s\n",argv[i],p->finL);
    }
    else if (0==strncmp(argv[i],"finG=",5)){
      if(strlen(argv[i])>5)
        p->finG = argv[i]+5;
      else
        p->finG = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, finG=%s\n",argv[i],p->finG);
    }
    else if (0==strncmp(argv[i],"fdet=",5)){
      if(strlen(argv[i])>5)
        p->fdet = argv[i]+5;
      else
        p->fdet = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, fdet=%s\n",argv[i],p->fdet);
    }
    else if (0==strncmp(argv[i],"fobs=",5)){
      if(strlen(argv[i])>5)
        p->fobs = argv[i]+5;
      else
        p->fobs = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, fobs=%s\n",argv[i],p->fobs);
    }
    else if((strcmp(argv[i],"--help")==0)
            ||(strcmp(argv[i],"-h")==0)
            ||(strcmp(argv[i],"help")==0)){
      printf( USAGE , argv[0],argv[0]);
      exit (-1);
    }
    else if((strcmp(argv[i],"--morehelp")==0)
            ||(strcmp(argv[i],"-mh")==0)
            ||(strcmp(argv[i],"morehelp")==0)){
      printf( USAGE "\n", argv[0],argv[0]);
      printf( MORE_HELP,argv[0]);
      exit (-1);
    }
    else{ /* unrecognized option */
      printf("# unrecognized parameter \"%s\" at position %d\n",argv[i],i);
      ERROR("try \"%s -h\" for help",argv[0]);
    }

  }
 
  if (p->seed == 0){
    p->seed=time(NULL)+1000000ul*getpid(); /* ensure a different seed */
    if((p->debug)&&(p->mode!=3))
      printf("# initializing seed=%d from time(NULL)+1000000ul*getpid()\n",p->seed);
    /** use `tinymt64_generate_double(&pp.tinymt)` for double [0,1] */
  }
  // srand(time(p->seed));
  tinymt64_init(&tinymt,p->seed);

  if((! p->fdem) && (! p->finH))
    ERROR("Please specify the DEM file 'fdem' or check matrix file 'finH'\n");
  if((p->fdem) && (p->finH))
    ERROR("Please specify fdem=%s OR finH=%s but not both\n",p->fdem, p->finH);

  if(p->fdem){
    if((p->finL)||(p->finP)||(p->finG)||(p->finH))
      ERROR("DEM file fdem=%s can not be specified together with \n"
	    " finH=%s or finL=%s or finP=%s\n", p->fdem,  
	    p->finH ? p->finH : "",  p->finL ? p->finL : "",  p->finP ? p->finP : "");

    /** read in the DEM file, initialize sparse matrices */
    void * ptrs[]= {p->mH,p->mL,p->vP};
    read_dem_file(p->fdem, ptrs, p->debug);
    p->mH=ptrs[0];
    p->mL=ptrs[1];
    p->vP=ptrs[2];
    //    csr_out(p->mH);
    p->nvar = p->mH->cols;
    p->nchk = p->mH->rows;
    p->ncws = p->mL->rows;
  }
  
  if(p->finH){
    if(! p->finL){
      if((p->fobs) || (p->fdet))
	ERROR("Without L matrix, cannot specify fdet=%s or fobs=%s\n",
	      p->fdet ? p->fdet : "",  p->fobs ? p->fobs : "");
      else
	p->internal=1;
    }
    if((! p->finL) && (! p->finG)){ /** only `H` matrix specified */
      p->classical=1;
      p->internal=1;
    }
    p->mH=csr_mm_read(p->finH, p->mH, 0, p->debug);
    p->nvar = p->mH->cols;
    p->nchk = p->mH->rows;
  }
  /** at this point we should have `H` matrix for sure */
  assert(p->mH);
  
  if(p->finL){
    p->mL=csr_mm_read(p->finL, p->mL, 0, p->debug);
    p->ncws = p->mL->rows;
    if(p->mL->cols != p->nvar)
      ERROR("column number mismatch L[%d,%d] and H[%d,%d]\n",
	    p->mL->rows, p->mL->cols, p->mH->rows, p->mH->cols);
  }

  if(p->finG){
    p->mG=csr_mm_read(p->finG, p->mG, 0, p->debug);
    if(p->mG->cols != p->nvar)
      ERROR("column number mismatch G[%d,%d] and H[%d,%d]\n",
	    p->mG->rows, p->mG->cols, p->mH->rows, p->mH->cols);

    /** verify row orthogonality of `G` and `H` */
    csr_t *Gt = csr_transpose(NULL, p->mG);
    mzd_t *MGt = mzd_from_csr(NULL,Gt); /* convert to dense matrix */
    csr_free(Gt);
    if((p->mH)&&(product_weight_csr_mzd(p->mH,MGt,0)))
      ERROR("rows of H=Hx and G=Hz should be orthogonal \n");

    if(!p->mL) /** create `Lx` */
      p->mL=Lx_for_CSS_code(p->mH,p->mG);
    p->ncws = p->mL->rows;

    /** verify row orthogonality of `G` and `L` */
    if((p->mL)&&(product_weight_csr_mzd(p->mL,MGt,0)))
      ERROR("rows of L=Lx and G=Hz should be orthogonal \n");
    mzd_free(MGt);        
  }

  if(!p->mL){
    p->mL = csr_identity(p->nvar, p->nvar);
    p->ncws = p->mL->rows;
    p->classical = 1;
  }

  if(p->finP){
    int nrows, ncols, siz=0;
    p->vP=dbl_mm_read(p->finP, &nrows, &ncols, &siz, NULL);
    if(nrows*ncols != p->nvar)
      ERROR("invalid dimension=%d of P vector for H[%d,%d]\n",
	    nrows*ncols,p->mH->rows, p->mH->cols);
  }

  if(p->useP > 0){/** override */
    if(!p->vP)
      p->vP=malloc(p->nvar * sizeof(double));
    if(!p->vP)
      ERROR("memory allocation failed, n=%d",p->nvar);
    double *ptr=p->vP;
    for(int i=0;i<p->nvar;i++, ptr++)
      *ptr = p->useP;    
  }
  if(!p->vP)
    ERROR("probabilities missing, specify 'fdem', 'finP', or 'useP'");
  
  // rci_t linedet=0, lineobs=0;
  if(((p->fdet)&&(!p->fobs))||((!p->fdet)&&(p->fobs)))
    ERROR("need both fdet=%s and fobs=%s or none\n",
	  p->fdet ? p->fdet : "",  p->fobs ? p->fobs : "");
  
  if(p->fdet){/** expect both `fdet` and `fobs` to be defined */
    p->file_det=fopen(p->fdet, "r");
    if(p->file_det==NULL)
      ERROR("can't open the (det) file %s for reading\n",p->fdet);
    p->file_obs=fopen(p->fobs, "r");
    if(p->file_obs==NULL)
      ERROR("can't open the (obs) file %s for reading\n",p->fobs);
  }
  else{
    p->internal=1;
  }

  if((p->debug & 64)&&(p->nvar < 128)){ /** print matrices */   
    assert(p->nvar != 0);
    if(p->mH){
      mzd_t *mH0 = mzd_from_csr(NULL,p->mH);
      printf("matrix mH0:\n");  mzd_print(mH0);
      mzd_free(mH0);
    }
    if(p->mG){
      mzd_t *mG0 = mzd_from_csr(NULL,p->mG);
      printf("matrix mG0:\n");  mzd_print(mG0);
      mzd_free(mG0);
    }
    if(p->mL){
      mzd_t *mL0 = mzd_from_csr(NULL,p->mL);
      printf("matrix mL0:\n");  mzd_print(mL0);
      mzd_free(mL0);
    }   
    if(p->vP){      
      for(int i=0; i < p->nvar; i++)
	printf(" P[%d]=%g \n",i,p->vP[i]);
    }
  }

  /** prepare transposed matrices */
  const int n = p->nvar;
  p->mHt = csr_transpose(p->mHt, p->mH);
  p->mLt = csr_transpose(p->mLt,p->mL);
  p->vLLR = malloc(n*sizeof(qllr_t));
  assert(p->vLLR !=0);
  p->LLRmin=1e9;
  p->LLRmax=-1e9;
  for(int i=0;  i < n; i++){
    qllr_t val=llr_from_P(p->vP[i]);
    p->vLLR[i] = val;
    if(val < p->LLRmin)
      p->LLRmin = val;
    if(val > p->LLRmax)
      p->LLRmax = val;
  }
  if(p->LLRmin<=0)
    ERROR("LLR values should be positive!  LLRmin=%g LLRmax=%g",
	  dbl_from_llr(p->LLRmin),dbl_from_llr(p->LLRmax));
  if(p->debug & 2){/** `print` out the entire error model ******************** */
    printf("# error model read: r=%d k=%d n=%d LLR min=%g max=%g\n",
           p->nchk, p->ncws, p->nvar, dbl_from_llr(p->LLRmin),dbl_from_llr(p->LLRmax));
  }

  /** initialize the counters */
  for(extr_t i=0; i< EXTR_MAX; i++){
    cnt[i]=0;
    iter1[i]=0;
    iter2[i]=0;
  }

  return 0;
}

void var_kill(params_t *p){
  if(p->file_det)
    fclose(p->file_det);
  if(p->file_obs)
    fclose(p->file_obs);
  
  if(p->vP)
    free(p->vP);
  if(p->vLLR)
    free(p->vLLR);
  p->vP = p->vLLR = NULL;
  p->mH =  csr_free(p->mH); /** OK if `NULL` */
  p->mHt = csr_free(p->mHt);
  p->mL =  csr_free(p->mL);
  p->mLt = csr_free(p->mLt);
  p->mG = csr_free(p->mG);
}

/** @brief Check the syndrome for the LLR vector given (`0`th row of `syndrome`).
 * @return 1 if codeword is valid, 0 otherwise */
int syndrome_check(const qllr_t LLR[], const mzd_t * const syndrome,
		   const csr_t * const H,
		   [[maybe_unused]] const params_t * const p){
  const int nchk = H->rows;
  //  const csr_t * const H = p->mH;
  for(int ic=0; ic<nchk; ic++){ /** target `check` node */
    int synd = mzd_read_bit(syndrome,0,ic);
    for(int j = H->p[ic]; j < H->p[ic+1] ; j++){
      const int iv = H->i[j]; /** origin `variable` node index */
      if(LLR[iv] < 0)
	synd ^= 1;
    }
    if(synd)
      return 0; /** invalid codeword */
  }
  return 1; /** valid codeword */
}

/** @brief baseline parallel BP algorithm (sum-product or max-product depending on LLRs used).  
 * 
 * Decode `e` vector using syndrome bits in `srow = e*Ht` and apriori
 * bit LLR values in `LLR`.  Assume `srow` non-zero (no check for a
 * trivial solution).
 * 
 * @param[out] outLLR the resulting LLR (whether converged or not)
 * @return 0 if failed to converge or (the number of steps) if converged successfully
 */
int do_parallel_BP(qllr_t * outLLR, const mzd_t * const srow,
		   const csr_t * const H, const csr_t * const Ht,
		   const qllr_t LLR[], const params_t * const p){
  assert(outLLR != NULL);
  const int nchk = p->nchk; /** number of check nodes */
  const int nvar = p->nvar;     /** number of variable nodes */
  const int nz = p->mH->p[p->nchk]; /** number of messages */
  qllr_t *mesVtoC = malloc(nz*sizeof(qllr_t));
  qllr_t *mesCtoV = malloc(nz*sizeof(qllr_t));
  //  double * xLLR = malloc(nvar*sizeof(double));
  qllr_t * const xLLR = outLLR; /** return soft vector of LLR here */
  qllr_t * aLLR = calloc(nvar,sizeof(qllr_t)); /** averaged LLR */
  if((!aLLR)||(!mesVtoC)||(!mesCtoV))
    ERROR("memory allocation failed");
  int succ_BP=0;
  
    /** init V->C messages to bare LLR */
  for(int ic=0; ic<nchk; ic++){ /** target `check` node index */
    for(int j = H->p[ic]; j < H->p[ic+1] ; j++){
      const int iv = H->i[j]; /** `variable` node index */
      mesVtoC[j]= LLR[iv];   /** `j` is the edge (message) index */
    }
  }
  
  for (int istep=1; istep <= p->steps  ; istep++){ /** main decoding cycle */
    //    cnt[TOTAL]++;
    /** C -> V messages */
    for(int iv=0; iv<nvar; iv++){ /** target `variable` node index */
      for(int j = Ht->p[iv]; j < Ht->p[iv+1] ; j++){
	/** TODO: optimize this loop as in `LDPC_Code::bp_decode()` of `itpp` package */
	const int ic = Ht->i[j];/** origin `check` node index */
	int sbit = mzd_read_bit(srow,0,ic); /** syndrome bit */
	qllr_t msg=0; 
	int new=1;
	for(int j1 = H->p[ic]; j1 < H->p[ic+1] ; j1++){
	  const int iv1 = H->i[j1]; /** aux `variable` node index */
	  if(iv1!=iv){
	    if(new){
	      msg = mesVtoC[j1];
	      new=0;
	    }
	    else
	      msg = boxplus(msg, mesVtoC[j1]);
	  }
	}
	mesCtoV[j] = (sbit ? -msg : msg);
      }
    }

        
    for(int iv = 0; iv< nvar; iv++){
      qllr_t val= LLR[iv]; /** HERE! (WARNING: do we need this?) */
      for(int j1 = Ht->p[iv]; j1 < Ht->p[iv+1] ; j1++){
	// const int ic1 = Ht->i[j1];/** aux `check` node index */
	val += mesCtoV[j1];
      }
      xLLR[iv] = val;
      /** TODO: play with the decay value `0.5` */
      aLLR[iv] = llr_from_dbl(0.5 * dbl_from_llr(aLLR[iv]) + dbl_from_llr(val)); 
    }
    if(p->debug & 8){
      for(int iv = 0; iv < nvar; iv++)
	printf(" %6.2g%s", dbl_from_llr(xLLR[iv]), iv+1< nvar ? "" : "\n");
      for(int iv = 0; iv < nvar; iv++)
	printf(" %6.2g%s", dbl_from_llr(aLLR[iv]), iv+1< nvar ? "" : "\n");
    }


    if(syndrome_check(xLLR,srow,p->mH, p)){
      //      outLLR = xLLR; 
      cnt_update(CONV_BP, istep);
      succ_BP=istep;
      break;
    }
    else if(syndrome_check(aLLR,srow,p->mH, p)){
      for(int iv=0; iv< p->nvar; iv++)
	outLLR[iv] = aLLR[iv];
      cnt_update(CONV_BP_AVG, istep);
      succ_BP=-istep;
      break;
    }

    
    /* V -> C messages */ 
    for(int ic=0; ic<nchk; ic++){ /** target `check` node */
      for(int j = H->p[ic]; j < H->p[ic+1] ; j++){
	const int iv = H->i[j]; /** origin `variable` node index */
#if 0 /** TODO: verify these give identical results! */
	qllr_t msg = LLR[iv];
	for(int j1 = Ht->p[iv]; j1 < Ht->p[iv+1] ; j1++){
	  const int ic1 = Ht->i[j1];/** aux `check` node index */
	  if(ic1 != ic)
	    msg += mesCtoV[j1];
	}
	mesVtoC[j] = msg;
#else /* simplified calculation */
	//	assert(ic == Ht->i[j]);
	for(int j1 = Ht->p[iv]; j1 < Ht->p[iv+1] ; j1++){
	  const int ic1 = Ht->i[j1];/** aux `check` node index */
	  if(ic1 == ic)
	    mesVtoC[j] = xLLR[iv] - mesCtoV[j1];
	  /** TODO: use binary search here instead */
	}
#endif /* 0 */	
      }
    }

  }
  if(!succ_BP)
    outLLR = xLLR; /** TODO: make an option to change this to aLLR */
  /** clean up ********************/
  //  free(xLLR);
  if(aLLR){
    free(aLLR);
    aLLR=NULL;
  }
  free(mesVtoC);
  free(mesCtoV);

  return succ_BP;
}

/** @brief this is the actual program ******************************************/
int main(int argc, char **argv){
  params_t * const p=&prm;
  var_init(argc,argv, p); /* initialize variables, read matrices, and open file (if any) */

  /** prepare error vectors ************************************************************/
  mzd_t *mHe = mzd_init(p->nchk, p->nvec); /** each column a syndrome vector `H*e` */
  mzd_t *mLe = mzd_init(p->ncws,  p->nvec); /** each column `L*e` vector */

  if(p->internal){ /** generate errors internally */
    do_errors(mHe,mLe,p->mHt, p->mLt, p->vP);
    if(p->debug&1)
      printf("# generated %d error/obs pairs\n",mHe->ncols);		 
  }
  else{
    rci_t il1=read_01(mHe,p->file_det, &p->line_det, p->fdet, p->debug);
    /** TODO: enable external processing of observables */
    rci_t il2=read_01(mLe,p->file_obs, &p->line_obs, p->fobs, p->debug);
    if(il1!=il2)
      ERROR("mismatched DET %s (line %d) and OBS %s (line %d) files!",
	    p->fdet,p->line_det,p->fobs,p->line_obs);
    //    if(il1==0)      break; /** no more rounds */
    if(p->debug&1)
      printf("read %d error/obs pairs\n",il1);		 
  }
  mzd_t *mHeT = mzd_transpose(NULL, mHe);
  mzd_t *mLeT = mzd_transpose(NULL, mLe);
  qllr_t *ans = calloc(p->nvar, sizeof(qllr_t));
  if(!ans) ERROR("memory allocation failed!"); 
  /** counters */
  for(int ierr = 0; ierr < mHeT->nrows; ierr++){ /** cycle over errors */
    cnt[TOTAL]++;
    int succ_BP = 0;
    assert(p->LLRmin>0); /** sanity check */
    if(mzd_row_is_zero(mHeT,ierr)){
      cnt_update(CONV_TRIVIAL,0); /** trivial convergence after `0` steps */
      if(mzd_row_is_zero(mLeT,ierr)){
	cnt[SUCC_TRIVIAL]++;
	cnt[SUCC_TOT]++;
      }
      continue ; /** next error / syndrome vector pair */     
    }
    else{ /** non-trivial syndrome */
      if(p->debug&8){
	printf("non-trivial error %d:\n",ierr);
	mzd_print_row(mHeT,ierr);
	mzd_print_row(mLeT,ierr);
	for(int iv = 0; iv < p->nvar; iv++)
	  printf(" %6.2g%s", dbl_from_llr(p->vLLR[iv]), iv+1< p->nvar ? "" : "\n");
      }
      mzd_t * const srow = mzd_init_window(mHeT, ierr,0, ierr+1,p->nchk); /* syndrome row */     
      succ_BP = do_parallel_BP(ans, srow, p->mH, p->mHt, p->vLLR, p);    
      mzd_free(srow);

      if(succ_BP){/* convergence success  */
	mzd_t * const obsrow = mzd_init_window(mLeT, ierr,0, ierr+1,mLeT->ncols);
	if(syndrome_check(ans, obsrow, p->mL, p)){
	  cnt[SUCC_BP]++;
	  cnt[SUCC_TOT]++;
	}
	mzd_free(obsrow);
      }
      if(p->debug&16)
	printf("i=%d of %d succ=%d\n",ierr,mHeT->nrows,succ_BP);
    }
    if((p->nfail) && cnt[TOTAL]-cnt[SUCC_TOT] >= p->nfail)
      break;
  }
  
  cnt_out(p->debug&1);
  /** clean-up ***********************************************************************/
  if(ans){
    free(ans);
    ans=NULL;
  }
  mzd_free(mHe);
  mzd_free(mLe);
  mzd_free(mHeT);
  mzd_free(mLeT);
  
  var_kill(p);
  return 0;
}
