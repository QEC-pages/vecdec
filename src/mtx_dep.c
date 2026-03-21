/**
 *  @file mtx_dep.c
 *
 * @brief mtx_dep - generate depolarizing noise matrices
 *
 * @author Leonid Pryadko (Google Quantum AI)
 *
 * Copyright (C) 2025 Leonid Pryadko
 * Google Quantum AI
 * All rights reserved.
 *
 *    z  x  y 
 * H=[Hx    Hx]
 *   [   Hz Hz]   rk H = rX+rZ = n - k
 *
 * L=[Lx    Lx]   rk L = k
 *
 * G=[   Hx   ]   rk G = n+rX+rZ+k = 2n
 *   [Hz      ]
 *   [I  I  I ]
 *   [   Lx   ]   tot  = 3*n, correct
 *
 * K=[Lz      ]   rk K = k 
 */

#include "utils.h"
#include "util_m4ri.h"

typedef struct PAR_T {
  char *finH, *finG, *finL, *out;
  int rH,rG,k,n;/** input matrices: H[rH,n], G[rG,n], L[k,n], K[k,n] */
  int rkH,rkG; /** computed */
  csr_t *matH, *matL, *matG, *matHH, *matLL;
  int debug;
} par_t;

#define USAGE								\
  "%s: generate repeated depolarizing code matrices\n"			\
  "  usage: %s param=value [[param=value] ... ]\n"                      \
  "\t Command line arguments are processed in the order given.\n"       \
  "\t Supported parameters:\n"                                          \
  "\t --help\t\t: give this help (also '-h' or just 'help')\n"          \
  "\t --morehelp\t\t: illustrate matrices created (or just 'morehelp')\n" \
  "\t finH=[string]\t: name of the input file with the H matrix\n"	\
  "\t finL=[string]\t: name of the input file with the L matrix (mm or alist)\n" \
  "\t finG=[string]\t: name of the input file with the H matrix\n"	\
  "\t finK=[string]\t: name of the input file with the L matrix (mm or alist)\n" \
  "\t out=[string]\t: base name of the output files (or 'stdout' by default)\n" \
  "\t debug=[integer]\t: bitmap for aux information to output (default: 3)\n" \
  "\t .1 (bit 0) output helpful information, e.g. ranks and code params\n" \
  "\t .2 (bit 1) output skeleton map of the code\n"                     \
  "\t .4 (bit 2) debug input parameters processing\n"                   \
  "\t .8 (bit 3) print out the input matrices\n"                        \
  "\t .16 (bit 4) print out the constructed matrices\n"                 \
  "\t Multiple 'debug' parameters are XOR combined except for 0.\n"     \
  "\t The program can take sparse binary matrices in MatrixMarket or alist formats.\n" \
  "\t Space is allowed between 'finH=' / 'finL=' /'out=' and the file name string.\n"

#define SAMPLE							\
  "%s: generate depolarizing noise code matrices\n"             \
  "\n"								\
  "\t\t     z x y \n"                                           \
  "\t\t HH=[H   H]\n"                                           \
  "\t\t    [  G G]   rk HH = rk H + rk G = n - k\n"             \
  "\t\t\n"                                                      \
  "\t\t LL=[L   L]   rk LL = k\n"                               \
  "\t\t\n"                                                      \
  "\t\t GG=[  H  ]   rk GG = n + rk H + rk G + k = 2*n\n"       \
  "\t\t    [G    ]\n"                                           \
  "\t\t    [I I I]\n"                                           \
  "\t\t    [  L  ]   total  = 3*n\n"                            \
  "\t\t\n"                                                      \
  "\t\t KK=[K    ]   rk KK = k \n"                              \
  "\n"                

par_t params={
  .finH=NULL,
  .finG=NULL,
  .finL=NULL,
  .out="std" "out",
  .matH=NULL,
  .matG=NULL,
  .matL=NULL,
  .matHH=NULL,
  .matLL=NULL,
  .rH=0,
  .rG=0,
  .k=0,
  .n=0,
  .rkH=0,
  .rkG=0,
  .debug=3
};

int var_init(int argc, char **argv, par_t *p){
  int cnt=0;

  if (argc==1)
    ERROR("try '%s --help'\n",argv[0]);

  for(int i=1; i<argc; i++){
    int dbg=0;
    if(sscanf(argv[i],"debug=%d",& dbg)==1){/** `debug` */
      cnt++;
      if(cnt==1)
	p->debug = dbg; /** just assign if encountered once */
      else
	p->debug ^= dbg; /** otherwise `XOR` */
      if(p->debug &4)
	printf("# read %s, debug=%d octal=%o\n",argv[i],p->debug,p->debug);
      
    }
    //    printf("# i=%d cnt=%d debug=%d octal=%o\n",i,cnt,p->debug,p->debug);

  };

  for(int i=1; i<argc; i++){   /** remaining arguments */
    int dbg;
    if(sscanf(argv[i],"debug=%d",& dbg)==1){
      /** do nothing, already processed all `debug` entries */
    }
    else if (0==strncmp(argv[i],"out=",4)){ /** `out` */
      if(strlen(argv[i])>4)
	p->out = argv[i]+4;
      else
	p->out = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, out=%s\n",argv[i],p->out);
    }
    else if (0==strncmp(argv[i],"finH=",5)){ /** `finH` */
      if(strlen(argv[i])>5)
        p->finH = argv[i]+5;
      else
        p->finH = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finH=%s\n",argv[i],p->finH);
      p->matH=csr_mm_read(p->finH, p->matH, 0, p->debug);
      p->n=p->matH->cols;
      p->rH=p->matH->rows;
      if(p->debug&1){
	p-> rkH = rank_csr(p->matH);
	printf("# read matrix H %d by %d , rank=%d\n", p->rH,p->n, p->rkH);
	if(p->debug&8){
	  if (p->n < 80){
	    printf("# matrix H:\n");
	    mzd_t *mmat = mzd_from_csr(NULL,p->matH);
	    mzd_print(mmat);
	    mzd_free(mmat);
	  }
	  else
	    csr_out(p->matH);
	}
      }
    }
    else if (0==strncmp(argv[i],"finG=",5)){ /** `finG` */
      if(strlen(argv[i])>5)
        p->finG = argv[i]+5;
      else
        p->finG = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finG=%s\n",argv[i],p->finG);
      p->matG=csr_mm_read(p->finG, p->matG, 0, p->debug);
      p->rG=p->matG->rows;
      if((p->n)&&(p->n != p->matG->cols))
        ERROR("mismatched matrices: n=%d G[%d,%d]\n",p->n,p->rG,p->matG->cols);
      if(p->debug&1){
	p-> rkG = rank_csr(p->matG);
	printf("# read matrix G %d by %d , rank=%d\n", p->rG,p->n, p->rkG);
	if(p->debug&8){
	  if (p->n < 80){
	    printf("# matrix G:\n");
	    mzd_t *mmat = mzd_from_csr(NULL,p->matG);
	    mzd_print(mmat);
	    mzd_free(mmat);
	  }
	  else
	    csr_out(p->matG);
	}
      }
    }
    else if (0==strncmp(argv[i],"finL=",5)){ /** `finL` */
      if(strlen(argv[i])>5)
        p->finL = argv[i]+5;
      else
        p->finL = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finL=%s\n",argv[i],p->finL);
      p->matL=csr_mm_read(p->finL, p->matL, 0, p->debug);
      p->n=p->matL->cols;
      p->k=p->matL->rows;
      int rkL=rank_csr(p->matL);
      if(p->debug&1){
	printf("# read matrix L %d by %d , rank=%d\n",
	       p->k,p->matL->cols, rkL);
	if(p->debug&8){
	  if (p->n < 80){
	    printf("# matrix L:\n");
	    mzd_t *mmat = mzd_from_csr(NULL,p->matL);
	    mzd_print(mmat);
	    mzd_free(mmat);
	  }
	  else
	    csr_out(p->matL);
	}
      }
      if (rkL != p->k)
	ERROR("invalid L matrix : mismatch rows=%d vs rank=%d\n",p->k,rkL);
    }
    else if((strcmp(argv[i],"--help")==0)
            ||(strcmp(argv[i],"-h")==0)
            ||(strcmp(argv[i],"help")==0)){
      printf( USAGE , argv[0],argv[0]);
      exit (-1);
    }
    else if((strcmp(argv[i],"morehelp")==0)
            ||(strcmp(argv[i],"-hm")==0)
            ||(strcmp(argv[i],"--morehelp")==0)){
      printf( SAMPLE , argv[0]);
      exit (-1);
    }
    else{
      printf("unrecognized parameter argv[%d]=%s\n",i,argv[i]);
      ERROR("try '%s --help'\n",argv[0]);
    }
  }

  if ((p->matL )&& (p->n != p->matL->cols))
    ERROR("matrix dimension mismatch: L %d x %d, n=%d", p->matL->rows,p->matL->cols,p->n);
  if ((p->matH )&& (p->n != p->matH->cols))
    ERROR("matrix dimension mismatch: H %d x %d, n=%d", p->matH->rows,p->matH->cols,p->n);
  if ((p->matG )&& (p->n != p->matG->cols))
    ERROR("matrix dimension mismatch: G %d x %d, n=%d", p->matG->rows,p->matG->cols,p->n);
  if((p->matG==NULL) || (p->matH==NULL))
    ERROR("need both H=Hx and G=Hz matrices");
  return 0;
}


int main(int argc, char **argv){
  
  par_t * const p = &params;
  var_init(argc,argv, p);

  if(p->debug&3)
    printf("# constructing HH\n");
  if(p->debug&2){
    printf("#\t\t     z x y \n"
           "#\t\t HH=[H   H]\n" 
           "#\t\t    [  G G]   %d by %d, rk HH = %d \n",
           p->rH+p->rG,3*p->n,p->rkH + p->rkG);

  }
  /** actually construct HH */
  int nzH=p->matH->p[p->matH->rows];
  int nzG=p->matG->p[p->matG->rows];
  if(p->debug&1){
      printf("# H:  rows=%d cols=%d nz=%d\n",p->matH->rows, p->matH->cols, nzH);
      printf("# G:  rows=%d cols=%d nz=%d\n",p->matG->rows, p->matG->cols, nzG);
  }      
  int nzHH = 2*nzH + 2*nzG;
  int rHH = p->rH + p->rG;
  int cHH = 3*p->n;
  if(p->debug&1)
    printf("HH: rows=%d cols=%d nz=%d\n",rHH, cHH, nzHH);
    
  p->matHH = csr_init(NULL, rHH, cHH, nzHH);
  int gr=0; /** global row */
  int ge=0;     /** global element index */
  for (int br=0; br< 2 ; br++){ /** block row index  */
    const csr_t *const mat = br == 0 ? p->matH : p->matG ;        
    for(int j=0; j< mat->rows; j++, gr++){ /** rows of `H` or `G` */
      p->matHH->p[gr]=ge;
      int idx_beg=mat->p[j], idx_max=mat->p[j+1];	
      for (int bc=0; bc < 2; bc++){ /** block columns */
        int cbeg=0;
        if(bc)
          cbeg=2*p->n;  /** last column */
        else
          cbeg=p->n * br; /** `H0` or `0G` */
        for(int idx=idx_beg; idx<idx_max; idx++) /** insert `H` or `G` */
          p->matHH->i[ge++] = cbeg + mat->i[idx];
      }
    }
    p->matHH->p[gr]=ge;     
  }
  p->matHH->nz=-1;     
  if(ge != nzHH)
    ERROR("this should not happen nzHH=%d expected %d ",ge, nzHH);
   

  if(p->finL){
    if(p->debug&3)
      printf("# constructing LL\n");
    if(p->debug&2) /** skeleton */
      printf("# LL = [ L . L ]\n");    
    /** actually construct LL */
    int nzL=p->matL->p[p->matL->rows];
    if(p->debug&1)
      printf("L:  rows=%d cols=%d nz=%d\n",p->matL->rows, p->matL->cols, nzL);
    int nzLL = 2 * nzL ;
    int rLL = p->matL->rows;
    int cLL = 3*p->n;
    if(p->debug&1)
      printf("LL: rows=%d cols=%d nz=%d\n",rLL, cLL, nzLL);

    p->matLL = csr_init(NULL, rLL, cLL, nzLL);
    int gr=0; /** global row */
    int ge=0;     /** global element index */
    for(gr=0; gr< p->matL->rows; gr++){/** rows of `L` */
      p->matLL->p[gr]=ge;
      int idx_beg=p->matL->p[gr], idx_max=p->matL->p[gr+1];	
      for (int bc=0; bc < 3 ; bc += 2){ /** block columns */
	int cbeg = bc * p->n; 	  
	for(int idx=idx_beg; idx<idx_max; idx++) /** insert `L` */
	  p->matLL->i[ge++] = cbeg + p->matL->i[idx];
      }
    }
    p->matLL->p[gr]=ge;         
    p->matLL->nz=-1;     
    if(ge != nzLL)
      ERROR("this should not happen nzLL=%d expected %d ",ge, nzLL);
  }

  const size_t siz=1000;
  char comment[siz+1], fnam[siz+1], *s;
  
  if(p->matHH){
    snprintf(comment,siz,"depol block matrix from H=%s G=%s",p->finH,p->finG);
    comment[siz]='\0';/** sanity check */
    if(strcmp(p->out,"stdout")!=0){    
      snprintf(fnam,siz,"%s_HH.mmx",p->out);
      fnam[siz]='\0';/** sanity check */
      s=fnam;
      if(p->debug&1)
	printf("# writing HH matrix to %s\n",s);
    }
    else
      s=p->out;
    csr_mm_write(s,"" /** no extension */, p->matHH, comment);

    if(p->debug&8){
      if(p->matHH->cols<80){
	mzd_t *dmat = mzd_from_csr(NULL,p->matHH);
	mzd_print(dmat);
	mzd_free(dmat);
      }
      else
	csr_out(p->matHH);    
    }
    csr_free(p->matHH);
  }

  if(p->matLL){
    snprintf(comment,siz,"depol block matrix from L=%s",p->finL);
    comment[siz]='\0';/** sanity check */
    if(strcmp(p->out,"stdout")!=0){    
      snprintf(fnam,siz,"%s_LL.mmx",p->out);
      fnam[siz]='\0';/** sanity check */
      s=fnam;
      if(p->debug&1)
	printf("# writing LL matrix to %s\n",s);
    }
    else
      s=p->out;
    csr_mm_write(s,"" /** no extension */, p->matLL, comment);

    if(p->debug&8){
      if(p->matLL->cols<80){
	mzd_t *dmat = mzd_from_csr(NULL,p->matLL);
	mzd_print(dmat);
	mzd_free(dmat);
      }
      else
	csr_out(p->matLL);    
    }
    csr_free(p->matLL);
  }
  if(p->matH) csr_free(p->matH);
  if(p->matG) csr_free(p->matG);
  if(p->matL) csr_free(p->matL);
  return 0;
}
