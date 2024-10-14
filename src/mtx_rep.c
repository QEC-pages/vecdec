/**
 *  @file mtx_rep.c
 *
 * @brief mtx_rep - generate repeated measurement matrices
 *
 * @author Leonid Pryadko (University of California, Riverside)
 *
 * Copyright (C) 2024 Leonid Pryadko
 * University of California, Riverside
 * All rights reserved.
 *
 *	    num=1: HH=[H I] 
 *	           LL=[L 0]
 *
 *	    num=2: HH=[H I 0]  
 *	              [0 I H]  
 *                 LL=[L 0 L]               
 *	   
 *	    num=3: HH=[H I      ]  
 *    	              [  I H I  ]
 *  	   	      [      I H]
 *                 LL=[L 0 L 0 L]
 */

#include "utils.h"
#include "util_m4ri.h"

typedef struct PAR_T {
  char *finH, *finL, *out;
  int num;  /* number of repeated measurements. */	   
  int r,k,n;/** dimensions of the input matrices */
  int rkH; /** computed */
  int rr, nn; /** dimensions of the output matrix */
  csr_t *matH, *matL, *matHH, *matLL;
  int debug;
} par_t;

#define USAGE								\
  "%s: generate repeated measurement matrices\n"			\
  "  usage: %s param=value [[param=value] ... ]\n"                      \
  "\t Command line arguments are processed in the order given.\n"       \
  "\t Supported parameters:\n"                                          \
  "\t --help\t\t: give this help (also '-h' or just 'help')\n"          \
  "\t --morehelp\t\t: illustrate matrices created (or just 'morehelp')\n" \
  "\t finH=[string]\t: name of the input file with the H matrix\n"	\
  "\t finL=[string]\t: name of the input file with the L matrix (mm or alist)\n" \
  "\t out=[string]\t: base name of the output files (or 'stdout' by default)\n" \
  "\t num=[int]\t: number of repeated measurements (0: original matrices).\n" \
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
  "%s: generate repeated measurement matrices\n"		\
  "\n"								\
  "\t\t     num=0: HH=[H]  LL=[L] \n"				\
  "\n"								\
  "\t\t     num=1: HH=[H I 0]  LL=[L 0 L] \n"                   \
  "\t\t               [0 I H]             \n"                

par_t params={
  .finH=NULL,
  .finL=NULL,
  .out="std" "out",
  .num=1,
  .matH=NULL,
  .matL=NULL,
  .matHH=NULL,
  .matLL=NULL,
  .r=0,
  .k=0,
  .n=0,
  .rr=0,
  .nn=0,
  .rkH=0,
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
      p->r=p->matH->rows;
      if(p->debug&1){
	p-> rkH = rank_csr(p->matH);
	printf("# read matrix H %d by %d , rank=%d\n", p->r,p->n, p->rkH);
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
    else if  (sscanf(argv[i],"num=%d",&dbg)==1){ /** `num` param */
      p -> num = dbg;
      if (p->debug&4)
	printf("# read %s, setting num=%d\n",argv[i],p-> num);
      if(p->num<=0)
	ERROR("must be positive num=%d\n",p->num);
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

  if (((p->matL )&& (p->n != p->matL->cols)) ||
      ((p->matH )&& (p->n != p->matH->cols)))
	ERROR("matrix dimension mismatch: H %d x %d, L %d x %d",
	      p->matH->rows,p->matH->cols,p->matL->rows,p->matL->cols);

  if((p->matL) && (p->matH==NULL))
    ERROR("can't use L without H matrix");
  return 0;
}


int main(int argc, char **argv){
  
  par_t * const p = &params;
  var_init(argc,argv, p);

  //  if(p->debug&1)
  //    printf("# generating %d x %d QC matrix from %d x %d blocks of size ell=%d nz=%d\n",
  //	   rc,n,p->rows,p->cols,p->ell,nz);


  if(p->finH){
    if(p->debug&3)
      printf("# constructing HH\n");
    if(p->debug&2){
      for (int br=0; br< p->num; br++){ /** block row index  */
	printf("#  ");
	for (int bc=0; bc < p->num; bc++){ /** block columns */
	  if(bc<br)
	    printf("%s. ", bc==0? ""  : "* ");
	  else if(br==bc)
	    printf("%sH ", bc==0? ""  : "I ");
	  else if (br+1 == bc)
	    printf("I . ");
	  else
	    printf("* . ");
	}
	printf("  #\n");
      }
    }
    /** actually construct HH */
    int nzH=p->matH->p[p->matH->rows];
    if(p->debug&1)
      printf("H:  rows=%d cols=%d nz=%d\n",p->matH->rows, p->matH->cols, nzH);
    int nzHH = p->num*nzH + 2*p->matH->rows*(p->num-1);
    int rows = p->num*p->matH->rows;
    int cols = (p->num-1)*p->matH->rows + p->num*p->matH->cols;
    if(p->debug&1)
      printf("HH: rows=%d cols=%d nz=%d\n",rows, cols, nzHH);
    
    p->matHH = csr_init(NULL, rows, cols, nzHH);
    int gr=0; /** global row */
    int ge=0;     /** global element index */
    for (int br=0; br< p->num; br++){ /** block row index  */
      for(int j=0; j< p->matH->rows; j++, gr++){/** rows of `H` */
	if(gr >= p->matHH->rows)
	  ERROR("this should not happen row=%d >= nrows=%d",gr, p->matHH->rows);	  
	int cbeg=0; /** starting column of this block */
	p->matHH->p[gr]=ge;
	int idx_beg=p->matH->p[j], idx_max=p->matH->p[j+1];	
	for (int bc=0; bc < p->num ; bc++){ /** block columns */
	  if(bc<br){
	    cbeg += (bc!=0 ? p->matH->rows : 0) + p->matH->cols;
	  }
	  else if(br==bc){
	    if(bc!=0){ /** insert `I` */
	      p->matHH->i[ge++] = cbeg + j;
	      if(cbeg+j>=p->matHH->cols)
		ERROR("this should not happen col=%d >= ncols=%d",cbeg+j, p->matHH->cols);
	      cbeg += p->matH->rows; 
	    }
	    for(int idx=idx_beg; idx<idx_max; idx++){ /** insert `H` */
	      p->matHH->i[ge++] = cbeg + p->matH->i[idx];
	      if(cbeg+p->matH->i[idx]>=p->matHH->cols)
		ERROR("this should not happen col=%d >= ncols=%d",cbeg+p->matH->i[idx], p->matHH->cols);
	    }
	    cbeg += p->matH->cols; 
	  }
	  else if (br+1 == bc){
	    p->matHH->i[ge++] = cbeg + j;
	    if(cbeg+j>=p->matHH->cols)
	      ERROR("this should not happen col=%d >= ncols=%d",cbeg+j, p->matHH->cols);	    
	  }
	}
      }
      p->matHH->p[gr]=ge;     
    }
    p->matHH->nz=-1;     
    if(ge != nzHH)
       ERROR("this should not happen nzHH=%d expected %d ",ge, nzHH);
   
  }
  if(p->finL){
    if(p->debug&3)
      printf("# constructing LL\n");
    if(p->debug&2){ /** skeleton */
      for (int bc=0; bc < p->num ; bc++) /** block columns */
	printf("%sL %s", bc==0? "#  ":"", bc+1< p->num? ". ":"  #\n");
    }
    /** actually construct LL */
    int nzL=p->matL->p[p->matL->rows];
    if(p->debug&1)
      printf("L:  rows=%d cols=%d nz=%d\n",p->matL->rows, p->matL->cols, nzL);
    int nzLL = p->num * nzL ;
    int rows = p->matL->rows;
    int cols = (p->num-1)*p->matH->rows + p->num*p->matH->cols;
    if(p->debug&1)
      printf("LL: rows=%d cols=%d nz=%d\n",rows, cols, nzLL);

    p->matLL = csr_init(NULL, rows, cols, nzLL);
    int gr=0; /** global row */
    int ge=0;     /** global element index */
    for(gr=0; gr< p->matL->rows; gr++){/** rows of `L` */
      if(gr >= p->matLL->rows)
	ERROR("this should not happen row=%d >= nrows=%d",gr, p->matLL->rows);	  
      int cbeg=0; /** starting column of this block */
      p->matLL->p[gr]=ge;
      int idx_beg=p->matL->p[gr], idx_max=p->matL->p[gr+1];	
      for (int bc=0; bc < p->num ; bc++){ /** block columns */
	for(int idx=idx_beg; idx<idx_max; idx++){ /** insert `H` */
	  p->matLL->i[ge++] = cbeg + p->matL->i[idx];
	  if(cbeg+p->matL->i[idx]>=p->matLL->cols)
	    ERROR("this should not happen col=%d >= ncols=%d",cbeg+p->matL->i[idx], p->matLL->cols);
	}
	cbeg += p->matL->cols; 	  
	if(bc!=p->num){ /** insert empty `*` matrix */
	  cbeg += p->matH->rows;
	}	  
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
    snprintf(comment,siz,"rep block matrix num=%d from H=%s",p->num,p->finH);
    comment[siz]='\0';/** sanity check */
    if(strcmp(p->out,"stdout")!=0){    
      snprintf(fnam,siz,"%s_%dH.mmx",p->out,p->num);
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
    snprintf(comment,siz,"rep block matrix num=%d from L=%s",p->num,p->finL);
    comment[siz]='\0';/** sanity check */
    if(strcmp(p->out,"stdout")!=0){    
      snprintf(fnam,siz,"%s_%dL.mmx",p->out,p->num);
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
  if(p->matL) csr_free(p->matL);
  return 0;
}
