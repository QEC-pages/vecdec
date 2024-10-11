/**
 *  @file mtx_sub.c
 *
 * @brief mtx_qc - generate a quasicyclic matrix from polynomials
 *
 * @author Leonid Pryadko (University of California, Riverside)
 *
 * Copyright (C) 2024 Leonid Pryadko
 * University of California, Riverside
 * All rights reserved.
 *
 *
 */

#include "utils.h"
#include "util_m4ri.h"

#define MAX_ROWS 2
#define MAX_COLS 5
#define MAX_W 10
#define WGHT(i,j) p->wght[p->rows*(i)+(j)]
#define COEF(i,j,k) p->poly[MAX_W*(p->cols*(i)+(j))+k]

typedef struct PAR_T {
  int rows; /** how many circulant rows, default 1 */
  int cols; /** how many circulant cols, default 1 */
  int ell; /** size of a circulant (must be set), default `0` */
  int r[MAX_ROWS]; /** num rows in each row block (values can't exceed `n`) */
  int poly[(MAX_ROWS)*(MAX_COLS)*(MAX_W)]; /** coefficients of polynomials */
  int wght[(MAX_ROWS)*(MAX_COLS)]; /** polynomial weight (number of non-zero terms) */
  char *out; /** name of the output file */
  int tot_cols,tot_rows; /** actual size of the matrix: `n=cols*ell`, `r=r[0]+r[1]+... r[rows]` */
  csr_t *mat; /** constructed matrix */
  int debug;
} par_t;

par_t params={
  .rows=1,
  .cols=0,  //! must specify 
  .ell=0,   //! must specify
  .out="std" "out",
  .mat=NULL,
  .debug=1
};

par_t *p = &params;

#define USAGE							       \
  "%s: mtx_qc - generate a binary quasicyclic matrix from polynomials\n" \
  "  usage: %s param=value [[param=value] ... ]\n"			\
  "\t Parameters are read in the order specified.\n"			\
  "\t Supported parameters:\n"						\
  "\t --help\t\t: give this help (also '-h' or just 'help')\n"		\
  "\t out=[string]\t: name of the output file (or 'stdout', default value)\n" \
  "\t rows=int\t: number of circulant rows (default: 1)\n"		\
  "\t cols=int\t: number of circulant columns (no default)\n" \
  "\t ell=int\t: cize of a circulant (no default)\n"			\
  "\t a_#=[int,int,...]\t degrees of 0th row polynomial in column '#' (same as 'a0_#')\n" \
  "\t a#_#=[int,int,...]\t degrees of polynomial (row,col)\n" \
  "\t r#=int\t: number of rows in row block '#' (cannot exceed ell, the default value)\n" \
  "\t\t Specify polynomials only for non-zero blocks\n" \
  "\t\t List polynomial deegrees in increasing order from 0 to (ell-1)\n" \
  "\t debug=[integer]\t: bitmap for aux information to output (default: 1)\n" \
  "\t Multiple 'debug' parameters are XOR combined except for 0.\n"	\
  "\t Space is allowed between 'out=' and the file name string.\n"	\
  "\t Output file name 'out' must be the last argument.\n"

int var_init(int argc, char **argv, par_t *p){
  
  memset(p->r,    '\0', sizeof(p->r));
  memset(p->poly, '\0', sizeof(p->poly));
  memset(p->wght, '\0', sizeof(p->wght));

  for(int i=1; i<argc; i++){   /** scan for debug first */
    int dbg;
    if(sscanf(argv[i],"debug=%d",& dbg)==1){/** `debug` */
      if(dbg==0)
	p->debug = 0;
      else{
	if(i==1)
	  p->debug = dbg; /** just assign if in the `1st position` */
	else
	  p->debug ^= dbg; /** otherwise `XOR` */
	if(p->debug &4)
	  printf("# read %s, debug=%d octal=%o\n",argv[i],p->debug,p->debug);
      }
    }
  }

  for(int i=1; i<argc; i++){   /** scan for rows, cols, ell second */
    int dbg;
    if (sscanf(argv[i],"ell=%d",& dbg)==1){ /** `ell` */
      p -> ell = dbg;
      if (p->debug&4)
	printf("# read %s, setting ell=%d\n",argv[i],p->ell);
    }
    else if (sscanf(argv[i],"rows=%d",& dbg)==1){ /** `rows` */
      if(dbg>=MAX_ROWS)
	ERROR("parameter rows=%d must be < MAX_ROWS=%d",dbg,MAX_ROWS);
      p -> rows = dbg;
      if (p->debug&4)
	printf("# read %s, setting rows=%d\n",argv[i],p->rows);
    }
    else if (sscanf(argv[i],"cols=%d",& dbg)==1){ /** `cols` */
      if(dbg>=MAX_COLS)
	ERROR("parameter cols=%d must be < MAX_COLS=%d",dbg,MAX_COLS);
      p -> cols = dbg;
      if (p->debug&4)
	printf("# read %s, setting cols=%d\n",argv[i],p->cols);
    }
  }
  
  if((p->ell==0) || (p->rows==0) || (p->cols==0))
    ERROR("must specify non-zero block rows=%d block cols=%d and"
	  " circulant size ell=%d",p->rows,p->cols,p->ell);

  for(int i=0; i< p->rows; i++)
    p->r[i]=p->ell; /** default values */
  
  for(int i=1; i<argc; i++){   /** remaining arguments */
    int pos, dbg, p1=0, p2=0;
    if((sscanf(argv[i],"debug=%d",& dbg)==1) ||
       (sscanf(argv[i],"rows=%d",& dbg)==1)  ||
       (sscanf(argv[i],"cols=%d",& dbg)==1)  ||
       (sscanf(argv[i],"ell=%d",& dbg)==1)){
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
    else if (sscanf(argv[i],"r%d=%d",&p1, & dbg)==2){ /** `r#=` */
      if((p1<0) || (p1>=p->rows))
	ERROR("arg[%d]='%s' : invalid row index %d, must be from 0 to rows-1 = %d\n",
	      i,argv[i],p1,p->rows -1);
      if((dbg<0) || (dbg>= p->ell))
	ERROR("arg[%d]='%s' invalid row num=%d in row block %d, "
	      "must be from 1 to ell-1 = %d",
	      i,argv[i],dbg,p1,p->ell -1);
      p -> r[p1] = dbg;
      if (p->debug&4)
	printf("# read %s, setting r[%d]=%d\n",argv[i],p1,p->r[p1]);
    }
    else if ((sscanf(argv[i],"a_%d=%d%n",&p2, & dbg, &pos)==2)||
	     (sscanf(argv[i],"a%d_%d=%d%n",&p1, &p2, &dbg, &pos)==3)){
      /* a[p1,p2]=... */
      int w=0;
      char *c = argv[i];
      if((p1<0)||(p1>=p->rows)||(p2<0)||(p2>=p->cols))
	ERROR("arg[%d]='%s' : a(%d,%d)=... invalid position rows=%d cols=%d\n",
	      i,argv[i],p1,p2,p->rows,p->cols);
      do{
	if((dbg<0)||(dbg>=p->ell))
	  ERROR("arg[%d]=%s : invalid degree(%d)=%d , ell=%d\n",
		i,argv[i],w,dbg,p->ell);
	if((w>0)&&(dbg<=COEF(p1,p2,w-1)))
	  ERROR("arg[%d]='%s' deg(%d)=%d deg(%d)=%d :"
		" must be an increasing sequence! ",
		i,argv[i],w-1,COEF(p1,p2,w-1),w,dbg);

	COEF(p1,p2,w++)=dbg;
	c+=pos;
      }
      while (sscanf(c,",%d%n",& dbg, &pos)==1);
	
      WGHT(p1,p2)=w;
      if(p->debug&4){
	printf("# read polynomial a(%d,%d)=[ ",p1,p2);
	for(int ii=0; ii< w; ii++)
	  printf("%s%d%s", COEF(p1,p2,ii)==0 ? "":"x^",
		 COEF(p1,p2,ii)==0 ? 1:COEF(p1,p2,ii),
		 ii+1==w ? " ] ":", ");
	printf(" wght=%d\n",w);
      }
    }
  }
  return 0;
}

int main(int argc, char **argv){
  var_init(argc,argv, p);

  int n = p->ell * p->cols; /* number of columns */
  int rc =0; /* row count */
  int nz=0; /* non-zero elements */
  for(int i=0; i< p->rows; i++){
    int ri= p->r[i]; /* rows in this block */
    rc += ri;
    for(int j=0; j< p->cols; j++){
      nz+= WGHT(i,j) * ri;
    }
  }
  p->mat = csr_init(NULL, rc, n, nz);
  int gr=0; /** global row */
  int ge=0;     /** global element index */
  for (int br=0; br< p->rows; br++){ /** block row index  */
    for(int ii=0; ii< p->r[br]; ii++, gr++){ /** rows in block */
      int cbeg=0; /** start column of this block */
      for(int bc=0; bc < p->cols; bc++, cbeg+=p->ell){ /** block columns */
	int wei=WGHT(br,bc);
	for(int idx=0; idx<wei; idx++){
	  int i_shift = (ii+wei-idx)%wei; /* shifted element index */
	  //	  int c_shift = (ii+p->ell+COEF(br,bc,i_shift))%(p->ell);
	  int c_shift = (ii + p->ell - COEF(br,bc,i_shift)) % (p->ell);
	  if(p->debug&32)
	    printf("br=%d ii=%d bc=%d idx=%d i_s=%d c_beg=%d c_s=%d  ge=%d COEF=%d\n",
		   br,ii,bc,idx,i_shift,cbeg,c_shift,ge,COEF(br,bc,i_shift));
	  /** write a list of pairs */
	  p->mat->p[ge] = gr;
	  p->mat->i[ge++] = cbeg + c_shift;
	}
      }
    }
  }
  if((nz!=ge)||(gr!=rc))
    ERROR("this should not happen: nz=%d ge=%d  rc=%d gr=%d\n",nz,ge,rc,gr);
  p->mat->nz=nz; /* non-zero elements */

  if(p->debug&64)
    csr_out(p->mat);
  csr_compress(p->mat);


  if(p->debug&128){
    mzd_t *dmat = mzd_from_csr(NULL,p->mat);
    mzd_print(dmat);
    mzd_free(dmat);
    
  }

  const size_t siz=1000;
  char comment[siz+1];
  snprintf(comment,siz,"QC block matrix [%d,%d] ell=%d\n",p->rows,p->cols, p->ell);
  comment[siz]='\0';/** sanity check */
  csr_mm_write(p->out,"" /** no extension */, p->mat, comment);
  if(p->mat) csr_free(p->mat);
  return 0;
}

      

  
