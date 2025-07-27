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

#define MAX_ROWS 20
#define MAX_COLS 20
#define MAX_W 40
#define BLOCK_INDEX(i,j) (p->cols*(i)+(j)) //! row=i, col=j 
#define WGHT(i,j) p->wght[BLOCK_INDEX(i,j)]
#define COEF(i,j,k) p->poly[(MAX_W)*BLOCK_INDEX(i,j)+k]

typedef struct PAR_T {
  int rows; /** how many circulant rows, default 1 */
  int cols; /** how many circulant cols, default 1 */
  int ell; /** size of a circulant (must be set), default `0` */
  int r[MAX_ROWS]; /** num rows in each row block (values can't exceed `n`) */
  int c[MAX_COLS]; /** num cols in each column block (values can't exceed `n`) */
  int poly[(MAX_ROWS)*(MAX_COLS)*(MAX_W)]; /** coefficients of polynomials */
  int wght[(MAX_ROWS)*(MAX_COLS)]; /** polynomial weight (number of non-zero terms) */
  char *out; /** name of the output file */
  int tot_cols,tot_rows; /** actual size of the matrix: `n=c[0]+c[1]+... +c[cols]`, `r=r[0]+r[1]+... r[rows]` */
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
  "\t b_#=... b#_#=...\t same but also revert the polynomial\n" \
  "\t r#=int\t: number of rows in row block '#' (cannot exceed ell, the default value)\n" \
  "\t c#=int\t: number of cols in col block '#' (cannot exceed ell, the default value)\n" \
  "\t\t Specify polynomials only for non-zero blocks\n" \
  "\t\t List polynomial deegrees in increasing order from 0 to (ell-1)\n" \
  "\t debug=[integer]\t: bitmap for aux information to output (default: 1)\n" \
  "\t .1 (bit 0) output helpful information, e.g. ranks and code params\n" \
  "\t .2 (bit 1) output skeleton map of the code\n"                     \
  "\t .4 (bit 2) echo input parameters being read\n"                   \
  "\t .8 (bit 3) print out the constructed matrix\n"			\
  "\t .16 (bit 4) extra processing information\n"			\
  "\t Multiple 'debug' parameters are XOR combined except for 0.\n"	\
  "\t Space is allowed between 'out=' and the file name string.\n"

void do_scan_a(int i, const char argvi[], int p1, int p2, int val, int pos, int reverse){

  int w=0, dbg=val;
  const char *c = argvi;
  if((p1<0)||(p1>=p->rows)||(p2<0)||(p2>=p->cols))
    ERROR("arg[%d]='%s' : a(%d,%d)=... invalid position rows=%d cols=%d\n",
	  i,argvi,p1,p2,p->rows,p->cols);
  do{
    if((dbg<0)||(dbg>=p->ell))
      ERROR("arg[%d]='%s' : invalid degree(%d)=%d , ell=%d\n",
	    i,argvi,w,dbg,p->ell);
    if((w>0)&&(dbg<=COEF(p1,p2,w-1)))
      ERROR("arg[%d]='%s' deg(%d)=%d deg(%d)=%d :"
	    " must be an increasing sequence! ",
	    i,argvi,w-1,COEF(p1,p2,w-1),w,dbg);
    COEF(p1,p2,w++)=dbg;
    if(w>MAX_W)
      ERROR("arg[%d]='%s' : too large polynomial weight>=%d MAX_W=%d\n", i,argvi,w,MAX_W);
    c+=pos;
  }
  while (sscanf(c,",%d%n",& dbg, &pos)==1);

  if (reverse)
    for(int i=0; i< w; i++)
      COEF(p1,p2,i) = (p->ell - COEF(p1,p2,i)) % p->ell;        
	
  WGHT(p1,p2)=w;
  if(p->debug&4){
    printf("# read %spolynomial a(%d,%d)=[ ",reverse? "reverse ":"",p1,p2);
    for(int ii=0; ii< w; ii++)
      printf("%s%d%s", COEF(p1,p2,ii)==0 ? "":"x^",
	     COEF(p1,p2,ii)==0 ? 1:COEF(p1,p2,ii),
	     ii+1==w ? " ] ":", ");
    printf(" wght=%d\n",w);
  }

}


int var_init(int argc, char **argv, par_t *p){
  
  if (argc==1)
    ERROR("try '%s --help'\n",argv[0]);

  memset(p->r,    '\0', sizeof(p->r));
  memset(p->poly, '\0', sizeof(p->poly));
  memset(p->wght, '\0', sizeof(p->wght));
  int cnt=0;

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
  };

  for(int i=1; i<argc; i++){   /** scan for rows, cols, ell second */
    int dbg;
    if (sscanf(argv[i],"ell=%d",& dbg)==1){ /** `ell` */
      p -> ell = dbg;
      if (p->debug&4)
	printf("# read %s, setting ell=%d\n",argv[i],p->ell);
    }
    else if (sscanf(argv[i],"rows=%d",& dbg)==1){ /** `rows` */
      if(dbg>MAX_ROWS)
	ERROR("parameter rows=%d must be <= MAX_ROWS=%d",dbg,MAX_ROWS);
      p -> rows = dbg;
      if (p->debug&4)
	printf("# read %s, setting rows=%d\n",argv[i],p->rows);
    }
    else if (sscanf(argv[i],"cols=%d",& dbg)==1){ /** `cols` */
      if(dbg>MAX_COLS)
	ERROR("parameter cols=%d must be <= MAX_COLS=%d",dbg,MAX_COLS);
      p -> cols = dbg;
      if (p->debug&4)
	printf("# read %s, setting cols=%d\n",argv[i],p->cols);
    }
    else if((strcmp(argv[i],"--help")==0)
            ||(strcmp(argv[i],"-h")==0)
            ||(strcmp(argv[i],"help")==0)){
      printf( USAGE , argv[0],argv[0]);
      exit (-1);
    }
  }
  
  if((p->ell==0) || (p->rows==0) || (p->cols==0)){
    printf("must specify non-zero block rows=%d block cols=%d and"
	  " circulant size ell=%d\n",p->rows,p->cols,p->ell);
    ERROR("run '%s -h' for help",argv[0]);
  }
  for(int i=0; i< p->rows; i++)
    p->r[i]=p->ell; /** default values */
  for(int i=0; i< p->cols; i++)
    p->c[i]=p->ell; /** default values */
  
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
      if((dbg<=0) || (dbg> p->ell))
	ERROR("arg[%d]='%s' invalid row num=%d in row block %d, "
	      "must be from 1 to ell = %d",
	      i,argv[i],dbg,p1,p->ell);
      p -> r[p1] = dbg;
      if (p->debug&4)
	printf("# read %s, setting r[%d]=%d\n",argv[i],p1,p->r[p1]);
    }
    else if (sscanf(argv[i],"c%d=%d",&p1, & dbg)==2){ /** `c#=` */
      if((p1<0) || (p1>=p->cols))
	ERROR("arg[%d]='%s' : invalid col index %d, must be from 0 to cols-1 = %d\n",
	      i,argv[i],p1,p->cols -1);
      if((dbg<=0) || (dbg> p->ell))
	ERROR("arg[%d]='%s' invalid col num=%d in column block %d, "
	      "must be from 1 to ell = %d",
	      i,argv[i],dbg,p1,p->ell);
      p -> c[p1] = dbg;
      if (p->debug&4)
	printf("# read %s, setting c[%d]=%d\n",argv[i],p1,p->c[p1]);
    }    
    else if ((sscanf(argv[i],"a_%d=%d%n",&p2, & dbg, &pos)==2)||
	     (sscanf(argv[i],"a%d_%d=%d%n",&p1, &p2, &dbg, &pos)==3)){
      /* a[p1,p2]=... */
      do_scan_a(i,argv[i],p1,p2,dbg,pos,0);
    }
    else if ((sscanf(argv[i],"b_%d=%d%n",&p2, & dbg, &pos)==2)||
	     (sscanf(argv[i],"b%d_%d=%d%n",&p1, &p2, &dbg, &pos)==3)){
      /* a[p1,p2]=... */
      do_scan_a(i,argv[i],p1,p2,dbg,pos,1); /** reverse polynomial */
    }
    else{
      printf("unrecognized parameter argv[%d]=%s\n",i,argv[i]);
      ERROR("try '%s --help'\n",argv[0]);	   
    }
  }
  return 0;
}

int main(int argc, char **argv){
  var_init(argc,argv, p);

  int n = 0; /* number of columns */
  for(int i=0; i< p->cols; i++){
    int ci= p->c[i]; /* cols in this block */
    n += ci;
  }

  int rc =0; /* row count */
  int nz=0; /* non-zero elements */
  for(int i=0; i< p->rows; i++){
    int ri= p->r[i]; /* rows in this block */
    rc += ri;
    for(int j=0; j< p->cols; j++){
      nz+= WGHT(i,j) * ri; //! OK for an estimate 
    }
  }
  if(p->debug&1)
    printf("# generating %d x %d QC matrix from %d x %d blocks of size ell=%d nz=%d\n",
	   rc,n,p->rows,p->cols,p->ell,nz);
  if(p->debug&2){ /** skeleton */
    for (int br=0; br< p->rows; br++){ /** block row index  */
      printf("# ");
      for (int bc=0; bc < p->cols; bc++){ /** block columns */
	if(WGHT(br,bc))
	  printf(" M[%d,%d]",br,bc);
	else
	  printf(" .     " );
      }
      printf("  # rows=%d\n",p->r[br]);
    }
  }
  p->mat = csr_init(NULL, rc, n, nz);
  int gr=0; /** global row */
  int ge=0;     /** global element index */
  for (int br=0; br< p->rows; br++){ /** block row index  */
    for(int ii=0; ii< p->r[br]; ii++, gr++){ /** rows in block */
      int cbeg=0; /** start column of this block */
      for(int bc=0; bc < p->cols; cbeg += p->c[bc], bc++){ /** block columns */
	int wei=WGHT(br,bc);
	for(int idx=0; idx<wei; idx++){
	  int i_shift = (ii+wei-idx)%wei; /* shifted element index */
	  //	  int c_shift = (ii+p->ell+COEF(br,bc,i_shift))%(p->ell);
	  int c_shift = (ii + p->ell - COEF(br,bc,i_shift)) % (p->ell);
	  if(p->debug&16)
	    printf("br=%d gr=%d gc=%d ii=%d bc=%d idx=%d i_s=%d c_beg=%d "
		   "c_shift=%d  ge=%d COEF=%d %s\n",
		   br,gr,cbeg+c_shift,ii,bc,idx,i_shift,cbeg,c_shift,ge,COEF(br,bc,i_shift),
		   c_shift >= p->c[bc] ? "skipping":"");
	  /** write a list of pairs */
	  if(c_shift < p->c[bc]){
	    p->mat->p[ge] = gr;
	    p->mat->i[ge++] = cbeg + c_shift;
	  }
	}
      }
    }
  }
  if((nz<ge)||(gr!=rc))
    ERROR("this should not happen: nz=%d ge=%d  rc=%d gr=%d\n",nz,ge,rc,gr);
  p->mat->nz=ge; /* actual count of non-zero elements */

  //  if(p->debug&64)    csr_out(p->mat);
  csr_compress(p->mat);


  if(p->debug&8){
    if(p->mat->cols<100){
      mzd_t *dmat = mzd_from_csr(NULL,p->mat);
      mzd_print(dmat);
      mzd_free(dmat);
    }
    else
      csr_out(p->mat);    
  }

  const size_t siz=1000;
  char comment[siz+1];
  int rank = rank_csr(p->mat);
  snprintf(comment,siz,"QC block matrix [%d,%d] ell=%d rank=%d",p->rows,p->cols, p->ell, rank);
  comment[siz]='\0';/** sanity check */
  csr_mm_write(p->out,"" /** no extension */, p->mat, comment);
  if(p->mat) csr_free(p->mat);
  if(p->debug&1)
    printf("# wrote %s to %s\n",comment,p->out);
  return 0;
}

      

  
