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
#define RS(i) (p->rs[((MAX_W)+1)*i])
#define CS(i) (p->cs[((MAX_W)+1)*i])
#define BLOCK_INDEX(i,j) (p->cols*(i)+(j)) //! row=i, col=j 
#define WGHT(i,j) p->poly[((MAX_W)+1)*BLOCK_INDEX(i,j)]
//! lists of non-zero coefficients indices `x+nx*y` 
#define INDEX(i,j,k) p->poly[((MAX_W)+1)*BLOCK_INDEX(i,j)+(k)+1]
//! index=`x+nx*y` of term `(x,y)` -- values out of range to be continued periodically 
#define IDX(x,y) (p->nx * (((y) + 100 * p->ny) % p->ny) + (((x) + 100 * p->nx) % p->nx))
//! x and y positions corresponding to the `index` value
#define DEG_X(index) ((index) % p->nx)
#define DEG_Y(index) ((index) / p->nx)

typedef struct PAR_T {
  int rows; /** how many circulant rows, default 1 */
  int cols; /** how many circulant cols, default 1 */
  int ell; /** =nx*ny, group order */
  int nx, ny; /** degrees of generators x and y, respectively; x^nx=y^ny=1 */
  /** nx must be bigger than 1 **/
  /** a polynomial is stored as weight `w` (in position 0), followed by `w` non-negative indices */
  int rs[(MAX_ROWS)*(MAX_W+1)]; /** rows to skip in each row block (x and y degrees) */
  int cs[(MAX_COLS)*(MAX_W+1)]; /** cols to skip in each col block (x and y degrees) */
  int r[MAX_ROWS]; /* how many rows in each row block - computed from `rs` */
  int c[MAX_ROWS]; /* how many cols in each col block - computed from `cs` */
  int poly[(MAX_ROWS)*(MAX_COLS)*(MAX_W+1)]; /** coefficients of polynomials */
  char *out; /** name of the output file */
  int tot_cols,tot_rows; /** actual matrix size: `n=c[0]+c[1]+... +c[cols]`, `r=r[0]+r[1]+... r[rows]` */
  int tot_nz; /** total number of non-zero elements */
  csr_t *mat; /** constructed matrix */
  int debug;
} par_t;

par_t params={
  .rows=1,
  .cols=0,  //! must specify
  .nx=0,    //! must specify
  .ny=0,    //! must specify
  .ell=0,
  .out="std" "out",
  .mat=NULL,
  .debug=1
};

par_t *p = &params;

#define USAGE								\
  "%s: mtx_qc - generate a quasi-abelian two-generator (bivariate-bicycle) matrix\n" \
  "  usage: %s param=value [[param=value] ... ]\n"			\
  "\t Parameters are read in the order specified.\n"			\
  "\t Supported parameters:\n"						\
  "\t --help\t\t: give this help (also '-h' or just 'help')\n"		\
  "\t out=[string]\t: name of the output file (or 'stdout', default value)\n" \
  "\t rows=int\t: number of row blocks (default: 1)\n"			\
  "\t cols=int\t: number of column blocks (no default)\n"		\
  "\t nx=int\t: degree of x-generator (no default)\n"			\
  "\t ny=int\t: degree of y-generator (no default)\n"			\
  "\t a_#=int_int[,int_int,...]\t 'x_y' degrees of row=0,col=# polynomial (same as 'a0_#')\n" \
  "\t a#_#=int_int[,int_int,...]\t 'x_y' degrees of polynomial in (row,col)\n" \
  "\t b_#=... b#_#=...\t same but also invert each monomial (flip degree signs)\n" \
  "\t\t It is ok to skip zero values: '0', '0_', '0_0', and '_0' are equivalent,\n" \
  "\t\t '1', '1_', and '1_0' give (1,0), '_3' and '0_3' both give (0,3)\n" \
  "\t rs#=int_int[,int_int,...]\t: 'x_y' rows to skip in row block '#' \n" \
  "\t cs#=int_int[,int_int,...]\t: 'x_y' cols to skip in col block '#' \n" \
  "\t\t Specify only non-zero polynomials\n"				\
  "\t\t Polynomial x,y deegrees will be adjusted mod nx, mod ny\n"	\
  "\t debug=[integer]\t: bitmap for aux information to output (default: 1)\n" \
  "\t .1 (bit 0) output helpful information, e.g. ranks and code params\n" \
  "\t .2 (bit 1) output skeleton map of the code\n"                     \
  "\t .4 (bit 2) echo input parameters being read\n"			\
  "\t .8 (bit 3) print out the constructed matrix\n"			\
  "\t .16 (bit 4) extra processing information\n"			\
  "\t Multiple 'debug' parameters are XOR combined except for 0.\n"	\
  "\t Space is allowed between 'out=' and the file name string.\n"


void print_as_poly(const int *const ptr, const int p1, const int p2, const int reverse, const int init){
  int as_poly = 0;
  if ((p1>=0)&&(p2>=0)){
    printf("%s%c(%d,%d) = ",init ? "# read polynomial ":"",reverse? 'b' : 'a',p1,p2);
      as_poly =1;
  }
    else if ((p1<0)&&(p2>=0))
      printf("%scs(%d) = ",init? "# skip list ":"",p2);
    else if ((p1>=0)&&(p2<0))
      printf("%srs(%d) = ",init ? "# skip list ":"",p1);
    else
      ERROR("invalid indices p1=%d p2=%d",p1,p2);
    const int w = ptr[0];
    for(int ii=1; ii<= w; ii++){
      int x=DEG_X(ptr[ii]);
      int y=DEG_Y(ptr[ii]);
      if(x==0){
	if(y==0)
	  printf("1");
	else if(y==1)
	  printf("y");
	else
	  printf("y^%d",y);
      }
      else if (x==1){
	if(y==0)
	  printf("x");
	else if(y==1)
	  printf("x*y");
	else
	  printf("x*y^%d",y);
      }
      else{
	if(y==0)
	  printf("x^%d",x);
	else if(y==1)
	  printf("x^%d*y",x);
	else
	  printf("x^%d*y^%d",x,y);
      }
      printf("%s", 	     ii==w ? "; " : as_poly ? " + " : " , " );
    }
    if (as_poly)
      printf(" wght=%d\n",w);
    else
      printf(" total %d %s removed, block %s remain %d\n",w,p1<0?"cols":"rows",p1<0?"cols":"rows",p->ell-w);
      

}

/** scan a list of coefficients from a string
 * format:= '(pair0)[,(pair1),...]', (pair) := 'x_y' or just 'x' for (x,0) 
 * @return the number of terms scanned 
 * @param[out] ptr place to put the coefficient indices 
 * @param[in] i position of the command line argument in `argv[]` (for debugging)
 * @param[in] argvi string to process (e.g., `argv[i]`)
 * @param pos current position in the string `argvi` 
 * @param p1 block `row` index 
 * @param p2 block `col` index
 * @param reverse invert each monomial if non-zero
 * @return the number of coefficients read 
 * 
 */
int do_scan_a(int * const ptr, int i, const char argvi[], int glbpos, int p1, int p2, int reverse){

  int w=0, dbgx, dbgy, num, pos1, pos;
  if((p1>=p->rows)||(p2>=p->cols))
    ERROR("arg[%d]='%s' : %c(%d,%d)=... invalid block position rows=%d cols=%d",
	  i,argvi,reverse?'b':'a',p1,p2,p->rows,p->cols);
  const char *c = argvi + glbpos;
  do{
    pos1=pos=-1;
    dbgx=-1;
    dbgy=-1;
    num=sscanf(c,"%d%n_%d%n",&dbgx,&pos1,&dbgy,&pos);
    switch(num){
    case 1:
      dbgy=0;
      pos=pos1;
      /* fall through */
    case 2: case 3:
      /*      printf("recording dbgx=%d dbgy=%d w=%d idx=%d\n",
	      dbgx,dbgy,w,reverse==0 ? IDX(dbgx,dbgy) : IDX(-dbgx,-dbgy)); */
      /** this is correct: reserve 0th position for the weight */
      ptr[++w] = reverse==0 ? IDX(dbgx,dbgy) : IDX(-dbgx,-dbgy);
      if(c[pos]==','){ /** we are looking at `#,` or `#_#,` format */
	//	printf("skip comma\n");
	pos++;
      }
      glbpos += pos;
      c += pos;
      break;
    case 0:
      ERROR("arg[%d]='%s' at position %d w=%d : trailing garbage %s",i,argvi,glbpos,w,c);
      break;
    default:
      ERROR("arg[%d]='%s' at position %d w=%d : invalid entry!",i,argvi,glbpos,w);
      break;
    }
    if(w>MAX_W)
      ERROR("arg[%d]='%s' pos=%d: too many terms: weight>=%d MAX_W=%d\n", i,argvi,pos,w,MAX_W);
  } while(c[0]);
  ptr[0]=w;
  qsort(ptr+1,w,sizeof(int),cmp_rci_t);
  for(int ii=1; ii< w; ii++){
    if(ptr[ii]==ptr[ii+1])
      ERROR("arg[%d]='%s' w=%d : two identical terms x^%d*y^%d ",
	    i,argvi,w,DEG_X(ptr[ii]),DEG_Y(ptr[ii]));
  }
  if(p->debug&4){
    print_as_poly(ptr,p1,p2,reverse,1);
  }
  return w;
}


int var_init(int argc, char **argv, par_t *p){
  
  if (argc==1)
    ERROR("try '%s --help'\n",argv[0]);

  memset(p->rs,    '\0', sizeof(p->rs));
  memset(p->cs,    '\0', sizeof(p->cs));
  memset(p->r,    '\0', sizeof(p->r));
  memset(p->c,    '\0', sizeof(p->c));
  memset(p->poly, '\0', sizeof(p->poly));

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
    if (sscanf(argv[i],"nx=%d",& dbg)==1){ /** `nx` */
      p -> nx = dbg;
      if (p->debug&4)
	printf("# read %s, setting nx=%d\n",argv[i],p->nx);
    }
    if (sscanf(argv[i],"ny=%d",& dbg)==1){ /** `ny` */
      p -> ny = dbg;
      if (p->debug&4)
	printf("# read %s, setting ny=%d\n",argv[i],p->ny);
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
  
  if((p->nx==0) || (p->ny==0) || (p->rows==0) || (p->cols==0)){
    printf("must specify non-zero block rows=%d, block cols=%d, and"
	   " generator orders nx=%d, ny=%d\n",p->rows,p->cols,p->nx, p->ny);
    ERROR("run '%s -h' for help" ,argv[0]);
  }
  p->ell = p->nx * p->ny;

  for(int i=1; i<argc; i++){   /** remaining arguments */
    int pos=0, dbg, p1=0, p2=0;
    if((sscanf(argv[i],"debug=%d",& dbg)==1) ||
       (sscanf(argv[i],"rows=%d",& dbg)==1)  ||
       (sscanf(argv[i],"cols=%d",& dbg)==1)  ||
       (sscanf(argv[i],"nx=%d",& dbg)==1)  ||
       (sscanf(argv[i],"ny=%d",& dbg)==1)){
      /** do nothing, already processed all `debug` entries */
    }
    else if (0==strncmp(argv[i],"out=",4)){ /** `out` */
      if(strlen(argv[i])>4){
	p->out = argv[i]+4;
	if (p->debug&4)
	  printf("# read %s, out=%s\n",argv[i],p->out);
      }
      else{
	p->out = argv[++i]; /**< allow space before file name */
	if (p->debug&4)
	  printf("# read %s=%s, out=%s\n",argv[i-1],argv[i],p->out);
      }
    }
    else if ((sscanf(argv[i],"a_%d=%n",&p2,&pos)==1)||
	     (sscanf(argv[i],"a%d_%d=%n",&p1, &p2, &pos)==2)){
      /* a[p1,p2]=... */
      do_scan_a(&WGHT(p1,p2),i,argv[i],pos,p1,p2,0);
    }
    else if ((sscanf(argv[i],"b_%d=%n",&p2,&pos)==1)||
	     (sscanf(argv[i],"b%d_%d=%n",&p1, &p2, &pos)==2)){
      /* b[p1,p2]=... */
      do_scan_a(&WGHT(p1,p2),i,argv[i],pos,p1,p2,1); /** transposed version */
    }
    else if (sscanf(argv[i],"rs%d=%n",&p1,&pos)==1){ /** `rs#=` */
      if((p1<0) || (p1>=p->rows))
	ERROR("arg[%d]='%s' : invalid row index %d, must be from 0 to rows-1 = %d\n",
	      i,argv[i],p1,p->rows -1);
      int w = do_scan_a(&RS(p1),i,argv[i],pos,p1,-1,0);
      if (p->debug&4)
	printf("# read %s, setting %d skip row%s in block row %d\n",argv[i],w,w>1?"s":"",p1);
    }
    else if (sscanf(argv[i],"cs%d=%n",&p1,&pos)==1){ /** `cs#=` */
      if((p1<0) || (p1>=p->cols))
	ERROR("arg[%d]='%s' : invalid col index %d, must be from 0 to cols-1 = %d\n",
	      i,argv[i],p1,p->cols -1);
      int w = do_scan_a(&CS(p1),i,argv[i],pos,-1,p1,0);
      if (p->debug&4)
	printf("# read %s, setting %d skip cols in block column %d\n",argv[i],w,p1);
    }    
    else{
      printf("unrecognized parameter argv[%d]=%s\n",i,argv[i]);
      ERROR("try '%s --help'\n",argv[0]);	   
    }
  }

  // set block row and column dimensionso
  for(int i=0; i < p->rows; i++){
    p->r[i] = p->ell - RS(i);
    if (p->r[i] <=0){
      print_as_poly(&RS(i),i,-1,0,0);
	ERROR("too many rows removed in row block %d, ell=%d removed %d",
	      i,p->ell,RS(i));
      }
    if((p->debug&2)&&RS(i))
      printf("# block row i=%d removed %d rows, remain r=%d\n",i,RS(i),p->r[i]);
  }
  
  for(int i=0; i < p->cols; i++){
    p->c[i] = p->ell - CS(i);
    if (p->c[i] <=0){
      print_as_poly(&CS(i),-1,i,0,0);
      ERROR("too many cols removed in col block %d, ell=%d removed %d",
	    i,p->ell,CS(i));
    }
    if((p->debug&2)&&CS(i))
      printf("# block col i=%d removed %d cols, remain c=%d\n",i,CS(i),p->c[i]);
  }

  int n = 0; /* number of columns */
  for(int i=0; i< p->cols; i++)
    n += p->c[i]; /* cols in this block */  
  p->tot_cols = n;

  int rc =0; /* row count */
  int nz=0; /* non-zero elements */
  for(int i=0; i< p->rows; i++){
    int ri= p->r[i]; /* rows in this block */
    rc += ri;
    for(int j=0; j< p->cols; j++){
      nz+= WGHT(i,j) * ri; //! OK for an estimate 
    }
  }
  p->tot_rows = rc;
  p->tot_nz = nz;
  
  return 0;
}

int main(int argc, char **argv){
  var_init(argc,argv, p);

  if(p->debug&1)
    printf("# generating %d x %d %s matrix from %d x %d blocks of size nx=%d ny=%d ell=%d nz <= %d\n",
	   p->tot_rows,p->tot_cols,p->ny>1 ? "BB" : "QC", p->rows,p->cols,p->nx,p->ny,p->ell,p->tot_nz);

  if(p->debug&2){ /** skeleton */
    for (int bc=0; bc<p->cols;bc++)
      printf("%s| %4d %s",bc==0?"# ":"", p->c[bc], bc+1 == p->cols ? "| # cols\n":"");
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
    for(int br=0; br< p->rows; br++){ /** block row index  */
      for (int bc=0; bc < p->cols; bc++){ /** block columns */
	if(WGHT(br,bc)){
	  printf(" M[%d,%d]: ",br,bc);
	  print_as_poly(&WGHT(br,bc),br,bc,0,0);
	}
      }
    }
  }
  
  /** actual matrix construction ***************************** */
  p->mat = csr_init(NULL, p->tot_rows, p->tot_cols, p->tot_nz);
  int true_row=0; /** global row */
  int ge=0;     /** global element index */
  for (int br=0; br< p->rows; br++){ /** block row index  */
    for(int y0=0; y0<p->ny; y0++){
      for (int x0=0; x0<p->nx; x0++){
	int idx0 = IDX(x0,y0);
	//	printf("br=%d x=%d y=%d idx=%d true_row=%d\n",br,x0,y0,idx0,true_row);
	if(RS(br) && bsearch(&idx0, (&RS(br))+1, RS(br), sizeof(int), cmp_rci_t)){
	  if(p->debug&16)
	    printf("# br=%d : skipping row idx=%d x=%d y=%d\n",br,idx0,x0,y0);
	}
	else{ /** this is global row `true_row` with local shift by (x0,y0) **************/
	  int cbeg=0; /** starting column of a current  block */
	  for(int bc=0; bc < p->cols; cbeg += p->c[bc], bc++){ /** block columns */
	    int wei=WGHT(br,bc); /** weight of this block polynomial */
	    for(int ii=0; ii<wei; ii++){
	      const int iidx = INDEX(br,bc,ii);
	      const int idx = IDX(x0 + DEG_X(iidx),y0 + DEG_Y(iidx)); /** index after (x0,y0) shift */
	      /** calculate the true column in this block with any */
	      int true_col = cbeg + idx;
	      for(int i=1; i<=CS(bc); i++){
		int * const skip = &CS(bc);
		if (skip[i] < idx)
		  true_col --;
		else if (skip[i] == idx){ /** we are skipping this column */
		  if(p->debug&16)
		    printf("# skipping br=%d bc=%d x=%d y=%d idx=%d true_col=%d\n",br,bc,DEG_X(idx),DEG_Y(idx),idx,true_col);
		  true_col = -1;
		  break;
		}
	      }
	      if(true_col >=0){ /** add to the list of pairs */
		p->mat->p[ge] = true_row;
		p->mat->i[ge++] = true_col;
	      }
	    }
	  }
	  true_row++;
	}
      }
    }
  }
  if((p->tot_nz < ge)||(true_row != p->tot_rows))
    ERROR("this should not happen: nz=%d ge=%d  rc=%d gr=%d\n",p->tot_nz,ge,p->tot_rows,true_row);
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
  snprintf(comment,siz,"ell=%d (%d x %d) %s [%d x %d] block matrix dims= %d x %d rank=%d",
	   p->ell, p->nx,p->ny,
	   p->ny>1 ? "BB" : "QC",
	   p->rows,p->cols, p->tot_rows,p->tot_cols,rank);
  comment[siz]='\0';/** sanity check */
  csr_mm_write(p->out,"" /** no extension */, p->mat, comment);
  if(p->mat) csr_free(p->mat);
  if(p->debug&1)
    printf("# wrote %s to %s\n",comment,p->out);
  return 0;
}

      

  
