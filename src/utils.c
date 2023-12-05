
/**
 * @file utils.c
 *
 * @brief Collection of utility functions
 *
 * @author Leonid Pryadko (University of California, Riverside)
 *
 * Copyright (C) 2022 Leonid Pryadko
 * University of California, Riverside
 * All rights reserved.
 *
 */
#include <strings.h>
#include <ctype.h>
#include <m4ri/m4ri.h>
#include "util_m4ri.h"
#include "utils.h"
#include "mmio.h"

tinymt64_t tinymt;

/** @brief read an `MMX` array of doubles
 * With a column of doubles, `nrows` and `ncols` point to zeros on return.
 * @param[in] siz number of entries in the pre-allocated array 
 * @param[out] siz actual number of entries read
 * @param *arr must be allocated to size `siz` or NULL if `siz==0`.
 * @return the pointer to the array with the data
 */
double * dbl_mm_read(const char * const fin, int *nrows, int *ncols, int *siz, double *  arr){
  MM_typecode matcode;
  FILE *f;
  int num_items;
  *nrows = 0;
  *ncols = 0;

  if ((f = fopen(fin, "r")) == NULL) 
    ERROR("can't open file %s",fin);

  if (mm_read_banner(f, &matcode) != 0)
    ERROR("Could not process Matrix Market banner.");

  if (!(mm_is_matrix(matcode) && mm_is_dense(matcode) && 
	mm_is_real(matcode) && mm_is_general(matcode) )){
    printf("Sorry, this application does not support ");
    printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
    ERROR("input file %s",fin);
    exit(1);
  }
  
  if(mm_read_mtx_array_size(f, nrows, ncols))
    ERROR("Can't read array size");
  if(arr != NULL){/** allocated array */
    if(*siz <= 0)
      ERROR("invalid input value siz=%d",*siz);
    if((*siz) < (*nrows) * (*ncols)){
      *siz = (*nrows) * (*ncols);
      arr = realloc(arr,sizeof(double)*(*siz));
      if(!arr)
	ERROR("memory allocation failed");
    }
  }
  else{
    if(*siz != 0)
      ERROR("invalid input value siz=%d, arr==NULL",*siz);
    *siz = (*nrows) * (*ncols);
    arr = malloc(sizeof(double)*(*siz));
    if(!arr)
      ERROR("memory allocation failed");
  }
  num_items = *siz;
  double *ptr = arr; 
  for(int i=0; i < num_items; i++, ptr++){
    if(1 != fscanf(f," %lg ", ptr))
      ERROR("failed to read entry %d of %d",i,num_items);
  }
  return arr;
}

void dbl_mm_write( char * const fout, const char fext[],
		   const int rows, const int cols, const double buf[],
		  const char comment[]){
  int result=0; /**< non-zero if write error */
  size_t len=strlen(fout)+strlen(fext)+1;
  char *str;
  
  FILE *f;
  if(strcmp(fout,"stdout")!=0){
    str=calloc(len, sizeof(char));
    sprintf(str,"%s%s",fout,fext);
    f=fopen(str,"w");
  }
  else{/** fout == "stdout" */
    str=fout;
    f=stdout;
  }
  
  if(!f)
    ERROR("can't open file '%s' for writing",str);
  if(fprintf(f,"%%MatrixMarket matrix array real general\n")<0)
    result++;
  if(comment!=NULL){
    if(fprintf(f,"%% %s\n",comment)<0)
      result++;
  }
  if(fprintf(f,"%d %d\n",rows,cols)<3)
    result++;
  const size_t num=rows*cols;
  for(size_t i=0;i<num;i++)
    if(fprintf(f,"%+18.12g\n", buf[i])<0)
      result++;
  if(result)
    ERROR("error writing to file '%s'",str);
  
  if(strcmp(fout,"stdout")!=0){
    fclose(f);
    free(str);
  }
}

/** @brief read detector error model (DEM) created by `stim`.
 * Immediately create CSR matrices `mH` and `mL` and vector `vP`; 
 * return the corresponding pointers via `ptrs` (in this order).
 * TODO: make it test for file type 
 * @param fnam file name for reading DEM from
 * @param p structure to store produced matrices
 */
void read_dem_file(char *fnam, void * ptrs[3], int debug){
  ssize_t linelen, col=0;
  size_t lineno=0, bufsiz=0; /**< buffer size for `readline` */
  char *buf = NULL;          /** actual buffer for `readline` */
  //  int nzH=0, nzL=0;  /** count non-zero entries in `H` and `L` */
  int maxH=100, maxL=100, maxN=100; 
  double *inP = malloc(maxN*sizeof(double));
  int_pair * inH = malloc(maxH*sizeof(int_pair));
  int_pair * inL = malloc(maxL*sizeof(int_pair));
  if ((!inP)||(!inH)||(!inL))
    ERROR("memory allocation failed\n");

  if(debug & 1)
    printf("# opening DEM file %s\n",fnam);
  FILE *f = fopen(fnam, "r");
  if(f==NULL)
    ERROR("can't open the (DEM) file %s for reading\n",fnam);

  int r=-1, k=-1, n=0;
  int iD=0, iL=0; /** numbers of `D` and `L` entries */
  do{ /** read lines one-by-one until end of file is found *************/
    lineno++; col=0; linelen = getline(&buf, &bufsiz, f);
    if(linelen<0)
      break;
    if(debug & 32) printf("# %s",buf);
    char *c=buf;
    double prob;
    int num=0, val;
    while(isspace(*c)){ c++; col++; } /** `skip` white space */
    if((*c != '\0')&& (*c != '#') &&(col < linelen)){
      if(sscanf(c,"error( %lg ) %n",&prob,&num)){
        if((prob<=0)||(prob>=1))
          ERROR("probability should be in (0,1) exclusive p=%g\n"
                "%s:%zu:%zu: '%s'\n", prob,fnam,lineno,col+1,buf);
        c+=num; col+=num;
        if(n>=maxN){
          maxN=2*maxN;
          inP=realloc(inP,maxN*sizeof(*inP));
        }
        inP[n]=prob;
        do{/** deal with the rest of the line */
          num=0;
          if(sscanf(c," D%d %n",&val, &num)){/** `D` entry */
            c+=num; col+=num;
            assert(val>=0);
            if(val>=r)
              r=val+1;  /** update the number of `D` pairs */
            if(iD>=maxH){
              maxH=2*maxH;
              inH=realloc(inH,maxH*sizeof(*inH));
            }
            inH[iD].a   = val;   /** add a pair */
            inH[iD++].b = n;
            if(debug & 32) printf("n=%d iD=%d val=%d r=%d\n",n,iD,val, r);
          }
          else if(sscanf(c," L%d %n",&val, &num)){/** `L` entry */
            c+=num; col+=num;
            assert(val>=0);
            if(val>=k)
              k=val+1;  /** update the number of `L` pairs */
            if(iL>=maxL){
              maxL=2*maxL;
              inL=realloc(inL,maxL*sizeof(*inL));
            }
            inL[iL].a   = val;   /** add a pair */
            inL[iL++].b = n;
            if(debug & 32) printf("n=%d iL=%d val=%d k=%d\n",n,iD,val,k);
          }
          else
            ERROR("unrecognized entry %s"
		  "%s:%zu:%zu: '%s'\n",c,fnam,lineno,col+1,buf);
        }
        while((c[0]!='#')&&(c[0]!='\n')&&(c[0]!='\0')&&(col<linelen));
        n++;
      }
      else if (sscanf(c,"detector( %d %n",&val,&num)){
        /** do nothing */
        //        printf("# ignoring row[%zu]=%s\n",lineno,c);
      }
      else if (sscanf(c,"shift_detectors( %d %n",&val,&num)){
        /** do nothing */
        //        printf("# ignoring row[%zu]=%s\n",lineno,c);
      }
      else
        ERROR("unrecognized DEM entry %s"
              "%s:%zu:%zu: '%s'\n",c,fnam,lineno,col+1,buf);

    }
    /** otherwise just go to next row */
  }
  while(!feof(f));

  if(debug &1)
    printf("# read DEM: r=%d k=%d n=%d\n",r,k,n);

  ptrs[0] = csr_from_pairs(ptrs[0], iD, inH, r, n);
  ptrs[1] = csr_from_pairs(ptrs[1], iL, inL, k, n);
  ptrs[2] = inP;
  if (buf)
    free(buf);
  free(inH);
  free(inL);
  fclose(f);

}


#ifdef __MINGW32__ /** windows compiler */
int getline(char **line, size_t *n, FILE *fp){
  size_t size, oldsize;
  size_t len=0;
  if ((line==NULL) || (n==NULL) || (fp== NULL))
    ERROR("Null input parameters!");
  if (ferror (fp))
    return -1;
  if (feof(fp))
    return -1;
  for(oldsize=0, size=BUFSIZ ; ; size += BUFSIZ){
    /* BUFSIZ is "the optimal read size for this platform" */
    char *block;
      if (*n<=size){
	block = (char *) realloc(*line,size+1);
	/* realloc(NULL,n) equiv malloc(n) */
	if (block !=NULL){  /* reallocation worked */
	  *n=size+1;
	  *line=block;
	}
	else	   /* the original pointer may have been retained */
	  ERROR("memory allocation in getline");	
      }
      block = *line + oldsize; 
      if(NULL==fgets(block,BUFSIZ+1,fp)){
	if (oldsize>0)
	  break;  /* there is a string present */
	else
	  return -1;  /* nothing can be read - read error or EOF*/
      }
      len = oldsize+strlen(block);      
      if ((*line)[len-1]=='\n'){
	(*line)[--len]='\0';	
	break;  /* we are done */
      }
      if (feof(fp))
	break;
      oldsize=size; /* fgets puts a terminal '\0' on the end of the
		       string, so we make sure to overwrite this */
  } 
  return len;   
}

#endif /* __MINGW32__ */

/** verify that line has only space */
int all_space(const char * str) {
  while (*str) 
    if (!isspace(*str++)) 
      return 0;    
  return 1;
}

/** @brief read up to `lmax` lines from a file in `01` format

 * read up to `lmax` binary vectors of length `m` from a `01` file `fin` open
 * for reading.  Place the vectors as columns of matrix `M` of size `m` rows by
 * `lmax` colums.  Lines starting with `#` are silently ignored; a non-`01`
 * line, or a `01` line of an incorrect length will give an error.
 *
 * @param M initialized output matrix with `lmax` rows and `m` columns
 * @param fin file with 01 data open for reading
 * @param[input,output] lineno current line number in the file.
 * @param fnam file name (for debugging purposes)
 * @param p Other parameters (only `p->debug` is used).
 * @return the number of rows actually read.
 *
 */
int read_01(mzd_t *M, FILE *fin, int *lineno, const char* fnam,
	      const int debug){
  if(!M)
    ERROR("expected initialized matrix 'M'!\n");
  else
    mzd_set_ui(M,0);
  int m   =M->nrows;
  int lmax=M->ncols, il=0;
  if(!fin)
    ERROR("file 'fin' named '%s' must be open for reading\n",fnam);
  if(debug&8) /** file io */
    printf("# about to read 01 data from line %d in file '%s'\n",
           *lineno,fnam);

  char *buf=NULL;
  size_t bufsiz=0;

  ssize_t linelen;
  while((il<lmax) && (!feof(fin)) &&
        ((linelen = getline(&buf, &bufsiz, fin))>=0)){
    (*lineno)++;
    switch(buf[0]){
    case '0': case '1':
      if(linelen<=m)
	ERROR("line is too short, expected %d 01 characters\n"
	      "%s:%d:1: '%s'\n", m,fnam,*lineno,buf);
      else{
	for(int i=0; i<m; i++){
	  if (buf[i]=='1')
	    mzd_write_bit(M,i,il,1); /** row `i`, col `il` */
	  else if (buf[i]!='0')
	    ERROR("invalid 01 line\n"
		  "%s:%d:%d: '%s'\n", fnam,*lineno,i+1,buf);
	}
	(il)++; /** success */
      }
      break;
    case '#':       /** do nothing - skip this line */
      break;
    default:
      if (!all_space(buf))
	ERROR("invalid 01 line\n"
	      "%s:%d:1: '%s'\n", fnam,*lineno,buf);
      break;
    }
  }
  if(debug&8) /** file io */
    printf("# read %d 01 rows from file '%s'\n",il,fnam);
  if(buf)
    free(buf);
  return il;
}

