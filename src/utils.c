
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
    arr = realloc(arr,sizeof(double)*(*siz));
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
