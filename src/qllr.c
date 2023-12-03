/**
 * @file qllr.h
 *
 * @brief vecdec - a simple vectorized decoder (quantized LLR unit)
 * 
 * Many of the functions here copied from `it++` library.
 *
 * @author Leonid Pryadko (University of California, Riverside)
 *
 * Copyright (C) 2022 Leonid Pryadko
 * University of California, Riverside
 * All rights reserved.
 *
 */
#ifdef __cplusplus
extern "C"{
#endif

#include "utils.h"  
#include "qllr.h"

#ifdef USE_QLLR 

  qllr_params_t * LLR_table = NULL;

  static inline double pow2i(int x){
    if(x>=0)
      return (double) (1<<x);
    return 1.0/(1<<(-x));
  }
  
  qllr_params_t * init_LLR_tables (const int d1, const int d2, const int d3){
    if(d1<d3)
      ERROR("invalid QLLR table parameters d1=%d should not be smaller than d3=%d", d1,d3);
    if((d2<0)||(d1<0))
      ERROR("invalid QLLR table parameters d1=%d d2=%d should both be non-negative", d1,d2);
    qllr_params_t * ans = malloc(sizeof(qllr_params_t) + sizeof(qllr_t)*d2);
    ans->Dint1=d1;
    ans->Dint2=d2;
    ans->Dint3=d3;
    const double delta = pow2i(d3 - d1);
    for (int i = 0; i < d2; i++) {
      double x = delta * i;
      ans->logexp_table[i] = llr_from_dbl(log(1.0 + exp(-x)));
    }
    return ans;
  } 

  qllr_t boxplus(const qllr_t x, const qllr_t y){
    const int Dint2 = LLR_table->Dint2;

    const qllr_t a = abs(x+y);
    const qllr_t b = abs(x-y);
    const qllr_t term1 = (a-b) >>1;

    if (Dint2 == 0) {  // logmax approximation - avoid looking into empty table
      // Don't abort when overflowing, just saturate the QLLR
      if (term1 > QLLR_MAX) {
	return QLLR_MAX;
      }
      if (term1 < -QLLR_MAX) {
	return -QLLR_MAX;
      }
      return term1;
    }

    qllr_t term2 = logexp(a);
    qllr_t term3 = logexp(b);
    qllr_t result = term1 + term2 - term3;

    if (result > QLLR_MAX) {
      return QLLR_MAX;
    }
    if (result < -QLLR_MAX) {
      return -QLLR_MAX;
    }
    return result;
}

  void out_LLR_params(qllr_params_t *lcu){
    printf( "---------- LLR calculation unit -----------------\n");
    printf( "LLR_calc_unit table properties:\n");
    printf( "The granularity in the LLR representation is %g \n", pow2i(- lcu->Dint1));
    printf( "The LLR scale factor is %d\n", 1 << lcu->Dint1);
    printf( "The largest LLR that can be represented is %g\n", dbl_from_llr(QLLR_MAX));
    printf( "The table resolution is %g\n", pow2i(lcu->Dint3 - lcu->Dint1));
    printf( "The number of entries in the table is %d\n", lcu->Dint2);
    printf( "The tables truncates at the LLR value %g\n",pow2i(lcu->Dint3 - lcu->Dint1) * lcu->Dint2);
    printf( "-------------------------------------------------\n");
}

#else /* not USE_QLLR */
  void out_LLR_params([[maybe_unused]] qllr_params_t *lcu){
    printf( "---------- LLR calculation unit uses exact double arithmetic ------------\n");
  }

  qllr_params_t * init_LLR_tables (const int d1, const int d2, const int d3){
    /** all must be zero */
    if((d1!=0)||(d2!=0)||(d3!=0))
      ERROR("define 'USE_QLLR' variable to set QLLR table parameters [%d,%d,%d]",d1,d2,d3);
    qllr_params_t * ans = malloc(sizeof(qllr_params_t) + sizeof(qllr_t)*d2);
    if(!ans)
      ERROR("memory allocation");
    ans->Dint1=d1;
    ans->Dint2=d2;
    ans->Dint3=d3;
    return ans;
  } 
  
#endif /* USE_QLLR */
  
#ifdef __cplusplus
}
#endif
