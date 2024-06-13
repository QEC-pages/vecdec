#ifndef QLLR_H
#define QLLR_H
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#ifdef USE_QLLR  

  typedef  signed int qllr_t;
  static const qllr_t QLLR_MAX = ((INT_MAX) >> 4);
  static const double MINPROB = 1.0e-15;

  typedef struct QLLR_PARAMS_T {
    short int Dint1, Dint2, Dint3;
    //! `1<<Dint1` determines how integral LLRs relate to real LLRs (to_double=(1<<Dint1)*int_llr)
    //! `Dint2` is number of entries in table for LLR operations
    //! table resolution is `2^(-(Dint1-Dint3))`
    int logexp_table[0]; /** variable size array of length `Dint2` */
  } qllr_params_t;

  extern qllr_params_t *LLR_table;
  
  qllr_params_t * init_LLR_tables (const int d1, const int d2, const int d3);

  
  static inline double dbl_from_llr(const qllr_t val){
    const int Dint1 = LLR_table->Dint1;
    return ((double) val) / (1 << Dint1);
  }

  static inline qllr_t llr_from_dbl(const double val){
    const int Dint1 = LLR_table->Dint1;
    const double QLLR_MAX_double = dbl_from_llr(QLLR_MAX);
   // Don't abort when overflow occurs, just saturate the QLLR
   if (val > QLLR_MAX_double) {
     //     it_info_debug("LLR_calc_unit::to_qllr(): LLR overflow");
     return QLLR_MAX;
   }
   if (val < -QLLR_MAX_double) {
     //     it_info_debug("LLR_calc_unit::to_qllr(): LLR overflow");
     return -QLLR_MAX;
   }
   return (qllr_t ) (floor(0.5 + (1 << Dint1) * val));
 }
 

  static inline qllr_t logexp(const qllr_t x){
    const int Dint2 = LLR_table->Dint2;
    const int Dint3 = LLR_table->Dint3;    

    assert(x>=0);
    
    int ind = x >> Dint3;
    if (ind >= Dint2) // outside table
      return 0;
    assert(ind>=0);
  
    // With interpolation
    // int delta=x-(ind<<Dint3);
    // return ((delta*logexp_table(ind+1) + ((1<<Dint3)-delta)*logexp_table(ind)) >> Dint3);
  
    // Without interpolation
    return LLR_table->logexp_table[ind];
  }
   
  static inline qllr_t jaclog(const qllr_t a, const qllr_t b){
    qllr_t x, maxab;
  
    if (a > b) {
      maxab = a;
      x = a - b;
    }
    else {
      maxab = b;
      x = b - a;
    }
  
    if (maxab >= QLLR_MAX)
      return QLLR_MAX;
    else
      return (maxab + logexp(x));
  }
  
  
  qllr_t boxplus(const qllr_t x, const qllr_t y);

  static inline qllr_t llr_from_P(const double P){
    assert(P<1);
    double val = P > MINPROB ? log((1.0/P -1.0)) : log(1.0/MINPROB - 1.0);
    return llr_from_dbl(val);
  }

  
#else /* not USE_QLLR */

  static const double MINPROB = 1.0e-15;

  typedef  double qllr_t;
  
/** @brief Hagenauer boxplus operator 
 * 
 * WARNING: this requires full `LLR` values `x` and `y`.  Resulting
 * LLR is calculated using the equivalent form (`**verified**`):
 * `arctahn(tanh(x)*tanh(y))=log[cosh[(x+y)/2]/cosh[(x-y)/2]]`, which
 * is further simplified using the fact that `cosh` is an even
 * function, by introducing `a=abs(x+y); b=abs(x-y);` and
 * `logexp(x)=log(1-exp(-x))`, to give, finally
 * `0.5*(a-b)+logexp(a)-logexp(b);`

 *  @return LLR of the probability to have only one non-zero
 */

static inline qllr_t boxplus(const qllr_t x, const qllr_t y){
  const double a = fabs(x+y);
  const double b = fabs(x-y);
  return 0.5*(a-b)+log((1+exp(-a))/(1+exp(-b)));
  //  return 0.5*(a-b); /** min-sum version */
}

  static inline qllr_t llr_from_P(const double P){
    return P > MINPROB ? log((1.0/P -1.0)) : log(1.0/MINPROB - 1.0);
  }

  static inline double dbl_from_llr(const qllr_t val){
    return val;
  }

  static inline qllr_t llr_from_dbl(const double val){
    return val;
  }

  typedef struct QLLR_PARAMS_T {
    short int Dint1, Dint2, Dint3; /** placeholder */
  } qllr_params_t;

  extern qllr_params_t *LLR_table;
  
  
    
#endif /* USE_QLLR */

  qllr_params_t * init_LLR_tables (const int d1, const int d2, const int d3);
  void out_LLR_params(_maybe_unused qllr_params_t *lcu);

  /** @brief calculate error probability from the LLR value */
  static inline double P_from_llr(const qllr_t llr){
    const double x = dbl_from_llr(llr);
    if (x>5){
      double emx=exp(-x);
      return emx/(1.0+emx);
    }
    else 
      return 1.0/(1.0+exp(x));
  }

 qllr_t mzd_row_energ(qllr_t *coeff, const mzd_t *A, const int i);
  
#ifdef __cplusplus
}
#endif


#endif /* QLLR_H */
