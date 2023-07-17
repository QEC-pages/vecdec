/**
 *  @file nbp_mat.cpp
 *
 * @brief estimate LER on a graph using non-backtracking paths
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
#include <time.h>
#include <unistd.h>
#include <m4ri/m4ri.h>
#include "utils.h"
#include "util_m4ri.h"
#include <complex>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SparseCore>
#include <vector>
#include <set>
#include "vecdec.h"

#include <utility>

typedef Eigen::SparseMatrix<double> SpMat;
typedef std::pair<int,int> iiPair;
typedef std::vector<iiPair> iivec;
//typedef Eigen::Triplet<int> T;

/** @brief Given a DEM `graph`, create weighted Hashimoto matrix
 *
 *
 */
void do_Hashimoto(params_t p){
  int nn = 2 * p->n; /** number of arcs */
  SpMat WH(nn,nn);
  Eigen::VectorXd b(m);
  /** locate count all leaves and classify them by connectivity */
  int n1=0;
  iivec leaves;
  leaves.reserve(p->n);
  std::set<int> used; /** leaves already processed */
  int type=0; /** edge type by connectivity */
  for(int a1 = 0; a1 < p->n; a1++){
    if(p->mHt->i[a1]+1 == p->mHt->i[a1+1]){ /** one detector event in this row */
      if(used.find(a1) != used.end()){ /** new element */   
        used.insert(a1);
        leaves.pushback(pair(a1,type));
        for(auto it = back(vector); it != end (vector); ++it) {
          if (used.find(*it) != used.end()){
            if(p->mHt->i[(it]+1 == p->mHt->i[*it+1]){ /** one detector event in this row */
            
          it->doSomething ();
        for(iivec. vector<int>::iterator iter 
        lea
        
      }
      leaves.push_back(a1);
  }
  std::sort(leaves.begin(), leaves.end());
  
}
