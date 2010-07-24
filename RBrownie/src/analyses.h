#ifndef _RBrownie_ANALYSES_H
#define _RBrownie_ANALYSES_H

#include <Rcpp.h>

// Run brownie file (fnamevect[0]) and return specified objects
//
RcppExport SEXP RunTest(SEXP fnamevect, bool retReturnTrees=false, bool retContinuous=false, bool retDiscrete=false, bool retTrees=false, bool retTaxa=false, bool retTaxasets=false);

#endif