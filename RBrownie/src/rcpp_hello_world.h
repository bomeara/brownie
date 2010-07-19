#ifndef _RBrownie_RCPP_HELLO_WORLD_H
#define _RBrownie_RCPP_HELLO_WORLD_H

#include <Rcpp.h>

RcppExport SEXP rcpp_hello_world();
RcppExport SEXP readBrownie(SEXP fnamevect);
RcppExport SEXP doASR(SEXP fnamevect);

#endif
