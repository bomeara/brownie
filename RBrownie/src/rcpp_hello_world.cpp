#include "rcpp_hello_world.h"
#include <strstream>
#include <fstream>
#include <iomanip>
#include <unistd.h>
#include <stdio.h>
#include <set>
#include <climits>
#include <cstring>
#include <memory>

#include "nexusdefs.h"
#include "xnexus.h"
#include "nexustoken.h"
#include "nexus.h"
#include "taxablock.h"
#include "assumptionsblock.h"
#include "treesblock.h"
#include "discretedatum.h"
#include "discretematrix.h"
#include "charactersblock.h"
#include "charactersblock2.h"
#include "gport.h"
#include "profile.h"
#include "nodeiterator.h"
#include "setreader.h"
#include "treeorder.h"
#include "treedrawer.h"
#include "ntree.h"
#include "stree.h"
#include "containingtree.h"
#include "quartet.h"
#include "version.h"
#include <gsl/gsl_sf_gamma.h>
#include "TreeLib.h"
#include "gtree.h"
#include "treereader.h"
#include "treewriter.h"
#include <time.h>
#include <map>
#include <limits>
#include "brownie.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_version.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_errno.h>
#include "optimizationfn.h"
#include "cdfvectorholder.h"
#include <sstream>
#include <iostream>



SEXP rcpp_hello_world()
{
	
	printf("Loading namespace");
	using namespace Rcpp;
	
	
	printf("...trying to load gsl");
	gsl_rng *r;
	
    gsl_rng_env_setup();
    double v;
    
    r = gsl_rng_alloc(gsl_rng_default);
		
	printf("Generator type: %s\n", gsl_rng_name(r));
	printf("Seed = %lu\n", gsl_rng_default_seed);
	v= gsl_rng_get(r);
	printf("First value = %.0f\n",v);

	nxsstring nx = "a";
	bool inputfilegiven=false;
	BROWNIE brownie;
	
	// TEST: Use some Rcpp classes
	CharacterVector x(2) ;
	x[0] = "seed"; 
	x[1] =  gsl_rng_name(r);
	
	NumericVector y(2) ;
	y[0] = gsl_rng_default_seed;
	y[1] = v;
	
	List z(2) ; 
	z[0] = x ; 
	z[1] = y ;
	
	gsl_rng_free(r);
	
	return z ;
	
}


