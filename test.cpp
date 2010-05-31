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


// globals and structures


// function declarations


//template <> class Profile<Tree>;

using namespace std;

int main()
{
	
    const gsl_rng_type * tt;
    gsl_rng_env_setup();
    tt = gsl_rng_mt19937;
    r = gsl_rng_alloc(tt);
    
	cout<<"Testing BROWNIE library"<<endl;
	
	return 0;
}
