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
#include <sstream>
#include <iostream>
#include "dl.h"
#include "dlInterface.h"


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
	brownie.Init();

	// load in text file
	printf("Executing text file...");
	cout << "preload status: "<< brownie.intrees.GetNumTrees()<<endl;
	std::string exstr = "execute parrot.nex\n";
	dlPipe(brownie,exstr);
	cout << " ... postload status: "<< brownie.intrees.GetNumTrees();
	printf(" ...done\n");

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

/* Method which executes a brownie file and returns a list
 * which can be construed into an Robject of class brownie
 *
 * @author Conrad Stack
 */
SEXP readBrownie(SEXP fnamevect)
{
	using namespace Rcpp;
	
	// Setup brownie object
	BROWNIE brownie;
	brownie.Init();
	
	// Execute the filename
	CharacterVector fname(fnamevect);
	std::string newstr = EXECUTE + fname[0];
	cout << "RCPP: " << newstr << endl;
	dlPipe(brownie, newstr);
	
	// Retrieve variables (might need to use CharacterVector or StringVector)
	int ntrees = brownie.intrees.GetNumTrees();
	List treelist(ntrees);
	for (int j=0; j<ntrees; j++)
	{
		treelist[j] = (*brownie.trees).GetTranslatedTreeDescription(j);
	}
	
	return treelist;
}
