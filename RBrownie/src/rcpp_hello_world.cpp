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
#include "dlInterface.h"


// TODO: Depreciate this class and replace it with a Module directly to dlInterface

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

	// create interface to brownie object
	dlInterface dli;

	// load in text file
	printf("Executing text file...");
	cout << "preload status: "<< dli.getNumLoadedTrees() << endl;
	
	std::string exstr = "execute parrot.nex\n";
	dli.pipe(exstr);
	
	cout << " ... postload status: "<< dli.getNumLoadedTrees();
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
 * which can be compiled into an Robject of class brownie
 *
 * @author Conrad Stack
 */
SEXP readBrownie(SEXP fnamevect)
{
	using namespace Rcpp;
	
	// Setup interface object
	dlInterface dli;
	
	// Execute the filename
	CharacterVector fname(fnamevect);
	std::string newstr = "" + fname[0];  // convery string_proxy to std::string
	cout << "RCPP: " << newstr << endl;
	dli.execute(newstr);
	
	
	// Number of loaded things:
	int ntrees = dli.getNumLoadedTrees();
	int nchard = dli.getNumDiscreteChars();
	int nchar = dli.getNumChars();
	int ntaxa = dli.getNumTaxa();
	int tset = dli.getNumTaxaSets();
// 	cout<<"Number of trees: "<<ntrees<<endl;
// 	cout<<"Number of taxa: "<< ntaxa <<endl;
// 	cout<<"Number of chars: "<< nchar<<endl;
// 	cout<<"Number of dchars: " << nchard <<endl;
// 	cout<<"Number of taxasets: " << tset <<endl;


	// retrieve TREES
	List treelist(ntrees);
	for (int j=0; j<ntrees; j++)
	{
		cout << "Tree weight: " << dli.getTreeWeight(j) << endl;
		treelist[j] = dli.getTree(j);
	}
	
	// retrieve TAXA
	//cout<<dli.getCharLabels()<<endl;
	std::vector<std::string> taxasets(dli.getTaxaSetNames());
	std::vector< std::vector<std::string> > taxasetfull(dli.getTaxaSets());
	
	std::vector< std::vector<char> > dchars(nchard, std::vector<char>(ntaxa,'A'));
	std::vector< std::vector<float> > cchars(nchar, std::vector<float>(ntaxa,0));
	

	
	
	// output trees to file (regular nexus format)
	// (an alternative way to retrieve info)
	//dli.writeTrees("asdf.txt");
	
	
	return List::create(Named("trees")=treelist,
						Named("taxasetnames")=wrap(taxasets),
						Named("taxasets")=taxasetfull);
	
}
