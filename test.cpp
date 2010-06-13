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

using namespace std;

int main()
{
	#ifndef PROFILE_H
		cout<<"Profile NOT defined."<<endl;
	#endif
	
	#ifndef TREELIB_H
		cout<<"TreeLib NOT defined."<<endl;	
	#endif
	
	// tyring using some important stuff
	Tree none;
	Profile<Tree> shit;
	TaxaBlock *taxa = new TaxaBlock();
	
	gsl_rng *r;
	gsl_rng_env_setup();
	double v;
	r = gsl_rng_alloc(gsl_rng_default);

	
	printf("Generator type: %s\n", gsl_rng_name(r));
	printf("Seed = %lu\n", gsl_rng_default_seed);
	v= gsl_rng_get(r);
	
	
	nxsstring nx = "a";
	bool inputfilegiven=false;
	printf("Initializing BROWNIE object...\n");
	BROWNIE brownie;
	brownie.Init();
	printf("done\n");
	
	// load in text file
	printf("Executing text file...");
	cout << "preload status: "<< brownie.intrees.GetNumTrees()<<endl;
	strcpy(brownie.next_command,"execute parrot.nex\n");
	brownie.PreprocessNextCommand();
	printf("\n .. conditioned command is: %s\n",brownie.next_command);
   	brownie.HandleNextCommand();
	cout << " ... postload status: "<< brownie.intrees.GetNumTrees();
	printf(" ...done\n");
	
	// creating log file:
	cout<<"Setting log file..."<<endl;
	strcpy(brownie.next_command,"log file=parrotlog.txt; echo file=parrotcommand.txt\n");
	brownie.PreprocessNextCommand();
   	brownie.HandleNextCommand();
	
   	// Censored rate test
   	cout<<"Running censored rate test..."<<endl;
	strcpy(brownie.next_command,"ratetest taxset=intrajoint reps=1000 charloop=yes  treeloop=no file=intrajointratetest.txt\n");
	brownie.PreprocessNextCommand();
   	brownie.HandleNextCommand();
   	
	// set model type
	cout<<"Setting log file..."<<endl;
	strcpy(brownie.next_command,"model type=BM1\n");
	brownie.PreprocessNextCommand();
   	brownie.HandleNextCommand();   	
   	
   	
	// non-censored rate test (still not working on windows)
	cout<<"Uncensored rate test..."<<endl;
	strcpy(brownie.next_command,"cont taxset=all treeloop=yes charloop=y file=parrotlog_feedingmodetest.txt\n");
	brownie.PreprocessNextCommand();
   	brownie.HandleNextCommand();   	
   	
   	   	
   	
	printf("First value = %.0f\n",v);
	cout<<"Testing BROWNIE library"<<endl;
	
	gsl_rng_free(r);

	return 0;
}



