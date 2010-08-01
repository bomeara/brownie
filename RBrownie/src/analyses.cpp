#include "analyses.h"
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

/*
 // Might use this function in the future:
RcppExport SEXP RunTest(SEXP fnamevect, bool retReturnTrees, bool retContinuous, bool retDiscrete, bool retTrees, bool retTaxa, bool retTaxasets)
{
	using namespace Rcpp;
	
	dlInterface dli;
	int nret = 0;
	int retcount = 0;
	
	// Certain number of objects to return:
	//
	nret += int(retReturnTrees);
	nret += int(retContinuous);
	nret += int(retDiscrete);
	nret += int(retTrees);
	nret += int(retTaxa);
	nret += int(retTaxasets);
	List returnlist(nret);
	
	// Return storage:
	
	
	// Check if file exists:
	CharacterVector fname(fnamevect);
	std::string filestr = "" + fname[0];  // convery string_proxy to std::string
	
	// Execute:
	dli.execute(filestr);
	
	// Retrieve data:
	int ntrees = dli.getNumLoadedTrees();
	int nchard = dli.getNumDiscreteChars();
	int nchar = dli.getNumContinuousChars();
	int ntaxa = dli.getNumTaxa();
	int tset = dli.getNumTaxaSets();
	int nrettrees = dli.getNumRetTrees();
	
	if(retTrees)
	{
		List treelist(ntrees);
		for (int j=0; j<ntrees; j++)
		{
			//cout << "Tree weight: " << dli.getTreeWeight(j) << endl;
			treelist[j] = dli.getTree(j);
		}
		returnlist[retcount++] = treelist;
	}
	
	if(retTaxasets)
	{
		std::vector<std::string> taxasets(dli.getTaxaSetNames());
		std::vector< std::vector<std::string> > taxasetfull(dli.getTaxaSets());		
		
		returnlist[retcount++] = List::create(Named("taxasetnames")=wrap(taxasets),Named("taxasets")=taxasetfull);
	}
	
	// Return trees (for ASR)
	if(retReturnTrees)
	{
		List rtreelist(nrettrees);
		
		// only one available for now, expand this later
		//
		if(nrettrees > 0)
		{
			rtreelist[0] = dli.getRetTree();
		}
		
		returnlist[retcount++] = rtreelist;
	}
	
	if(retDiscrete)
	{
		std::vector< std::vector<char> > dchars(nchard);
		
		// load discrete characters
		for(int ii = 0; ii < nchard; ii++)
		{
			dchars[ii] = dli.getDiscreteChar(ii);
		}
		returnlist[retcount++] = dchars;
	}
	
	if(retContinuous)
	{
		std::vector< std::vector<float> > cchars(nchar);
		// load discrete characters
		for(int ii = 0; ii < nchar; ii++)
		{
			cchars[ii] = dli.getContChar(ii);
		}
		returnlist[retcount++] = cchars;
	}
		
	if(retTaxa)
	{
		// TODO
		//returnlist[retcount++]
	}
	
	return returnlist;
}
*/


/* Run any analysis and hope that something is returned:
 */
SEXP doAnalysis(SEXP fnamevect)
{
	using namespace Rcpp;
	
	// Setup interface object
	dlInterface dli;
	
	// Execute the filename
	CharacterVector fname(fnamevect);
	std::string newstr = "" + fname[0];  // convery string_proxy to std::string
	cout << "RCPP: " << newstr << endl;
	dli.execute(newstr);
		
	List returnlist = List::create(Named("textout")=wrap(dli.getReturnStrings()),
						Named("treesout")=wrap(dli.getReturnTrees()));	
	
	return returnlist;
}



