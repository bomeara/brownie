#include "dlInterface.h"
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
#include <algorithm>

/* Default constructor
 *
 *
 */
dlInterface::dlInterface()
{
	//Initialize brownie object
	brownie.Init();
	
}

/* Default destructor
 */
dlInterface::~dlInterface()
{
	
}


/* Method which sends a command to the brownie object and processes it
 * @author Conrad Stack
 * 
 */
void dlInterface::pipe(std::string browniecmd)
{
	
	printf("Piping commands to brownie object");
	strcpy(brownie.next_command, browniecmd.c_str());
	brownie.PreprocessNextCommand();
	printf("\n .. conditioned command is: %s\n", brownie.next_command);
   	brownie.HandleNextCommand();
   	
}

/* Execute provided nexus file 
 *
 */
void dlInterface::execute(std::string browniefile)
{
	// TODO: check if file exists
	std::string newstr = EXECUTE + browniefile;
	pipe(newstr);
	
}


/* Method which returns the number of TREES which have been loaded
 *
 */
int dlInterface::getNumLoadedTrees()
{
	return brownie.intrees.GetNumTrees();
}

/* Method which returns the number of TAXA which have been loaded
 * @TODO: All taxa need labels, right?  If not then this function
 *			will not work
 */
int dlInterface::getNumTaxa()
{
	return (*brownie.taxa).GetNumTaxonLabels();
}


/* Method which returns the number of CHARACTERS (discrete) which have been loaded
 *
 */
int dlInterface::getNumDiscreteChars()
{
	int retval = 0;
	if (brownie.discretecharloaded)
		retval = (*brownie.discretecharacters).GetNCharTotal();
	
	return retval;
}

/* Method which returns the number of CHARACTERS which have been loaded (number of columns)
 * @TODO: need to make sure that brownie.characters is initialized, 
 *   		otherwise this method will crash the program.
 *
 */
int dlInterface::getNumContinuousChars()
{
	return (*brownie.continuouscharacters).GetNCharTotal();
}


/* Method to return the names of each taxa (if specified)
 * @TODO: Should this be continuous? or Discrete?
*/
std::vector<std::string> dlInterface::getCharLabels(bool cont)
{
	int nchar = 0;
	if(cont)
		nchar = getNumContinuousChars();
	else
		nchar = getNumDiscreteChars();
	
	std::vector<std::string> charnames(nchar);
	
	for(int i=0;i<nchar;i++)
	{
		
		if(cont)
			charnames[i] = (*brownie.continuouscharacters).GetCharLabel(i);
		else
			charnames[i] = (*brownie.discretecharacters).GetCharLabel(i);
	}
	
	return charnames;
}

/* Method to return DISCRETE character states (if they are loaded)
 * also makes the assumption that all taxa have a state to return
 *
 */
std::vector<char> dlInterface::getDiscreteChar(int colindex)
{
	int nchar = getNumDiscreteChars();
	int ntaxa = getNumTaxa();
	std::vector<char> charvals(ntaxa);
	
	// if colindex is within range and there are discrete characters loaded
	//
	if( colindex <= (nchar-1) && nchar > 0 ){
		for(int j =0; j < ntaxa; j++)
		{
			charvals[j] = (*brownie.discretecharacters).GetState(j,colindex);
		}
	}
	
	return charvals;
}


/* Method to return CONTINUOUS character states (if they are loaded)
 * also makes the assumption that all taxa have a state to return
 *
 */
std::vector<float> dlInterface::getContChar(int colindex)
{
	int nchar = getNumContinuousChars();
	int ntaxa = getNumTaxa();
	std::vector<float> fvals(ntaxa);
	
	// if colindex is within range and there are continuous characters loaded
	//
	if( colindex <= (nchar-1) && nchar > 0 ){
		for(int j =0; j < ntaxa; j++)
		{
			fvals[j] = (*brownie.continuouscharacters).GetValue(j,colindex,false);
		}
	}
	return fvals;
	
}


/* Method which returns the number of ASSUMPTIONS which have been loaded
 * (note: there should always be one, unless this method is called before a file is executed)
 * (		that one taxaset which is always present it called ALL, I think it is there by default)
 */
int dlInterface::getNumTaxaSets()
{
	return (*brownie.assumptions).GetNumTaxSets();
}


/* Method which returns taxa set names.  It only shows those (presumably) entered by the user,
 * leaving out those automatically generated by brownie.
 * 
 * @return taxaset names
 * @TODO Make sure that the implicit memory allocation done by push_back plays well with R
 *
 */
std::vector<std::string> dlInterface::getTaxaSetNames()
{
	int ntaxa = getNumTaxaSets();
	vector<string> taxanames;
	
	if(ntaxa > 0)
	{
		LabelList assnames(ntaxa);  // std::vector<nxsstring>
		(*brownie.assumptions).GetTaxSetNames(assnames);
		
		for(int j=0; j < ntaxa; j++)
		{
			
			if(assnames[j].substr(0,3).compare("NOT")!=0 && assnames[j].compare("ALL")!=0)
			{	
				taxanames.push_back((std::string)assnames[j]); // this dynamic allocation might cause problems... hard code the number...
			}
		}
	}

	// TODO: throw error if 
	return taxanames;
}


/* Method which returns taxa sets.  It only shows those (presumably) entered by the user,
 * leaving out those automatically generated by brownie.
 * 
 * @return taxaset names
 * @TODO Make sure that indexing and resizing is done correctly and isn't causing memory leaks
 *
 */
std::vector< std::vector<std::string> > dlInterface::getTaxaSets()
{
	int ntaxa = getNumTaxaSets();
	int currindex = 0; // number that current taxset should be stored at.
	vector< vector<string> > taxasets(ntaxa);
	
	if(ntaxa > 0)
	{
		
		LabelList assnames(ntaxa);  // std::vector<nxsstring>
		(*brownie.assumptions).GetTaxSetNames(assnames);		
		IntSet* isets = new IntSet[ntaxa];  // std::set<int,<less>>
				
		for(int j=0; j < ntaxa; j++)
		{			
			if(assnames[j].substr(0,3).compare("NOT")!=0 && assnames[j].compare("ALL")!=0)
			{
				// IntSets:
				isets[j] = (*brownie.assumptions).GetTaxSet(assnames[j]);
				
				if(!isets[j].empty())
				{
					for(IntSet::iterator kk=isets[j].begin(); kk != isets[j].end(); kk++)
					{
						cout << " " << (*brownie.taxa).GetTaxonLabel(*kk); // *kk is an int, gettaxonlabel returns nxsstring
						taxasets[currindex].push_back((*brownie.taxa).GetTaxonLabel(*kk));
					}
					currindex++;
					
					//cout << isets[j].size() << endl;
				}
			} 
		}		
		delete [] isets;
		
		// resize:
		if( currindex < ntaxa )
			taxasets.resize(currindex);
	}
	return taxasets;
}



/* Method which gets the ith loaded tree from BROWNIE object
 * 
 * @parameter i is the zero-based index of the tree to be returned
 * @parameter translated if true, then get translated trees (trees with actual taxa names, not numbers)
 * @return string object containing the newick-formatted tree.
 *
 */
std::string dlInterface::getTree(int i,bool translated)
{
	std::string retstr;
	
	if(i < getNumLoadedTrees())
	{	
		if(translated)
			retstr = (std::string)(*brownie.trees).GetTranslatedTreeDescription(i);
		else
			retstr = (std::string)(*brownie.trees).GetTreeDescription(i);
	} else {
		// if you are trying to index one larger
		retstr = "NA";
	}
	
	return retstr;
}

/* Overloaded.  Beware, this implementation currently gets trees from a
 * different brownie->object than the other.
 *
 */
bool dlInterface::getTree(int i, std::ostream &f ,bool translated)
{
	bool retval=false;
	
	if(i < getNumLoadedTrees())
	{	
		if(f.good() && translated){
			brownie.intrees.GetIthTree(i).WriteNoQuote(f);
			retval = true; //TODO: use Tree.GetError() to validate
		}
	}
	
	return retval;
}


/* Get tree weight (should come from comment in nexus file) [&W=...]
 *
 * 
 */
float dlInterface::getTreeWeight(int i)
{
	return (*brownie.trees).GetTreeWeight(i);
}


/* Write trees to file (should be temporary file)
 *
 * @author Conrad Stack
 * @parameter outfile the file to write to
 * @return whether or not file was written to.
 *
 */
bool dlInterface::writeTrees(std::string outfile)
{
	bool retbool = false;
	if(getNumLoadedTrees() >= 1)
	{
		std::ofstream output;
		try
		{
			output.open(outfile.c_str());
			retbool = brownie.intrees.WriteTrees(output);
			output.close();
		} 
		catch (std::ios_base::failure const& ex) 
		{
			if( output.is_open() )
				output.close();
			
			std::cerr << "dlInterface::writeTrees -> I/O error.  Could not write to file " << outfile << std::endl;
			throw;
		}
	}
	
	return retbool;
}


int dlInterface::getNumRetTrees()
{
	// only allow 1 for now:
	int maxtree=1;
	return( min( (int)brownie.rettree.str().length() , maxtree) );
}


/* Check if return tree stream has any data:
 *
 */
bool dlInterface::hasRetTrees()
{
	return(getNumRetTrees()==0);
}


/* Get reconstructed tree, if there is one available
 *
 *
 */
std::string dlInterface::getRetTree(int index)
{
	std::string retstr;
	// force only one for now:
	index = 0;
	
	if(hasRetTrees())
	{
		retstr = brownie.rettree.str();
	}
	return (retstr);
	
}

