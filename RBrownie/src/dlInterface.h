#ifndef _RBrownie_DLINTERFACE_H
#define _RBrownie_DLINTERFACE_H

#include "dl.h"

const std::string EXECUTE = "execute ";

class dlInterface
{
	private:
	BROWNIE brownie;
	
	
	public:
	dlInterface();
	~dlInterface();
	void pipe( std::string );
	void execute(std::string);
	int getNumLoadedTrees();
	int getNumTaxa();
	int getNumDiscreteChars();
	int getNumChars();
	std::vector<std::string> getCharLabels();
	std::vector<char> getDiscreteChar(int colindex=0);
	std::vector<float> getContChar(int colindex=0);
	//TODO: add indexing by character name?
		
	int getNumTaxaSets();
	std::vector<std::string> getTaxaSetNames();
	std::vector< std::vector<std::string> > getTaxaSets();
	
	std::string getTree(int,bool=true);  // these two might not return the same trees.
	bool getTree(int, std::ostream &f, bool=true);
	float getTreeWeight(int);
	bool writeTrees(std::string);
};



#endif
