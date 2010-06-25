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
	//int getCharNumber();
	void pipe( std::string );
	void execute(std::string);
	int getNumLoadedTrees();
	std::string getTree(int);
};



#endif
