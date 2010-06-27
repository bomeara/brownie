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

int dlInterface::getNumLoadedTrees()
{
	return brownie.intrees.GetNumTrees();
}


std::string dlInterface::getTree(int i)
{
	std::string retstr = (std::string)(*brownie.trees).GetTranslatedTreeDescription(i);
	return retstr;
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


