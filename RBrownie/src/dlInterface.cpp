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


/* Method which sends a command to the brownie object
 * @author Conrad Stack
 * TODO: handle errors (by capturing cerr?)
 */
void dlPipe(BROWNIE & brau, std::string browniecmd)
{
	
	printf("Piping commands to brownie object");
	strcpy(brau.next_command, browniecmd.c_str());
	brau.PreprocessNextCommand();
	printf("\n .. conditioned command is: %s\n", brau.next_command);
   	brau.HandleNextCommand();
   	
}


