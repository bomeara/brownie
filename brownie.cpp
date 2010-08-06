#include <strstream>
#include <fstream.h>
#include <iomanip.h>
#include <unistd.h>
#include <stdio.h>
#include <set>

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
#include "superdouble.h"

//took out this section since GTP is built in
//extern "C" {
//#include "gtp.c"
//    static long gtp(node n);
//    static long strongDup(node n);
//    static long numStrongDup(node n);
//    static void initgTree(node n, Hash h);
//    static void makeTipNameHashHelper(node n);
//    Hash makeTipNameHash(node root);
//    static void initAncArray(node n);
//    double ReturnScore(char *buffer, int unrooted);
//    void printUsage(void);
//}

using std::string;
using std::ostringstream;
gsl_rng * r;  /* global generator */


/**
*  Brownie
 *  Copyright Brian O'Meara, 2006
 *  http://wwww.brianomeara.info/brownie
 *  Released under the GNU Public License
 *  There is no warranty. However, please let
 *  me know about any bugs, and I'll try to
 *  squash them.
 */

/**
* This file is a heavily modified version of the
 * basicccmdline.cpp file that comes with the
 * Nexus class library.
 * Brian O'Meara. February 2006.
 */

/**
* This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */



/**
* @class      BROWNIE
 * @file       brownie.h
 * @file       brownie.cpp
 * @variable   command_maxlen [static int] maximum length of a command entered
 */

/**
* @constructor
 *
 * Initializes the id data member to "BROWNIE" and calls the FactoryDefaults
 * member function to perform the remaining initializations.
 */
BROWNIE::BROWNIE()
{
    id = "BROWNIE";
    FactoryDefaults();
}

/**
* @destructor
 *
 * Closes logf or echof if it is open.
 */
BROWNIE::~BROWNIE()
{
    if( logf_open )
        logf.close();
    if (echof_open) {
        echof<<"end;";
        echof.close();
    }
}

/**
* @method CharLabelToNumber [int:protected]
 * @param s [nxsstring] the character label to be translated to character number
 *
 * The code here is identical to the base class version (simply returns 0),
 * so the code here should either be modified or this derived version
 * eliminated altogether.  Under what circumstances would you need to
 * modify the default code, you ask?  This function should be modified
 * to something meaningful if this derived class needs to construct and
 * run a SetReader object to read a set involving characters.  The SetReader
 * object may need to use this function to look up a character label
 * encountered in the set.  A class that overrides this method should
 * return the character index in the range [1..nchar]; i.e., add one to the
 * 0-offset index.
 */
int BROWNIE::CharLabelToNumber( nxsstring /*s*/ )
{
    return 0;
}

/**
* @method EnteringBlock [virtual void:public]
 * @param blockName [nxsstring] the name of the block just entered
 *
 * Called by the Nexus object when a block named blockName is entered.
 * Allows program to notify user of progress in parsing the NEXUS file.
 * Virtual function that overrides the pure virtual function in the
 * base class Nexus.
 */
void BROWNIE::EnteringBlock( nxsstring blockName )
{
    message = "Reading ";
    message += blockName;
    message += " block...";
    PrintMessage();
}

/**
* @method ExitingBlock [virtual void:public]
 * @param blockName [nxsstring] the name of the block just exited
 *
 * Called by the Nexus object when exiting a block named blockName.
 * Allows program to notify user of progress in parsing the NEXUS file.
 * Virtual function that overrides the pure virtual function in the
 * base class Nexus.
 */
void BROWNIE::ExitingBlock( nxsstring /*blockName*/ )
{
}

/**
* @method FactoryDefaults [void:protected]
 *
 * Sets all data members to their factory default settings:
 * <table>
 * <tr><th align="left">Variable <th> <th align="left"> Initial Value
 * <tr><td> inf_open          <td>= <td> false
 * <tr><td> logf_open         <td>= <td> false
 * <tr><td> quit_now          <td>= <td> false
 * <tr><td> message           <td>= <td> ""
 * <tr><td> next_command[0]   <td>= <td> '\0'
 * <tr><td> trees             <td>= <td> NULL
 * <tr><td> taxa              <td>= <td> NULL
 * <tr><td> assumptions       <td>= <td> NULL
 * <tr><td> characters        <td>= <td> NULL
 * </table>
 */
void BROWNIE::FactoryDefaults()
{
    inf_open = false;
    logf_open = false;
    echof_open=false;
    quit_now = false;
	quit_onerr = false;
    message = "";
    next_command[0] = '\0';
    gslseedtoprint=time(NULL);
    gslseed=gslseedtoprint;
    gsl_rng_set(r,gslseed);
    trees = NULL;
    taxa = NULL;
    assumptions = NULL;
    characters = NULL;
    discretecharacters = NULL;
    discretecharloaded = false;
    chosenchar=1;
	discretechosenchar=0; //starts at index=0;
	ratematfixedvector.push_back(0); 
	ratematassignvector.push_back(0); 
	freerateletterstring="";
	negbounceparam=-1;
	nonnegvariables=true;
	gsl_vector *userstatefreqvector=gsl_vector_calloc(1);
	gsl_matrix *optimaldiscretecharQmatrix=gsl_matrix_calloc(1,1);
	gsl_vector *optimaldiscretecharstatefreq=gsl_vector_calloc(1);
	gsl_matrix *currentdiscretecharQmatrix=gsl_matrix_calloc(1,1);
	gsl_vector *currentdiscretecharstatefreq=gsl_vector_calloc(1);	
	discretechosenmodel=1;
	bestdiscretelikelihood=GSL_POSINF;
	optimizationalgorithm=1;
	discretechosenstatefreqmodel=3;
	allchar=false;
	globalstates=false;
	variablecharonly=false;
	numbercharstates=0;
	localnumbercharstates=0;
	numberoffreeparameters=0;
	numberoffreerates=0;
	numberoffreefreqs=0;
    chosentree=1;
    tipvariancetype=0;
    progressbartotal=0;
    progressbarcount=0;
    progressbarprinted=0;
    debugmode=false;
    maxiterations=1000;
    stoppingprecision=1e-7;
    randomstarts=15;
    treefilename="besttrees.tre";
	useCOAL=false;
	useMS=false;
	msbasereps=100000;
	contourBrlenToExport=2;
	contourMaxRecursions=20;
	contourstartingnumbersteps=15;
	contourstartingwidth=10;
	exportalltrees=true;
	COALaicmode=4;
	markedmultiplier=5.0;
	brlensigma=1.0;
	numbrlenadjustments=20;
    badgtpcount=0;
    stepsize=.1;
	npercent=0.95;
    detailedoutput=false;
	redobad=false;
	giveupfactor=100;
    outputallatonce=true;
    chosenmodel=1;
	gsl_vector *optimalvaluescontinuouschar = gsl_vector_calloc(1);
	gsl_matrix *optimalVCV = gsl_matrix_calloc(1,1);
	//optimalvalueslabels.clear();
	gsl_vector *optimalTraitMeans = gsl_vector_calloc(1);
	globalchosentaxset="ALL";
    citationarray[0]=true;
    maxnumspecies=100;
    minnumspecies=1;
	minsamplesperspecies=3;
	maxstartstops=10;
	rearrlimit=-1;
    steepest=false;
	exhaustive=false;
    status=true;
	jackknifesearch=false;
    nreps=5;
	jreps=100;
	pctdelete=1.0/3.0;
	jackrep=0;
    showtries=false;
    structwt=0.5;
	triplettoohigh=false;
	gtptoohigh=false;
	infinitescore=false;
	tripletdistthreshold=0.2; //Sets how often to use NJ tree distances for starting assignments (higher number=more often) and how often to use triplet support
    pthreshold=1;
	chosensubsampling=2.0;
    movefreqvector.push_back(0.80); //initial values of the movefreqvector; set up so first try doing swaps and rerootings, then reassignments
    movefreqvector.push_back(0.01);
    movefreqvector.push_back(0.01);
    movefreqvector.push_back(0.01);
    movefreqvector.push_back(0.17);
	movefreqvector.push_back(0.0); //do no brlen optimization
    sppnumfixed=false;
    unrooted=0; //by default, use rooted gene trees
    bestscore=GSL_POSINF;
	bestscorelocal=GSL_POSINF;
    for (int arrayloc=1;arrayloc<20;arrayloc++) {
        citationarray[arrayloc]=false;
    }
    for (int state=0;state<=(maxModelCategoryStates-1);state++) {
        staterestrictionvector.push_back(state);
        timeslicetimes.push_back(GSL_POSINF);
        timeslicemodels.push_back(0);
    }
    //timeslicetimes[maxModelCategoryStates]=GSL_POSINF; //Use the BMC and BMS code; the size of this vector is the maximum number of time intervals minus 1, and the max num of time intervals is the max number of models
    //timeslicemodels[maxModelCategoryStates]=0; //Since we use BMC and BMS code, the size of this vector is the maximum number of model categories
	CDFvectorholder bob;
	CDFvector=bob.Initialize();
}
/**
* @method FileExists [bool:protected]
 * @param fn [const char*] the name of the file to check
 *
 * Returns true if file named fn already exists,
 * false otherwise.
 */
bool BROWNIE::FileExists( const char* fn )
{
    bool exists = false;

    //	if( access( fn, 0 ) == 0 )
    //		exists = true;

    FILE* fp = fopen( fn, "r" );
    if( fp != NULL ) {
        fclose(fp);
        exists = true;
    }

    return exists;
}

/**
* @method GetFileName [nxsstring:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called whenever a file name needs to be read from either
 * the command line or a file.  Expects next token to be "="
 * followed by the token representing the file name.  Call
 * this function after, say, the keyword "file" has been
 * read in the following LOG command:
 * <pre>
 * log file=doofus.txt start replace;
 * </pre>
 * Note that this function will read only "=doofus.txt "
 * leaving "start replace;" in the stream for reading
 * at a later time.
 */
nxsstring BROWNIE::GetFileName( NexusToken& token )
{
    // Eat the equals sign
    //
    token.GetNextToken();

    if( !token.Equals("=") ) {
        errormsg = "Expecting an equals sign, but found ";
        errormsg += token.GetToken();
        errormsg += " instead";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    // Now get the filename itself
    //
    token.GetNextToken();

    return token.GetToken();
}


gsl_rng * BROWNIE::ReturnR()
{
    return r;
}

void BROWNIE::ProgressBar(int total)
{
    if (total>0) { //so to start it, give a value of total >0 (# reps); after that, feed it ProgressBar(0);
        progressbartotal=total;

        cout<<"\nProgress:\n0%     10%     20%     30%     40%     50%     60%     70%     80%     90%     100%\n|"<<flush;
        //cout<<"|....|....|....|....|....|....|....|....|....|....|"
    }
    else {
        progressbarcount++;
        double sampleratio=(1.0*progressbarcount)/(1.0*progressbartotal); // convert to floating point division
        double printratio=progressbarprinted/80.0;
        //cout<<sampleratio<<" "<<printratio<<endl;
        while (sampleratio>printratio) {
            cout<<"*"<<flush;
            progressbarprinted++;
            printratio=(1.0*progressbarprinted)/80.0;
        }
        if (progressbarcount==progressbartotal) { //stop it, reinitialize
            cout<<"|\n\n"<<flush;
            progressbarcount=0;
            progressbartotal=0;
            progressbarprinted=0;
        }
    }
}


/**
* @method GetNumber [nxsstring:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called whenever a number needs to be read from either
 * the command line or a file.  Expects next token to be "="
 * followed by the token representing the number.  Call
 * this function after, say, the keyword "tree" has been
 * read in the following LOG command:
 * <pre>
 * choose tree=5 start replace;
 * </pre>
 * Note that this function will read only "=5 "
 * leaving "start replace;" in the stream for reading
 * at a later time.
 */
nxsstring BROWNIE::GetNumber( NexusToken& token )
{
    // Eat the equals sign
    //
    token.GetNextToken();

    if( !token.Equals("=") ) {
        errormsg = "Expecting an equals sign, but found ";
        errormsg += token.GetToken();
        errormsg += " instead";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    // Now get the number itself
    //
    token.GetNextToken();

    return token.GetToken();
}

/**
* @method GetNumberOnly [nxsstring:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called whenever a number needs to be read from either
 * the command line or a file.  Does not expext an equals
 */
nxsstring BROWNIE::GetNumberOnly( NexusToken& token )
{

    token.GetNextToken();
    return token.GetToken();
}

/**
* @method HandleEndblock [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the END or ENDBLOCK command needs to be parsed
 * from within the BROWNIE block.  Basically just checks to make
 * sure the next token in  the data file is a semicolon.
 */
void BROWNIE::HandleEndblock( NexusToken& token )
{
    // get the semicolon following END or ENDBLOCK token
    //
    token.GetNextToken();

    if( !token.Equals(";") ) {
        errormsg = "Expecting ';' to terminate the END or ENDBLOCK command, but found ";
        errormsg += token.GetToken();
        errormsg += " instead";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }
}

/**
* @method HandleExecuteCmdLine
 * used to execute a datafile from the command line when starting the program
 */
void BROWNIE::HandleExecuteCmdLine(nxsstring fn)
{
    if( FileExists( fn.c_str() ) )
    {
        cerr << endl;
        cerr << "Opening " << fn << "..." << endl;

        ifstream inf( fn.c_str(), ios::binary | ios::in );

        inf_open = true;

        NexusToken ftoken(inf);

        inf_open = true;
        ifstream intreefile;
        intreefile.open(fn.c_str(),ios::in);
        if (!intrees.ReadTrees(intreefile))
        {
            message="No trees read from file\n";
            PrintMessage();
        }
        intreefile.close();

        try {
            Execute( ftoken );
        }
        catch( XNexus x )
        {
            NexusError( errormsg, x.pos, x.line, x.col );
            Reset();
        }
		assumptions->MakeTaxsetAll();


        //  if( !taxa->IsEmpty() ) {
        //      cerr << "  TAXA block found" << endl;
        //      if( logf_open )
        //          taxa->Report(logf);
        //  }

        if( inf_open )
            inf.close();
        inf_open = false;

        //     if( !trees->IsEmpty() ) {
        //         cerr << "  TREES block found" << endl;
        //         if( logf_open )
        //            trees->Report(logf);

		//need to call this here in case not already set by reading  Brownie block
		if(!characters->IsEmpty() ) {
			if (characters->GetDataType()==6) {
				continuouscharacters=characters;
				//cout<<"Found continuous characters\n";
			}
			else {
				discretecharacters=characters;
				numbercharstates=discretecharacters->GetMaxObsNumStates();	
				localnumbercharstates=numbercharstates;
				discretecharloaded=true;
				//cout<<"Found discrete characters\n";
			}
		}
		
		if(!characters2->IsEmpty() ) {
			if (characters2->GetDataType()==6) {
				continuouscharacters=characters2;
				//cout<<"Found continuous characters\n";
			}
			else {
				discretecharacters=characters2;
				numbercharstates=discretecharacters->GetMaxObsNumStates();	
				localnumbercharstates=numbercharstates;
				discretecharloaded=true;
				//cout<<"Found discrete characters\n";
			}
		}

        //     }



	
        //        cerr << "  ASSUMPTIONS block found" << endl;
        //        if( logf_open )
        //             assumptions->Report(logf);
        //    }

        //    if( !characters->IsEmpty() ) {
        //        cerr << "  CHARACTERS block found" << endl;
        //        if( logf_open )
        //            characters->Report(logf);
        //    }
    }
    else
    {
        cerr << endl;
        cerr << "Oops! Could not find specified file: " << fn << endl;
    }
}



/**
* @method HandleExecute [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Handles everything after the EXecute keyword and the terminating
 * semicolon.  Flushes all blocks before executing file specified,
 * and no warning is given of this
 */
void BROWNIE::HandleExecute( NexusToken& token )
{
    // Issuing the EXECUTE command from within a file is a no-no
    //
    if( inf_open ) {
        errormsg = "Cannot issue execute command from within a BROWNIE block";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    // Get the file name to execute
    //
    token.GetNextToken();
    if( token.Equals("?") ) {
        message="Usage: Exe <file-name>\n\n";
        message+="Executes nexus-formatted input file.\nThis file may contain continuous character data, trees, and Brownie blocks.\nThe file name, if put in single quotes, may contain the path to the file.\n\n";
        PrintMessage();
        token.GetNextToken(); //to remove the remaining token
    }
    else {
        nxsstring fn = token.GetToken();

        // get the semicolon terminating the EXECUTE command
        //
        token.GetNextToken();

        if( !token.Equals(";") ) {
            errormsg = "Expecting ';' to terminate the EXECUTE command, but found ";
            errormsg += token.GetToken();
            errormsg += " instead";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }

        // Before going through with this, make sure we're not going to overwrite
        // any stored blocks
        bool stuff_stored = !taxa->IsEmpty();
        stuff_stored = ( stuff_stored || !trees->IsEmpty() );
        stuff_stored = ( stuff_stored || !assumptions->IsEmpty() );
        stuff_stored = ( stuff_stored || !characters->IsEmpty() );
        if( stuff_stored && UserSaysOk( "Ok to delete?", "Data has already been read and stored" ) )
            PurgeBlocks();
        else if( stuff_stored ) {
            message = "\nExecute command aborted.";
            PrintMessage();
            return;
        }

        if( FileExists( fn.c_str() ) )
        {
            cerr << endl;
            cerr << "Opening " << fn << "..." << endl;

            ifstream inf( fn.c_str(), ios::binary | ios::in );

            inf_open = true;
            ifstream intreefile;
            intreefile.open(fn.c_str(),ios::in);
            if (!intrees.ReadTrees(intreefile))
            {
                message="No trees read from file\n";
                PrintMessage();
            }
            intreefile.close();

            NexusToken ftoken(inf);

            try {
                Execute( ftoken );
            }
            catch( XNexus x )
            {
                NexusError( errormsg, x.pos, x.line, x.col );
                Reset();
            }
		//need to call this here in case not already set by reading  Brownie block
			assumptions->MakeTaxsetAll();
			if(!characters->IsEmpty() ) {
				if (characters->GetDataType()==6) {
					continuouscharacters=characters;
					//cout<<"Found continuous characters\n";
				}
				else {
					discretecharacters=characters;
					numbercharstates=discretecharacters->GetMaxObsNumStates();	
					localnumbercharstates=numbercharstates;
					discretecharloaded=true;
					//cout<<"Found discrete characters\n";
				}
			}
			
			if(!characters2->IsEmpty() ) {
				if (characters2->GetDataType()==6) {
					continuouscharacters=characters2;
					//cout<<"Found continuous characters\n";
				}
				else {
					discretecharacters=characters2;
					numbercharstates=discretecharacters->GetMaxObsNumStates();	
					localnumbercharstates=numbercharstates;
					discretecharloaded=true;
					//cout<<"Found discrete characters\n";
				}
			}
			

            //    if( !taxa->IsEmpty() ) {
            //        cerr << "  TAXA block found" << endl;
            //        if( logf_open )
            //            taxa->Report(logf);
            //    }

            if( inf_open )
                inf.close();
            inf_open = false;

            //    if( !trees->IsEmpty() ) {
            //        cerr << "  TREES block found" << endl;
            //        if( logf_open )
            //             trees->Report(logf);




            //   }


            // if( !assumptions->IsEmpty() ) {
            //   cerr << "  ASSUMPTIONS block found" << endl;
            // if( logf_open )
            //   assumptions->Report(logf);
            //    }

            //  if( !characters->IsEmpty() ) {
            //       cerr << "  CHARACTERS block found" << endl;
            //       if( logf_open )
            //           characters->Report(logf);
            //   }
        }
        else
        {
            cerr << endl;
            cerr << "Oops! Could not find specified file: " << fn << endl;
        }
    }
}

/**
* @method HandleHelp [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the HELP command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleHelp( NexusToken& token )
{
    // Retrieve all tokens for this command, stopping only in the event
    // of a semicolon or an unrecognized keyword
    //
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading HELP command.\n\nTry typing \"";
            errormsg +=token.GetToken();
            errormsg +=" ?\" instead [without the quotes] if you want help for the command.";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }

	message = "\nGeneral commands:";
    message += "\n  help             -> shows this message";
    message += "\n  exe              -> executes nexus file";
    message += "\n  log              -> log output";
    message += "\n  echo             -> copies your commands into a batch file";
    //message += "\n  gettree          -> loads tree file";
    message += "\n  blocks           -> reports on blocks currently stored";
    message += "\n  showtree         -> displays currently loaded tree(s)";
	message += "\n  choose           -> chooses tree or char for analysis";
    message += "\n  taxset           -> stores a taxset";
	message += "\n  citation         -> outputs list of relevant papers for your analyses";
	message += "\n  tipvalues        -> return list of tip values";
	message += "\n  quit             -> terminates application";
	message += "\n\nCharacter evolution:";
    message += "\n  ratetest         -> does censored rate test (original Brownie function)";
    message += "\n  vcv              -> outputs a variance-covariance matrix";
	message += "\n  (discrete)       -> implements discrete character models and reconstructions";
    message += "\n  [tipvariance]    -> allows program to deal with variance in taxon means";
    message += "\n  model            -> sets model of continuous character evolution (OU, BM, etc)";
    message += "\n  (continuous)     -> gets score for chosen taxset for chosen model";
    message += "\n  [export]         -> exports a tree and data in deprecated Pagel format";
	message += "\n  (simulate)       -> simulate discrete or continuous character matrices";
	message += "\n  loss             -> estimate rates of binary character loss on branches";
	message += "\n\nSpecies delimitation and  tree search:";
	message += "\n  hs               -> perform a heuristic search";
	//message += "\n  [jackknife]      -> perform a jackknife search";
	//message += "\n  [exhaustive]     -> perform an exhaustive search";
	message += "\n  compare          -> compare triplet overlap for coalescent trees";
	message += "\n  assign           -> assign samples to species";
	message += "\n  accuracy         -> compute accuracy of reconstruction";
	message += "\n\nNumerical optimization settings:";
	message += "\n  set              -> sets options";
    message += "\n  numopt           -> sets parameters for numerical optimization functions";	
	message += "\n\nMiscellaneous:";
	message += "\n  orderbytree      -> reorders a datamatrix by order of taxa in a tree";
	message += "\n  printedgelength  -> prints branch lengths";
	message += "\n  (partitionededge)-> outputs all trees one NNI move away for NNIBS analysis";
	message += "\n\nIn development:";
    message += "\n  [nast]";	
	message += "\n  [timeslice]";
	message += "\n  [debug]";
	message += "\n  [Garland]";

    message += "\n\nType \"commandname ?\" [without the quotes]\nfor help on any command.\n\n*** IMPORTANT ***\nCommands in brackets (\"[]\") should not be used for published results yet\nCommands in parentheses (\"()\") should be used after checking with Brian O'Meara (omeara.brian@gmail.com) -- they may reflect things in review which might change pending reviewers' comments, for example";
    PrintMessage();
}

/**
* @method HandleBlocks [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the Blocks command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleBlocks( NexusToken& token )
{
    bool finishexecuting=true;
    // Retrieve all tokens for this command, stopping only in the event
    // of a semicolon or an unrecognized keyword
    //
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else if ( token.Equals("?") ) {
            message="Usage: Blocks\n\n";
            message+="Reports status of loaded nexus blocks\n";
            PrintMessage();
            finishexecuting=false;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading Blocks command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
    if (finishexecuting) {
        message = "\nNexus blocks currently stored:";
        PrintMessage();
        if( !taxa->IsEmpty() ) {
            cerr << "\n  TAXA block found" << endl;
            taxa->Report(cerr);
            if( logf_open )
                taxa->Report(logf);
        }

        if( !trees->IsEmpty() ) {
            cerr << "\n  TREES block found" << endl;
            trees->Report(cerr);
            if( logf_open )
                trees->Report(logf);
        }
        if( !assumptions->IsEmpty() ) {
            cerr << "\n  ASSUMPTIONS block found" << endl;
            assumptions->Report(cerr);
            if( logf_open )
                assumptions->Report(logf);
        }

        if( !characters->IsEmpty() ) {
            cerr << "\n  CHARACTERS block found" << endl;
            characters->Report(cerr);
            if( logf_open )
                characters->Report(logf);
        }
		
		if( !characters2->IsEmpty() ) {
            cerr << "\n  CHARACTERS2 block found" << endl;
            characters2->Report(cerr);
            if( logf_open )
                characters2->Report(logf);
        }
		
    }
}


/**
* @method HandleDebug [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 */
void BROWNIE::HandleDebug( NexusToken& token )
{
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            debugmode=true;
            break;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading Debug command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
    message="Now entering debug mode\n";
    PrintMessage();
}


/**
* @method HandleNoQuitOnErr [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 */
void BROWNIE::HandleNoQuitOnErr( NexusToken& token )
{
    for(;;)
    {
        token.GetNextToken();
		
        if( token.Equals(";") ) {
            quit_onerr=false;
            break;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading NoQuitOnErr command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
    message="Will not quit on error\n";
    PrintMessage();
}

/**
* @method HandleNoQuitOnErr [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 */
void BROWNIE::HandleQuitOnErr( NexusToken& token )
{
    for(;;)
    {
        token.GetNextToken();
		
        if( token.Equals(";") ) {
            quit_onerr=true;
            break;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading QuitOnErr command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
    message="Will instantly quit on error\n";
    PrintMessage();
}

double BROWNIE::AIC(double neglnL, int K) {
	double AIC = 2.0 * (1.0*K+neglnL);
	return AIC;
}

double BROWNIE::AICc(double neglnL, int K, int N) {
	double AICc=GSL_POSINF;
	if ((N-K-1)>0) {
		AICc= 2.0 * (1.0*K+neglnL) + 2.0*K*(K+1.0)/(1.0*(N-K-1.0));
	}
	return AICc;
}



/**
* @method HandleGettrees [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the GETTREES command needs to be parsed
 * from within the BROWNIE block.
 *
 * If there are trees already in the executed Nexus block, no need
 *   to run this.
 */
void BROWNIE::HandleGettrees( NexusToken& token )
{
    if( !trees->IsEmpty() ) {
        if( UserSaysOk( "Ok to delete?", "Trees have already been read and stored" ) ) {
            Detach( trees );
            delete trees;
            trees = new TreesBlock(*taxa);
            Add( trees );
        }
        else {
            message = "\nGetTrees command aborted.";
            PrintMessage();
            return;
        }
    }

    nxsstring fn;
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else if( token.Abbreviation("File") ) {
            fn=GetFileName(token);
            break;
        }
        else {
            fn=token.GetToken();
        }
    }
    //input stream
    ifstream intreefile;
    intreefile.open(fn.c_str(),ios::in);




    if (!intrees.ReadTrees(intreefile))
    {
        errormsg="Failed to read trees";
        throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }
    intreefile.close();
}

/**
* @method HandleShowtree [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the SHOWTREE command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleShowtree( NexusToken& token )
{
    nxsstring numbernexus;
    bool finishexecuting=true;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (finishexecuting) {
                intrees.ShowTrees(cout);
            }
            break;
        }
        else if( token.Abbreviation("?") ) {
            message="Usage: ShowTree\nDisplays all loaded trees. No option yet for just displaying the chosen tree.";
            PrintMessage();
            finishexecuting=false;
        }
    }
}


/**
* @method HandleHeuristicSearch [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 */
void BROWNIE::HandleHeuristicSearch( NexusToken& token )
{
	citationarray[3]=true;
    nxsstring numbernexus;
    bool finishexecuting=true;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (finishexecuting) {
				if(exhaustive) {
					DoExhaustiveSearch();
				}
				else if (!jackknifesearch) {
					DoHeuristicSearch();
				}
				else {
					ofstream jacktreef;
					nxsstring jacktreename=treefilename;
					jacktreename+=".jack.tre";
					jacktreef.open(jacktreename.c_str());
					jacktreef<<"#nexus\nbegin trees;\n"; //change these couts to jacktreef
					jacktreef.close();
					for (jackrep=1;jackrep<=jreps;jackrep++) { //a jackknife search consists of several heuristic searches, with the weights of several input trees set to zero
																//Each gene tree's weight is used with the jackknife proportion to get a deletion vector; using the weight ensures that the genes have the same expectation of weighted representation in the final set (think of combining  bootstrap samples for one gene with regular samples for another).
//NOTE: Need to deal with tree weights properly in GetGTPScoreNew; have not added anything to GetTripletScore; need to deal with how tree weights are used to calculate when trees come from different genes [perhaps add a new weight vector?]

//Use this in place of the original for ..intrees to only deal with selected trees
//	for (int i=0;i<intrees.GetNumTrees();i++) {
//		bool usethistree=true;
//		if (jackknifesearch) {
//			if (jackknifevector[i]==0) {
//				usethistree=false;
//			}
//		}
//		if (usethistree) {
//          STUFF originally just under the intrees for loop
//      }
//  }
						DoHeuristicSearch();	
						jacktreef.open(jacktreename.c_str(), ios::out | ios::app);
						jacktreef<<jackknifetreestooutput;
						jacktreef.close();

					}
					jacktreef.open(jacktreename.c_str(), ios::out | ios::app);
					jackknifesearch=false;
					jacktreef<<"end;";
					jacktreef.close();
				}
            }
            break;
        }
        else if( token.Abbreviation("?") ) {
			if (!jackknifesearch) {
				message="Usage: HSearch [options]\nDoes a heuristic search ";
			}
			else {
				message="Usage: Jackknife [options]\nDoes a jackknife search, deleting input trees ";
			}
            message+="Available options:\n\n";
            message+="Keyword ---- Option type --------------------------- Current setting --";
			if (jackknifesearch) {
				message+="\nJReps        <integer-value>                         ";
				message+=jreps;
				message+="\nPctDelete    <real-value>                            ";
				message+=pctdelete;
			}
            message+="\nNReps        <integer-value>                         ";
            message+=nreps;
            //message+="\nTimeLimit    <integer-value>|None                    None";
            //message+="\nClock        Wall|CPU                                Wall";
            message+="\nRearrLimit   <integer-value>|None                    ";
			if (rearrlimit<0) {
					message+="None";
			}
			else {
				message+=rearrlimit;
			}
            //message+="\nAssignFixed  No|Yes                                  Yes";
            //message+="\nSppNumFixed  No|Yes                                  Yes";
            message+="\nMaxNumSpp    <integer-value>                         ";
            message+=maxnumspecies;
            message+="\nMinNumSpp    <integer-value>                         ";
            message+=minnumspecies;
			message+="\nMinSamp      <integer-value>                         ";
            message+=minsamplesperspecies;
            message+="\nStructWt     <double>                                ";
            message+=structwt;
            message+="\nPThreshold   <double>                                ";
            message+=pthreshold;
			message+="\nSubsample    <double>                                ";
			message+=chosensubsampling;
            message+="\nMoveFreq     (number number number number number [number])    (";
            for (int i=0;i<6;i++) {
                message+=movefreqvector[i];
                message+=" ";
            }
            message+=")";
            message+="\nFile         <file-name>                             ";
            message+=treefilename;
            //message+="\nRecordSearch No|Yes                                  No";
            //message+="\nSearchFile   <file-name>                             ";
            //message+="\nStatus       No|Yes                                  Yes";
           // message+="\nSteepest     No|Yes                                  No";
			//message+="\nExhaustive   No|Yes                                  No";
            message+="\nShowTries    No|Yes                                  ";
            if (showtries) {
				message+="Yes";
			}
			else {
				message+="No";
			}			
			message+="\nCOAL:        No|Yes                                  ";
			if (useCOAL) {
				message+="Yes";
			}
			else {
				message+="No";
			}
			message+="\nMS:          No|Yes                                  ";
			if (useMS) {
				message+="Yes";
			}
			else {
				message+="No";
			}
			
			message+="\nAIC_mode      0|1|2|3|4                              ";
			message+=COALaicmode;
			message+="\nGridWidth     <double>                               ";
			message+=contourstartingwidth;
			message+="\nGridSize      <integer-value>                        ";
			message+=contourstartingnumbersteps;
			message+="\nMaxRecursions <integer-value>                        ";
			message+=contourMaxRecursions;
			message+="\nBranch_export 0|1|2                                  ";
			message+=contourBrlenToExport;
			
            message+="\n\nNReps: Number of random starting species trees to use";
            //message+="\nTimeLimit: Limit search to X seconds";
            //message+="\nClock: Count seconds for time limit using actual elapsed time ('Wall'->Clock on a wall) or CPU time"; //NOte to self: see discussion online of time() fn and clock() fn in C++
            message+="\nRearrLimit: Limit search to X rearrangements for each nrep";
            //message+="\nAssignFixed: Assignment of gene samples to species is not optimized during a search";
            //message+="\nSppNumFixed: The total number of species is not optimized during the search (if set to No, then AssignFixed is also set to No)";
            message+="\nMaxNumSpp: The maximum number of species to split the samples into (only relevant if SppNumFixed==No)";
			message+="\nMinSamp: The minimum  number of samples per species (NOT the number of species)";
            message+="\nMoveFreq: Sets the relative proportion of times to try\n\t1) Species tree branch swaps\n\t2) Moving samples from one species to another\n\t3) Increasing the number of species\n\t4) Decreasing the number of species\n\t5) Attempt to reroot the species tree\n  If these don't sum to one, the program will automatically correct this.";
            //message+="\nRecordSearch: Output each step to a file: assignments go to a .txt file, and trees go to a .tre file";
            //message+="\nSearchFile: If RecordSearch==Yes, the prefix to use for the output files.";
           // message+="\nStatus: Output status of search to screen (and a log file, if open).";
           // message+="\nSteepest: Whether to look at all rearrangements and then take the best one or just take the first better one.";
			message+="\nSubsample: How extensively to try taxon reassignments on leaf splits. \n\tA value of 1 means try all of the possible reassignments, \n\ta value of 2 means try the square root of all the possible assignments,\n\t3 means the cube root, etc. A higher number means a faster but less effective search.\n\tThe program won't let you try fewer than 10 assignments on average.";
			message+="\nMS: Use Hudson's program MS to simulate gene trees to estimate probabilities. MAKE SURE MS IS INSTALLED SOMEWHERE ACCESSIBLE TO BROWNIE.";
			message+="\nCOAL: Use Degnan's program COAL to optimize the likelihood of the species delimitation and tree rather than the semiparametric penalty function. MAKE SURE COAL IS INSTALLED SOMEWHERE ACCESSIBLE TO BROWNIE.";
			message+="\nAIC_mode: When using COAL, use the 0: likelihood as the penalty term, 1: AIC value (k=number of species), 2: AICc with n=number of genes, 3: AICc with n=number of samples, 4: AICc with n=(number of genes) * (number of samples)";
			message+="\nGridWidth, GridSize, MaxRecursions all affect grid search";
			message+="\nBranch_export: if set to 0, returns a single estimate of the best branch lengths on the species tree. If set to 1, returns a table of equally-good branch lengths. If set to 2, returns a table of the best and neighboring branch lengths, suitable for doing a contour plot of score versus branch lengths";
            PrintMessage();
            finishexecuting=false;
        }
        else if(token.Abbreviation("NReps")) {
            nxsstring numbernexus;
            numbernexus = GetNumber(token);
            nreps=atoi( numbernexus.c_str() ); //convert to int
        }
        else if(token.Abbreviation("MAXNumspp")) {
            nxsstring numbernexus;
            numbernexus = GetNumber(token);
            maxnumspecies=atoi( numbernexus.c_str() ); //convert to int
            message="There must be no more than ";
			message+=maxnumspecies;
			message+=" species total.";
			PrintMessage();

        }
		else if(token.Abbreviation("MINSamp")) {
            nxsstring numbernexus;
            numbernexus = GetNumber(token);
            minsamplesperspecies=GSL_MIN(atoi( numbernexus.c_str() ),taxa->GetNumTaxonLabels()); //convert to int
			message="All species must now have at least ";
			message+=minsamplesperspecies;
			message+=" samples in them.";
			PrintMessage();
        }
		else if(token.Abbreviation("REArrlimit")) {
            nxsstring numbernexus;
			message="Rearrangement limit set to ";
            numbernexus = GetFileName(token);
			if (numbernexus[0]=='n' || numbernexus[0]=='N') {
					rearrlimit=-1;
				message+="None";
			}
			else {
				rearrlimit=atoi(numbernexus.c_str() );
				message+=rearrlimit;
			}
			PrintMessage();
        }		
        else if( token.Abbreviation("STEepest") ) {
            nxsstring yesnoreplace=GetFileName(token);
            if (yesnoreplace[0] == 'n' || yesnoreplace[0] == 'N') {
                steepest=false;
            }
            else {
                steepest=true;
                cout<<"%%%%%%%%%%%%%% Warning: steepest hasn't been tested everywhere %%%%%%%%%%%%%%%%%%%"<<endl;
            }
        }
		else if( token.Abbreviation("EXhaustive") ) {
            nxsstring yesnoreplace=GetFileName(token);
            if (yesnoreplace[0] == 'n' || yesnoreplace[0] == 'N') {
                exhaustive=false;
            }
            else {
                exhaustive=true;
                message="Wow, an exhaustive search.";
				PrintMessage();
            }
        }
        else if( token.Abbreviation("SHOwtries") ) {
            nxsstring yesnoreplace=GetFileName(token);
            if (yesnoreplace[0] == 'n' || yesnoreplace[0] == 'N') {
                showtries=false;
            }
            else {
                showtries=true;
            }
        }
		else if( token.Abbreviation("COAL") ) {
            nxsstring yesnoreplace=GetFileName(token);
            if (yesnoreplace[0] == 'n' || yesnoreplace[0] == 'N') {
                useCOAL=false;
            }
            else {
                useCOAL=true;
				useMS=false;
				message="This will use COAL for the search";
				PrintMessage();
				
            }
        }
		else if( token.Abbreviation("MS") ) {
            nxsstring yesnoreplace=GetFileName(token);
            if (yesnoreplace[0] == 'n' || yesnoreplace[0] == 'N') {
                useMS=false;
            }
            else {
                useMS=true;
				useCOAL=false;
				message="This will use ms for the search";
				PrintMessage();
				if (movefreqvector[5]==0) {
					movefreqvector[0]=0.1;
					movefreqvector[1]=0.1;
					movefreqvector[2]=0.1;
					movefreqvector[3]=0.1;
					movefreqvector[4]=0.1;
					movefreqvector[5]=0.5;
					message="Now allocating 50% of moves to branch length changes";
					PrintMessage();
				}
				
            }
        }
		
		else if (token.Abbreviation("AIC_mode")) {
			nxsstring numbernexus;
            numbernexus = GetNumber(token);
            int AICmoderaw=atoi( numbernexus.c_str() ); //convert to int
			if (AICmoderaw==0) {
				message="Using likelihood as penalty function";
				COALaicmode=AICmoderaw;
				PrintMessage();
			}
			else if (AICmoderaw==1) {
				message="Using AIC as penalty function";
				COALaicmode=AICmoderaw;
				PrintMessage();
			}
			else if (AICmoderaw==2) {
				message="Using AICc as penalty function, with n=number of genes";
				COALaicmode=AICmoderaw;
				PrintMessage();
			}
			else if (AICmoderaw==3) {
				message="Using AICc as penalty function, with n=number of samples";
				COALaicmode=AICmoderaw;
				PrintMessage();
			}
			else if (AICmoderaw==4) {
				message="Using AICc as penalty function, with n=(number of genes) * (number of samples)";
				COALaicmode=AICmoderaw;
				PrintMessage();
			}
			else {
				errormsg="You entered an unrecognized option for AIC_mode: should be an integer 0-4, but you entered ";
				errormsg+=AICmoderaw;
				throw XNexus( errormsg);
			}
			
		}
		else if (token.Abbreviation("BRANch_export")) {
			nxsstring numbernexus;
            numbernexus = GetNumber(token);
            int contourBrlenToExportraw=atoi( numbernexus.c_str() ); //convert to int
			if (contourBrlenToExportraw==0) {
				message="Just export a single branch length for each best species topology";
				contourBrlenToExport=contourBrlenToExportraw;
				PrintMessage();
			}
			else if (contourBrlenToExportraw==1) {
				message="Export a table of equally-good branch lengths for each best species topology";
				contourBrlenToExport=contourBrlenToExportraw;
				PrintMessage();
			}
			else if (contourBrlenToExportraw==2) {
				message="Export a table from a grid of possible branch lengths, including suboptimal ones (good for plotting contours)";
				contourBrlenToExport=contourBrlenToExportraw;
				PrintMessage();
			}
			
			else {
				errormsg="You entered an unrecognized option for Branch_export: should be an integer 0-2, but you entered ";
				errormsg+=contourBrlenToExportraw;
				throw XNexus( errormsg);
			}
		}
		else if(token.Abbreviation("MAXRecursions")) {
            nxsstring numbernexus;
			 numbernexus = GetFileName(token);
			int contourMaxRecursionsraw=atoi(numbernexus.c_str() );
			if (contourMaxRecursionsraw>0) {
				message="Max recursions per grid search set to ";	
				contourMaxRecursions=contourMaxRecursionsraw;
				message+=contourMaxRecursions;
				PrintMessage();
			}
			else {
				errormsg="You entered an inappropriate option for MaxRecursions: should be an integer >0, but you entered ";
				errormsg+=contourMaxRecursionsraw;
				throw XNexus( errormsg);
			}
        }		
		else if(token.Abbreviation("GRIDSize")) {
            nxsstring numbernexus;
			numbernexus = GetFileName(token);
			int contourstartingnumberstepsraw=atoi(numbernexus.c_str() );
			if (contourstartingnumberstepsraw>0) {
				message="Number of steps on each axis in grid search set to ";	
				if (GSL_IS_ODD (contourstartingnumberstepsraw)==0) {
					contourstartingnumberstepsraw++;
				}
				contourstartingnumbersteps=1.0*contourstartingnumberstepsraw; //is actually a double
				message+=contourstartingnumbersteps;
				message+=" even numbers increased to next odd number (so the current best value is the center of the grid)";
				PrintMessage();
			}
			else {
				errormsg="You entered an inappropriate option for GridSize: should be an integer >0, but you entered ";
				errormsg+=contourstartingnumberstepsraw;
				throw XNexus( errormsg);
			}
        }		
		else if(token.Abbreviation("GRIDWidth")) {
            nxsstring numbernexus;
			numbernexus = GetFileName(token);
			double contourstartingwidthraw=atof(numbernexus.c_str() );
			if (contourstartingwidthraw>0) {
				message="Starting maximum branch length for each branch in grid search set to ";	
				contourstartingwidth=contourstartingwidthraw;
				message+=contourstartingwidth;
				message+=" times the best length";
				PrintMessage();
			}
			else {
				errormsg="You entered an inappropriate option for GridWidth: should be a number >0, but you entered ";
				errormsg+=contourstartingwidthraw;
				throw XNexus( errormsg);
			}
        }	
        else if( token.Abbreviation("MoveFreq") ) {
            token.GetNextToken();
            token.GetNextToken(); //eat the equals sign
            vector<double> temporarymovefreqvector;
            if (!token.Equals("(")) {
                errormsg="Expecting next token to be a \'(\' but instead got ";
                errormsg+=token.GetToken();
                throw XNexus( errormsg);
            }
            int inputcount=0;
            while (!token.Equals(")")) {
                nxsstring numbernexus;
                numbernexus=GetNumberOnly(token);
                if (numbernexus!=")") {
                    temporarymovefreqvector.push_back(atof( numbernexus.c_str() ));
                    inputcount++;
                }
                else {
                    break;
                }
            }
            if (inputcount!=5 && inputcount!=6) {
                errormsg="You should have entered five (or six) frequencies, you entered ";
                errormsg+=inputcount;
                throw XNexus( errormsg);
            }
            else {
                double sumoffreqs=temporarymovefreqvector[0]+temporarymovefreqvector[1]+temporarymovefreqvector[2]+temporarymovefreqvector[3]+temporarymovefreqvector[4];
				if (inputcount==6) {
					sumoffreqs+=temporarymovefreqvector[5];
				}
				else {
					temporarymovefreqvector.push_back(0);
				}
                movefreqvector.clear();
                for (int i=0; i<temporarymovefreqvector.size(); i++) {
                    movefreqvector[i]=(temporarymovefreqvector[i])/sumoffreqs;
                }

            }
        }

        else if(token.Abbreviation("MINNumspp") || token.Abbreviation("MINUMSP") ) { //to deal with typos
            nxsstring numbernexus;
            numbernexus = GetNumber(token);
            minnumspecies=atoi( numbernexus.c_str() ); //convert to int
            message="There must be at least ";
			message+=minnumspecies;
			message+=" species total.";
			PrintMessage();

        }
        else if(token.Abbreviation("FIlename")) {
            treefilename=GetFileName(token);
        }
        else if(token.Abbreviation("STRuctwt")) {
            nxsstring numbernexus;
            numbernexus = GetNumber(token);
            double proposedwt=atof( numbernexus.c_str() ); //convert to double
            if ((proposedwt<0) || (proposedwt>1)) {
                errormsg="Error: The weight for population substructure must be between 0 and 1 (inclusive)";
                throw XNexus( errormsg);
            }
            else {
                structwt=proposedwt;
                message="You have chosen to have the gene duplication cost to make up ";
                double gtpwt=(1-structwt)*100;
                message+=gtpwt;
                message+="% and the substructure cost to make up the remaining ";
                message+=100-gtpwt;
                message+="%.";
                PrintMessage();
            }
        }
        else if(token.Abbreviation("PThreshold")) {
            nxsstring numbernexus;
            numbernexus = GetNumber(token);
            double proposedthresh=atof( numbernexus.c_str() ); //convert to double
            if ((proposedthresh<0) || (proposedthresh>1)) {
                errormsg="Error: The p-value threshold must be between 0 and 1 (inclusive)";
                throw XNexus( errormsg);
            }
            else {
                pthreshold=proposedthresh;
                message="You have chosen to have structure p values > ";
                message+=pthreshold;
                message+=" ignored";
                PrintMessage();
            }
        }
		else if(token.Abbreviation("SUBsample")) {
			nxsstring numbernexus;
            numbernexus = GetNumber(token);
            double proposedsubsample=atof( numbernexus.c_str() ); //convert to double
			if (proposedsubsample<1) {
				errormsg="Error: The subsample threshold must be greater than one";
                throw XNexus( errormsg);
			}
			else {
				chosensubsampling=proposedsubsample;
				message="You have chosen a split leaf subsample value of ";
				message+=chosensubsampling;
				PrintMessage();
			}
			
			
		}
    }
}

vector<double> BROWNIE::GetCombinedScore(ContainingTree *SpeciesTreePtr)
{
	vector<double> scorevector;
	bool calculatescore=true;
	bool infinitescore=false;
	if (useCOAL) {
		int maxspecies=0;
		int currentnumberofspecies=0;
		float neglnlikelihood=GSL_POSINF;
		
		for (int i=0;i<convertsamplestospecies.size();i++) {
			maxspecies=GSL_MAX(maxspecies,convertsamplestospecies[i]);
		}
		vector<int> speciesvector(maxspecies,0);
		for (int i=0;i<convertsamplestospecies.size();i++) {
			speciesvector[(convertsamplestospecies[i])-1]++;
		}
		for (int i=0;i<speciesvector.size();i++) {
			if (speciesvector[i]<minsamplesperspecies) {
				calculatescore=false;
			}
		}		
		if (calculatescore && maxspecies==1) { //COAL can't calculate probability of gene tree on species tree with one species, so use Harding 1971 equation 5.3 to do it manually
			neglnlikelihood=0;
			for (int chosentreenum=0; chosentreenum<trees->GetNumTrees(); chosentreenum++) { //loop over all the gene trees
				Tree CurrentGeneTreeTreeFmt=intrees.GetIthTree(chosentreenum);
				CurrentGeneTreeTreeFmt.Update(); //sets weights
				NodeIterator <Node> n (CurrentGeneTreeTreeFmt.GetRoot()); //starts at leaves and works down
				double probability=1;
				cur = n.begin();
				while (cur) {
					if (!(cur->IsLeaf())) {
						int Q=(cur->GetChild())->GetWeight();
						int R=((cur->GetChild())->GetSibling())->GetWeight();
						int N=R+Q;
						probability*=2*(gsl_sf_fact(R))*(gsl_sf_fact(Q))/((N-1.0)*gsl_sf_fact(N)); //in Harding's equation, we calculate [Pl{Q] and Pl[R] when we are at the child nodes
						//cout<<"probability now "<<probability<<" with Q="<<Q<<" and R="<<R<<endl;
					}
					//else {
					//	cout<<"Leaf "<<cur->GetLabel()<<endl;
					//}
					cur = n.next();
				}
				neglnlikelihood+=-1.0*log(probability);
				//cout<<"neglnlikelihood now "<<neglnlikelihood<<endl<<endl;
			}
		}
		else if (calculatescore && maxspecies>1) {
		//export the species tree
			nxsstring coalspecies="coalspecies.txt";

			ofstream coalspeciesstream;
			coalspeciesstream.open( coalspecies.c_str() );
			ContainingTree CurrentTree=*SpeciesTreePtr;
			CurrentTree.FindAndSetRoot();
			CurrentTree.Update();
			CurrentTree.GetNodeDepths();
		//	if(logf_open) {
		//		logf<<"Now trying with species tree: ";
		//		CurrentTree.Write(logf);
		//	}
		//CurrentTree.Draw(cout);
			NodeIterator <Node> n (CurrentTree.GetRoot());
			cur = n.begin();
			while (cur) {
				if (cur->IsLeaf()) {
					currentnumberofspecies++;
				//nxsstring newlabel=PipeLeafFinalSpeciesTree();
				//cur->SetLabel(newlabel);
				}
				else {
					cur->SetLabel("");
				}
				cur = n.next();
			}
			(CurrentTree).Write(coalspeciesstream);
			coalspeciesstream<<endl;
			coalspeciesstream.close();
			
		//export the gene trees
			nxsstring coalgenes="coalgenes.txt";
			ofstream coalgenesstream;
			coalgenesstream.open( coalgenes.c_str());
			for (int chosentreenum=0; chosentreenum<trees->GetNumTrees(); chosentreenum++) { //loop over all the gene trees
				Tree CurrentGeneTreeTreeFmt=intrees.GetIthTree(chosentreenum);
				NodeIterator <Node> n (CurrentGeneTreeTreeFmt.GetRoot());
				cur = n.begin();
				while (cur) {
					if (cur->IsLeaf()) {
						int SampleNumber=taxa->FindTaxon(cur->GetLabel());
						nxsstring newlabel="taxon";
						newlabel+=convertsamplestospecies[SampleNumber];
						newlabel+="-";
						newlabel+=SampleNumber;
						cur->SetLabel(newlabel);
					//cout<<"current label is "<<cur->GetLabel()<<" perhaps with quotes"<<endl;
					}
					else {
						cur->SetLabel("");
					}
					cur = n.next();
				}
				(CurrentGeneTreeTreeFmt).WriteNoQuote(coalgenesstream);
				coalgenesstream<<endl;
				
				
			}
			coalgenesstream.close();
			
		//export the infile
			nxsstring coalinfile="infile";
			ofstream coalinfilestream;
			coalinfilestream.open(coalinfile.c_str() );
			coalinfilestream<<"begin coal;"<<endl;
			coalinfilestream<<"nstaxa ";
			coalinfilestream<<currentnumberofspecies<<";"<<endl;
			coalinfilestream<<"ngtaxa";
			nxsstring speciesnames="taxa names = ";
			for (int i=1; i<=currentnumberofspecies; i++) {
				int samplecount=0;
				for (int j=0; j<=convertsamplestospecies.size();j++) {
					if (convertsamplestospecies[j]==i)  {
						samplecount++;
					}
				}
				speciesnames+=" taxon";
				speciesnames+=i;
				coalinfilestream<<" "<<samplecount;
			}
			coalinfilestream<<";"<<endl;
			coalinfilestream<<speciesnames<<";"<<endl;
			coalinfilestream<<"intra yes;"<<endl;
			coalinfilestream<<"gene tree file=coalgenes.txt;"<<endl;
			coalinfilestream<<"species tree file=coalspecies.txt;"<<endl;
			coalinfilestream<<"nstrees 1;"<<endl;
			coalinfilestream<<"ngtrees "<<trees->GetNumTrees()<<";"<<endl;
			nxsstring coaloutputfile="coalout.txt";
			coalinfilestream<<"outfile = "<<coaloutputfile<<"  / probs   ;"<<endl;		
			coalinfilestream<<"end;"<<endl;
			coalinfilestream.close();
			
		//run coal
		//cout<<"Starting to call coal"<<endl;
			int coalreturncode=-1;
			if (currentnumberofspecies>0) { //just to head off a weird error
				int returncodefailurecount=-1;
				while (coalreturncode!=0) {
					returncodefailurecount++;
					coalreturncode=system("coal > /dev/null");
					if (returncodefailurecount>100) {
						errormsg="Had over 100 failures of COAL";
						throw XNexus( errormsg);
						break;
					}
				}
			}
		//cout<<"Done calling  coal, now parsing"<<endl;
		//cout <<"coal returncode is "<<coalreturncode<<endl;
			if (coalreturncode==0) {
				int returncount=0;
				neglnlikelihood=0;
				ifstream coalin;
				coalin.open( coaloutputfile.c_str(), ios::binary | ios::in );
				if (coalin) {
					char inputitem [COMMAND_MAXLEN];
					coalin>>inputitem;
					bool nextisprob=true; //since the input file alternates gene_tree_number SPACES prob
					while (!coalin.eof()) {
						coalin>>inputitem;
						if (nextisprob) {
							double probability=atof(inputitem);
							if (probability==0 || probability !=probability) { //test for it not being a number or for it being nan
								neglnlikelihood=GSL_POSINF;
								break;
							}
							neglnlikelihood+=-1.0*log(probability);
							returncount++;
						//cout<<"prob is "<<inputitem<<endl;
							nextisprob=false;
						}
						else {
							nextisprob=true;
						}
					}
				}
				coalin.close();
				if (returncount!=trees->GetNumTrees()) {
					neglnlikelihood=GSL_POSINF; //THere's an error, we returned the wrong number of likelihood scores
				}
				
			}
		//cout<<"Done parsing"<<endl<<endl;
		//process the output
		}
		assert(neglnlikelihood>0);
		double score=neglnlikelihood;
		if (COALaicmode==0) {
			//do nothing, we're just using raw lnL
		}
		else if (COALaicmode==1) {
			score=AIC(neglnlikelihood,currentnumberofspecies);
		}
		else if (COALaicmode==2) {
			score=AICc(neglnlikelihood,currentnumberofspecies,trees->GetNumTrees());
		}
		else if (COALaicmode==3) {
			score=AICc(neglnlikelihood,currentnumberofspecies,taxa->GetNumTaxonLabels());
		}
		else if (COALaicmode==4) {
			score=AICc(neglnlikelihood,currentnumberofspecies,(trees->GetNumTrees())*(taxa->GetNumTaxonLabels()));
		}
		scorevector.push_back(neglnlikelihood);
		scorevector.push_back(0);
		scorevector.push_back(0);
		//cout<<"neglnlikelihood = "<<neglnlikelihood<<endl;
		//if(logf_open) {
		//	logf<<"\t[neglnlikelihood = "<<neglnlikelihood<<" ]"<<endl;
		//}
		
	}
	else if (useMS) {
		double neglnlikelihood=0.0;
		//ms $nsamples $ntrees -T -I 2 $nsamplessp1 $nsamplessp2 -ej $splittime 1 2
		//echo "cat" | grep "\(cat\|dog\)"
		//Convert species tree to MS format
		ContainingTree CurrentTree=*SpeciesTreePtr;
		CurrentTree.FindAndSetRoot();
		CurrentTree.Update();
		CurrentTree.GetNodeDepths();
		int numspecies=CurrentTree.GetNumLeaves();		
		int ntax=taxa->GetNumTaxonLabels();
		nxsstring msstring="ms ";
		msstring+=ntax;
		msstring+=" ";
		msstring+=msbasereps; //number of trees to generate
		msstring+=" -T -I ";
		msstring+=numspecies;
		msstring+=" ";
		vector <int> samplesperspecies;
		vector <nxsstring> labelswithinspecies;
		double numberofpermutations=1;
		int previoustotal=0;
		for (int currentspecies=1; currentspecies<=numspecies; currentspecies++) {
			int samplecountthisspecies=0;
			nxsstring labelregex="\\(";
			for (int currentsample=0; currentsample<ntax; currentsample++) {
				if (convertsamplestospecies[currentsample]==currentspecies) {
					samplecountthisspecies++;
				}
			}
			msstring+=samplecountthisspecies;
			samplesperspecies.push_back(samplecountthisspecies);
			for (int i=1; i<=samplecountthisspecies; i++) {
				if (i>1) {
					labelregex+="\\|";
				}
				labelregex+=i+previoustotal;
				
			}
			labelregex+="\\):[0-9]*.[0-9]*"; //ms exports brlen, too
			labelswithinspecies.push_back(labelregex);
			previoustotal+=samplecountthisspecies;
			numberofpermutations*=gsl_sf_fact(samplecountthisspecies);
			msstring+=" ";
		}
		NodeIterator <Node> n (CurrentTree.GetRoot());
		//loop once to get MS string
		cur = n.begin();
		while (cur) {
			if (cur->IsLeaf()) {
				int speciesnumber;
				nxsstring specieslabel=cur->GetLabel();
				string speciesstring=specieslabel.c_str();
				size_t index = speciesstring.find('t');
				speciesstring.erase(index,5); //erase "taxon"
				cur->SetLabel(speciesstring);
			}
			else {
				//assumes binary, ultrametric tree
				if (cur->GetChild()!=NULL) { //will only really occur if the tree has just a root node
					if ((cur->GetChild())->GetSibling()!=NULL) {
						int childid=atoi(((cur->GetChild())->GetLabel()).c_str());
						int sibid=atoi((((cur->GetChild())->GetSibling())->GetLabel()).c_str());
						nxsstring newlabel="";
						newlabel+=GSL_MIN(childid,sibid);
						cur->SetLabel(newlabel);
						double splittime=0;
						NodePtr childnode=cur->GetChild();
						while (childnode!=NULL) {
							splittime+=childnode->GetEdgeLength();
							childnode=childnode->GetChild();
						}
						msstring+="-ej ";
						msstring+=splittime;
						msstring+=" ";
						msstring+=GSL_MAX(childid,sibid);
						msstring+=" ";
						msstring+=GSL_MIN(childid,sibid);
						msstring+=" ";
					}
				}
			}
			cur = n.next();
		}
		
		//rather than looking for exact match, get probabilities of given topology with all possible permutations of labels of samples from a given species onto ms taxon numbers, then divide by number of such permutations
		
		
		//loop over all the gene trees
		for (int chosentreenum=0; chosentreenum<trees->GetNumTrees(); chosentreenum++) { //loop over all the gene trees
			nxsstring grepstring="";
			int numberintraspecificcherries=0;
			Tree CurrentGeneTreeTreeFmt=intrees.GetIthTree(chosentreenum);
			NodeIterator <Node> n (CurrentGeneTreeTreeFmt.GetRoot());
			//traverse once to count intraspecific sister cherries
			cur = n.begin();
			while (cur) {
				if (cur->IsLeaf()) {
					NodePtr sister=cur->GetSibling();
					if (sister!=NULL) {
						if (sister->IsLeaf()) {
							//cout<<cur->GetLabel()<<" with "<<sister->GetLabel()<<endl;
							if (convertsamplestospecies[taxa->FindTaxon(cur->GetLabel())]==convertsamplestospecies[taxa->FindTaxon(sister->GetLabel())]) { //they have the same species (though node labels are changed, we do the one that's a child before its sib)
								numberintraspecificcherries++;
							}
						}
					}
				}
				cur = n.next();
			}
			
			//now get the grep lines
			cur = n.begin();
			while (cur) {
				if (cur->IsLeaf()) {
					cur->SetLabel(labelswithinspecies[-1+convertsamplestospecies[taxa->FindTaxon(cur->GetLabel())]]); //gets regex for tip labels
				}
				else {
					nxsstring combinedstring="";
					nxsstring leftchild=(cur->GetChild())->GetLabel();
					nxsstring rightchild=((cur->GetChild())->GetSibling())->GetLabel();
					nxsstring additionalregex="(\\(";
					additionalregex+=leftchild;
					additionalregex+=",";
					additionalregex+=rightchild;
					additionalregex+="\\|";
					additionalregex+=rightchild;
					additionalregex+=",";
					additionalregex+=leftchild;
					additionalregex+="\\))";
					if (cur->GetAnc()!=NULL) { // not root
						additionalregex+=":*[0-9]*.[0-9]*";
					}
					else { //originally, did this at each node, with the idea that from MS, you first filter out the lines that don't match a cherry, then filter from that filtered set, etc., thinking it would be faster than just using the regex at the root. Actually, just using the regex at the root was faster, so I do that.
						grepstring+=" | grep -c \"";
						grepstring+=additionalregex;
						grepstring+="\"";
					}
					cur->SetLabel(additionalregex);
				}
				cur = n.next();
			}
			nxsstring finalsystemcall=msstring;
			//cout<<msstring<<endl;
			finalsystemcall+=grepstring;
			nxsstring msinputfile="mscount.txt";
			finalsystemcall+=" > ";
			finalsystemcall+=msinputfile;
			//system("rm mscount.txt");
			int returncode=system(finalsystemcall.c_str());
			if (debugmode) {
				cout<<"return code is "<<returncode<<" (divided by 256, is "<<returncode/256<<")"<<endl;
			}
			//cout<<finalsystemcall<<endl;
		/*	int returnattempts=-1;
			int returncode=-1;
			while (returncode!=0) {
				if (returnattempts>0) {
					system("rm mscount.txt");
				}
				returncode=system(finalsystemcall.c_str());
				returnattempts++;
				if (returnattempts>5) {
					message="Had over 5 failures of MS, calling\n\n";
					message+=finalsystemcall;
					message+="\n";
					PrintMessage();
					//throw XNexus( errormsg);
					break;
				}
				if (returncode!=0) {
					message="Warning: ms returned error code ";
					message+=returncode;
					PrintMessage();
				}
				
			}*/
			//	if (returncode==0) {
			ifstream msin;
			msin.open( msinputfile.c_str(), ios::binary | ios::in );
			if (msin) {
				char inputitem [COMMAND_MAXLEN];
				msin>>inputitem;
				double numbermatches=atof(inputitem);
				if (numbermatches==0) {
					infinitescore=true;
				}
				neglnlikelihood+=-1.0*(log(GSL_MAX(numbermatches,0.01))-log(msbasereps)-log((1.0*numberofpermutations)/(1.0*numberintraspecificcherries))); //numbermatches is an integer, but log(0) is infinite. Idea is to make this a really big but not infinite number. Probability of a gene tree is #of times it was observed (with all permutations of tip labels within species) divided by the number of trees returned, all divided by the number of possible permutations (since this gene tree is just one realization of that), correcting for the number of intraspecific cherries (otherwise, for example, with a single species with a gene with two samples, (1,2) might be the gene tree, but we say the number of perms=2, and so though we find ((1|2),(1|2)) in all the ms returns, we would divide this by 2, when the real probability is one.
				
				//		}
				msin.close();
			}
			else {
				message="\nWarning: got return code of ";
				message+=returncode;
				message+=" for ms string of \n\n\t";
				message+=msstring;
				message+="\n\nand grep string of \n\n";
				message+=grepstring;
				message+="\n\n";
			
				neglnlikelihood=GSL_POSINF;
			}
			
		}
		double score=neglnlikelihood;
		if (COALaicmode==0) {
			//do nothing, we're just using raw lnL
		}
		else if (COALaicmode==1) {
			score=AIC(neglnlikelihood,numspecies);
		}
		else if (COALaicmode==2) {
			score=AICc(neglnlikelihood,numspecies,trees->GetNumTrees());
		}
		else if (COALaicmode==3) {
			score=AICc(neglnlikelihood,numspecies,taxa->GetNumTaxonLabels());
		}
		else if (COALaicmode==4) {
			score=AICc(neglnlikelihood,numspecies,(trees->GetNumTrees())*(taxa->GetNumTaxonLabels()));
		}
		scorevector.push_back(neglnlikelihood);
		scorevector.push_back(0);
		scorevector.push_back(0);
	}
	else {
    //We do this so we don't bother doing the GTP or triplet calculations if they're not needed.
		triplettoohigh=false;
		gtptoohigh=false;
		double combinedscore=0;
		double tripletscore=0;
		double gtpscore=0;
		if (minsamplesperspecies>1) {
			int maxspecies=0;
			for (int i=0;i<convertsamplestospecies.size();i++) {
				maxspecies=GSL_MAX(maxspecies,convertsamplestospecies[i]);
			}
			vector<int> speciesvector(maxspecies,0);
			for (int i=0;i<convertsamplestospecies.size();i++) {
				speciesvector[(convertsamplestospecies[i])-1]++;
			}
			for (int i=0;i<speciesvector.size();i++) {
				if (speciesvector[i]<minsamplesperspecies) {
					calculatescore=false;
					combinedscore=GSL_POSINF;
					tripletscore=GSL_POSINF;
					gtpscore=GSL_POSINF;
				}
			}		
		}
		if (calculatescore) {
			if (structwt>0) {
				tripletscore=double(GetTripletScore(SpeciesTreePtr));
        //	(*SpeciesTreePtr).Write(cout);
        //	cout<<endl;
        //	for (int i=0;i<convertsamplestospecies.size();i++) {
        //			cout<<convertsamplestospecies[i]<<" ";
        //	}
        //	cout<<endl;
        //	cout<<"\ntriplet score = "<<newscore<<endl<<endl;
        //  combinedscore+=newscore;
				
			}		
			if (structwt<1 && (triplettoohigh==false)) {
        //combinedscore+=(1.0-structwt)*GetGTPScore(SpeciesTreePtr);
        // cout<<"Mike gave "<<combinedscore/(1.0-structwt)<<"\nI gave "<<GetGTPScoreNew(SpeciesTreePtr)<<endl;
        //cout<<"SpeciesTree"<<endl;
        //(*SpeciesTreePtr).Write(cout);
        //cout<<endl;
				gtpscore=double(GetGTPScoreNew(SpeciesTreePtr));
			}
			combinedscore=((1.0-structwt)*gtpscore)+(structwt*tripletscore);
		}
		scorevector.push_back(combinedscore);
		scorevector.push_back(gtpscore);
		scorevector.push_back(tripletscore);
	}
    return scorevector;
}

//checks to make sure that the convertsamplestospecies vector is not missing assignments to any species (i.e, doesn't consist only of species 1, 3 and 4). If fix==true, it'll fix this. It will return true if check
bool BROWNIE::CheckConvertSamplesToSpeciesVector(bool fix)
{
    int maxspecies=0;
    bool goodshape=true;
    int firstbadspecies=0;
    for (int i=0;i<convertsamplestospecies.size();i++) {
        maxspecies=GSL_MAX(maxspecies,convertsamplestospecies[i]);
    }
    vector<int> speciesvector(maxspecies,0);
    for (int i=0;i<convertsamplestospecies.size();i++) {
        speciesvector[(convertsamplestospecies[i])-1]++;
    }
    for (int i=0;i<speciesvector.size();i++) {
        if (speciesvector[i]==0) {
            goodshape=false;
            firstbadspecies=i+1;
            break;
        }
    }
    if (!goodshape && fix) {
        for (int i=0;i<convertsamplestospecies.size();i++) {
            if (convertsamplestospecies[i]>firstbadspecies) {
                convertsamplestospecies[i]=convertsamplestospecies[i]-1;
            }
        }
        goodshape=CheckConvertSamplesToSpeciesVector(true);

    }
    return goodshape;
}

bool BROWNIE::CombineSpeciesWithTooFewSamples(bool fix) { //The idea of this is to fix the assignment vector if there are too few samples per species by combining the too small species with the next species and then fixing the resulting hole
	bool goodsamplenumber=true;
	if (minsamplesperspecies>1) {
		int maxspecies=0;
		for (int i=0;i<convertsamplestospecies.size();i++) {
			maxspecies=GSL_MAX(maxspecies,convertsamplestospecies[i]);
		}
		vector<int> speciesvector(maxspecies,0);
		for (int i=0;i<convertsamplestospecies.size();i++) {
			speciesvector[(convertsamplestospecies[i])-1]++;
		}
		for (int i=0;i<speciesvector.size();i++) {
			if (speciesvector[i]<minsamplesperspecies) {
				goodsamplenumber=false;
				if (fix) {
					for (int j=0;j<convertsamplestospecies.size();j++) {
						if (convertsamplestospecies[j]==i+1) {
							if ((i+1)<maxspecies) { //so there's a species to the right to combine with.
								convertsamplestospecies[j]=convertsamplestospecies[j]+1;
							}
							else { 
								convertsamplestospecies[j]=convertsamplestospecies[j]-1;
							}
						}
					}
					bool goodshape=CheckConvertSamplesToSpeciesVector(true);
					goodsamplenumber=CombineSpeciesWithTooFewSamples(true);
				}
			}
		}		
	}
	return goodsamplenumber;
}

bool BROWNIE::CheckConvertSamplesToSpeciesVectorSpNum(int actualmaxspeciesnum)
{
    int maxspecies=0;
    bool goodshape=true;
    int firstbadspecies=0;
    for (int i=0;i<convertsamplestospecies.size();i++) {
        maxspecies=GSL_MAX(maxspecies,convertsamplestospecies[i]);
    }
	if (maxspecies!=actualmaxspeciesnum) {
		goodshape=false;
	}
	if (goodshape) {
		vector<int> speciesvector(maxspecies,0);
		for (int i=0;i<convertsamplestospecies.size();i++) {
			speciesvector[(convertsamplestospecies[i])-1]++;
		}
		for (int i=0;i<speciesvector.size();i++) {
			if (speciesvector[i]==0) {
				goodshape=false;
				firstbadspecies=i+1;
				break;
			}
		}
	}
	return goodshape;
}

bool BROWNIE::CheckConvertSamplesToSpeciesTooManySpecies()
{
    int observedmaxspecies=0;
    bool goodshape=true;
    for (int i=0;i<convertsamplestospecies.size();i++) {
        observedmaxspecies=GSL_MAX(observedmaxspecies,convertsamplestospecies[i]);
    }
	if (observedmaxspecies>maxnumspecies) {
		goodshape=false;
	}
	return goodshape;
}

bool BROWNIE::CheckConvertSamplesToSpeciesTooFewSpecies()
{
    int observedmaxspecies=0;
    bool goodshape=true;
    for (int i=0;i<convertsamplestospecies.size();i++) {
        observedmaxspecies=GSL_MAX(observedmaxspecies,convertsamplestospecies[i]);
    }
	if (observedmaxspecies<minnumspecies) {
		goodshape=false;
	}
	return goodshape;
}

bool BROWNIE::MoveSamples(vector<int> Originalconvertsamplestospecies)
{
    bool morereassignments=true;
    int sampletomove=SamplesToMove.back();
    int destination=SampleDestinations.back();
    for(int i=0;i<(CladeVector[sampletomove]).size();i++) {
        convertsamplestospecies[(CladeVector[sampletomove][i])]=destination;
    }
    SampleDestinations.pop_back();
    if (SampleDestinations.size()==0) { //we've finished trying to assign all destinations
        convertsamplestospecies.swap(Originalconvertsamplestospecies); //need to use the original vector; if this try results in a better vector, we'll throw out the Sample___ vectors anyway; if not, we'll go back to the Originalconvertsamplestospecies vector, so we have to make sure any moves are allowed under that vector
        int nspecies=0;
        for (int i=0;i<convertsamplestospecies.size();i++) {
            nspecies=GSL_MAX(convertsamplestospecies[i],nspecies);
        }
        vector<int> TempSampleDestinations;
        if (SamplesToMove.size()<2) { //we were already on the last sample, too
            morereassignments=false;
        }
        else {
            bool enoughdestinations=false;
            while (!enoughdestinations && SamplesToMove.size()>1) {
                SamplesToMove.pop_back(); //we're done trying to move that sample
                SampleDestinations.clear();
                for (int i=0;i<nspecies;i++) {
                    if(TestMoveSamples(SamplesToMove.back(),i+1)) {
                        TempSampleDestinations.push_back(i+1);
                        enoughdestinations=true;
                    }
                }
            }
            morereassignments=enoughdestinations;
            if(enoughdestinations) {
                gsl_permutation * v = gsl_permutation_alloc (TempSampleDestinations.size());
                gsl_permutation_init (v);
                gsl_ran_shuffle (r, v->data, TempSampleDestinations.size(), sizeof(size_t));
                for (int i=0; i<TempSampleDestinations.size(); i++) {
                    SampleDestinations.push_back(TempSampleDestinations[gsl_permutation_get (v,i)]);
                }
                gsl_permutation_free(v);
            }
        }
        convertsamplestospecies.swap(Originalconvertsamplestospecies);
    }
    return morereassignments;
}

bool BROWNIE::TestMoveSamples(int Sample, int Destination)
{
    bool validmove=true;
    int maxspeciesorig=0;
    for (int i=0; i<convertsamplestospecies.size();i++) {
        maxspeciesorig=GSL_MAX(maxspeciesorig,convertsamplestospecies[i]);
    }
    vector<int> origconvertsamplestospecies=convertsamplestospecies;
    for(int i=1;i<=(CladeVector[Sample]).size();i++) {
		//convertsamplestospecies[(CladeVector.at(Sample)).at(i)]=Destination;
		int convertsamplestospeciessize=convertsamplestospecies.size();
		int cladevectorsamplesize=(CladeVector[Sample]).size();
		int cladevectorsize=CladeVector.size();
		if(cladevectorsamplesize>convertsamplestospeciessize || cladevectorsamplesize<0) {
			nxsstring outputname=treefilename;
			outputname+=".CladeVector.txt";
			ofstream cladevector;
			cladevector.open( outputname.c_str(), ios::out | ios::app );
			for (int k=0;k<CladeVector.size();k++) {
				cladevector<<endl;
				for (int j=0;j<(CladeVector[k]).size();j++) {
					cladevector<<" "<<CladeVector[k][j];
				}
			}
			cladevector<<endl;
			cladevector<<endl<<"convertsamplestospeciessize = "<<convertsamplestospeciessize<<endl;
			cladevector<<"cladevectorsamplesize = "<<cladevectorsamplesize<<endl;
			cladevector<<"cladevectorsize = "<<cladevectorsize<<endl;
			cladevector<<"Sample = "<<Sample<<endl<<"Destination = "<<Destination<<endl<<"i = "<<i<<endl;
			cladevector.close();
		}
        convertsamplestospecies[(CladeVector[Sample][i-1])]=Destination; // so it doesn't run if CladeVector[Sample].size()=0
    }
    validmove=CheckConvertSamplesToSpeciesVector(false);
    if(validmove) {
        if(convertsamplestospecies==origconvertsamplestospecies) {
            validmove=false; //since this induces no change
        }
    }
    if(validmove) {
        int maxspeciesfinal=0;
        for (int i=0; i<convertsamplestospecies.size();i++) {
            maxspeciesfinal=GSL_MAX(maxspeciesfinal,convertsamplestospecies[i]);
        }
        if (maxspeciesfinal!=maxspeciesorig) {
            validmove=false; //need to have same number of species at start and end
        }
    }
    convertsamplestospecies.swap(origconvertsamplestospecies);
    return validmove;
}

void BROWNIE::DoExhaustiveSearch()
{
	message="This will do an exhaustive search for up to 6 species. This will take some time.";
	PrintMessage();
	bestscore=GSL_POSINF;
	vector<double> nextscorevector;
	double nextscore;
	message="Now starting with 1 species";
	PrintMessage();
	ContainingTree CurrentTree;
	CurrentTree.RandomTree(1);
	CurrentTree.ConvertTaxonNamesToRandomTaxonNumbers();
	convertsamplestospecies.clear();
	for (int i=0;i<taxa->GetNumTaxonLabels();i++) {
		convertsamplestospecies.push_back(1);
		cout<<" 1";
	}
	if (useCOAL || useMS) {
		CurrentTree.InitializeMissingBranchLengths();
	}
	nextscorevector=GetCombinedScore(&CurrentTree);
	nextscore=nextscorevector[0]; 
	bestscore=nextscore; //this must be the best score, as it's the first tree
	//TotalScores.push_back(nextscorevector[0]);
	//GTPScores.push_back(nextscorevector[1]);
	//StructScores.push_back(nextscorevector[2]);
	FormatAndStoreBestTree(&CurrentTree,nextscorevector);
	
	message="Score: ";
	message+=nextscore;
	PrintMessage();
	message="Now starting with 2 species";
	ContainingTree TreeTwo;
	int maxspecies=2;
	TreeTwo.RandomTree(maxspecies);
	TreeTwo.ConvertTaxonNamesToRandomTaxonNumbers();
	bestscore=DoAllAssignments(bestscore,maxspecies,&TreeTwo);
	message="Now starting with 3 species";
	ContainingTree TreeThree;
	maxspecies=3;
	TreeThree.RandomTree(maxspecies);
	TreeThree.ConvertTaxonNamesToRandomTaxonNumbers();
	bestscore=DoAllAssignments(bestscore,maxspecies,&TreeThree);
	cout<<endl<<"Best trees overall"<<endl<<endl;
	for (int i=0; i<FormattedBestTrees.size(); i++) {
		(FormattedBestTrees[i]).Update();
		(FormattedBestTrees[i]).GetNodeDepths();
		(FormattedBestTrees[i]).Draw(cout);
		cout<<endl;
	}
	message="\n\n#nexus\nbegin trees;\n";
	for (int i=0; i<FormattedBestTrees.size(); i++) {
		message+="tree sptree";
		message+=i+1;
		message+=" = ";
		if (unrooted==1) {
			message+="[&U] ";
		}
		else {
			message+="[&R] ";
		}
		message+=ReturnFinalSpeciesTree(FormattedBestTrees[i]);
		message+="\n";
	}
	message+="end;\n\n";
	PrintMessage();
	message="Best score = ";
	message+=bestscore;
	PrintMessage();
	
}

double BROWNIE::DoAllAssignments(double bestscore, int maxspecies, ContainingTree *SpeciesTree )
{
	vector<double> nextscorevector;
	double nextscore;
	vector <int> countpertaxon(taxa->GetNumTaxonLabels(),1);
	for (int i=0;i<taxa->GetNumTaxonLabels();i++) {
		convertsamplestospecies[i]=1;
		cout<<" "<<convertsamplestospecies[i];
	}
	for (int assignmentrep=0;assignmentrep<pow(maxspecies,taxa->GetNumTaxonLabels());assignmentrep++) {
		if (CheckConvertSamplesToSpeciesVectorSpNum(maxspecies)) { //we have the right number of species
			nextscorevector=GetCombinedScore(SpeciesTree);
			nextscore=nextscorevector[0]; 
			cout<<"  Score: "<<nextscore<<" Best: "<<GSL_MIN(nextscore,bestscore);
			if (nextscore==bestscore) {
			//	TotalScores.push_back(nextscorevector[0]);
			//	GTPScores.push_back(nextscorevector[1]);
			//	StructScores.push_back(nextscorevector[2]);				
				FormatAndStoreBestTree(SpeciesTree,nextscorevector);
			}
			else if (nextscore<bestscore) {
				bestscore=nextscore;
				FormattedBestTrees.clear();
				TotalScores.clear();
				GTPScores.clear();
				StructScores.clear();
				BestConversions.clear();
				ContourSearchDescription.clear();
				//TotalScores.push_back(nextscorevector[0]);
				//GTPScores.push_back(nextscorevector[1]);
				//StructScores.push_back(nextscorevector[2]);
				
                FormatAndStoreBestTree(SpeciesTree,nextscorevector);
			}
		}
		else {
				cout<<"  Invalid assignment for "<<maxspecies<<" species";
		}
		cout<<endl;
	//	for (int i=0;i<countpertaxon.size();i++) {
			//cout<<" "<<countpertaxon[i];
	//	}
		//cout<<endl;
		for (int i=0;i<countpertaxon.size();i++) {
			countpertaxon[i]++;
		//	cout<<" maxspecies = "<<maxspecies<<" i = "<<i<<" pow = "<<pow(maxspecies,i)<<" countpertaxon = "<<countpertaxon[i];
			if (countpertaxon[i]>pow(maxspecies,i)) {
			//	cout<<" must change ";
				convertsamplestospecies[i]++;
				countpertaxon[i]=1;
				if (convertsamplestospecies[i]>maxspecies) {
					convertsamplestospecies[i]=1;
				}
			}
			cout<<" "<<convertsamplestospecies[i];
		}

	}	
	return bestscore;
}

void BROWNIE::DoHeuristicSearch()
{
	int chosenmove=0;
	if (jackknifesearch) {
			message="\n---------- Now starting jackknife search replicate ";
		message+=jackrep;
		message+=" ----------";
		PrintMessage();
		int includedtrees=0;
		while (includedtrees<2 || includedtrees==trees->GetNumTrees()) {
			jackknifevector.clear();
			geneidvector.clear();
			double currenttotalwt=0;
			int genenumber=1;
			for (int curtreenum=0; curtreenum<trees->GetNumTrees(); curtreenum++) {
				double newweight=trees->GetTreeWeight(curtreenum);
				currenttotalwt+=newweight;
				if (currenttotalwt>1) {
					currenttotalwt=newweight;
					genenumber++;
				}
				geneidvector.push_back(genenumber);
				jackknifevector.push_back(gsl_ran_bernoulli(r, newweight*(1-pctdelete))); //we randomly select trees
				if (jackknifevector.back()==1) {
					includedtrees++;
				}
			}
		}
		message="----------   ";
		message+=includedtrees;
		message+=" out of ";
		message+=trees->GetNumTrees();
		message+=" input trees are included.  ----------\n\n";
		PrintMessage();
	}
    message="Creating initial neighbor-joining tree for samples, based on triplet overlap. Please be patient.\n\nNow getting distances...";
    PrintMessage();

    GetTaxonTaxonTripletDistances();
    message="Computing tree...";
    PrintMessage();
    ContainingTree NJBrlenTree=ComputeTripletNJTree();
	ContainingTree TripletSupportBrlenTree=NJBrlenTree;
	TripletSupportBrlenTree.Update();
	//message="Tree of samples with triplet support brlen: ";
	//PrintMessage();
	//TripletSupportBrlenTree.Draw(cout);
	NodeIterator <Node> l (TripletSupportBrlenTree.GetRoot()); //goes from tips down
	NodePtr currentnode=l.begin();
	int nodecount=0;
	while (currentnode) {
		if (currentnode->IsLeaf()) {
			int leafnumtrans=(currentnode->GetLeafNumber())-1;
			currentnode->SetEdgeLength(CladeVectorTripletSupport[leafnumtrans]);
		}
		else if(currentnode!=TripletSupportBrlenTree.GetRoot()) {
			int nodenumber=nodecount+TripletSupportBrlenTree.GetNumLeaves();
			currentnode->SetEdgeLength(CladeVectorTripletSupport[nodenumber]);
			nodecount++;
		}
		currentnode=l.next();
	}
	if( logf_open ) {
		TripletSupportBrlenTree.Draw(logf);
		logf << endl;
		logf << "#nexus\nbegin trees;\ntree GuideTreeTripletSupportBrlen = ";
		TripletSupportBrlenTree.Write(logf);
		logf << "\ntree GuideTreeTripletSupportBrlen = ";
		NJBrlenTree.Write(logf);
		logf<<endl<<"end;"<<endl;
	}
	
	message="Proportion & type of moves:\n\t";
	message+=movefreqvector[0];
	message+="\tSubtree pruning and regrafting\n\t";
	message+=movefreqvector[1];
	message+="\tMove samples from one species to another\n\t";
	message+=movefreqvector[2];
	message+="\tIncrease the number of species\n\t";
	message+=movefreqvector[3];
	message+="\tDecrease the number of species\n\t";
	message+=movefreqvector[4];
	message+="\tReroot the species tree\n\t";
	message+=movefreqvector[5];
	message+="\tChange species tree branch lengths\n";
	PrintMessage();

	
    bestscore=GSL_POSINF;
    vector<int> intialconvertsamplestospeciesvector=convertsamplestospecies;
    double nextscore;
    double nextgtpscore;
    double nexttripletscore;
    vector<double> nextscorevector;
    nxsstring scoretype;
	message="Now starting the search proper.\nA \">\" before a score indicates that calculation of that score was aborted once the score for that move exceeded the best local score\n";
	if (!jackknifesearch) {
		if (!useCOAL && !useMS) {
			message+="\nRep\tMoves\t#Spp\tType\tQual\tCombScore\t      GTP\t   Struct\t    Local\t   Global\tNTrees\tRemaining";
		}
		else {
			if (COALaicmode==0) {
				message+="\nRep\tMoves\t#Spp\tType\tQual\t     NegLnL\t      Local\t   Global\tNTrees\tRemaining";
			}
			else if (COALaicmode==1) {
				message+="\nRep\tMoves\t#Spp\tType\tQual\t     AIC\t      Local\t   Global\tNTrees\tRemaining";
			}
			else  {
				message+="\nRep\tMoves\t#Spp\tType\tQual\t     AICc\t      Local\t   Global\tNTrees\tRemaining";
			}
		}
	}
	else {
	    message+="\nJackRep\tRep\tMoves\t#Spp\tType\tQual\tCombScore\t      GTP\t   Struct\t    Local\t   Global\tNTrees\tRemaining";
	}
    PrintMessage();

    for (int replicate=1;replicate<=nreps;replicate++) {
        convertsamplestospecies=intialconvertsamplestospeciesvector;
        ContainingTree StartingTree;
        //cout<<"Starting vector = "<<endl;
        //convertsamplestospecies=intialconvertsamplestospeciesvector;
        //for (int i=0;i<convertsamplestospecies.size();i++) {
        //		cout<<convertsamplestospecies[i]<<" ";
        //}
        //cout<<endl;
        bool assignmentbasedontriplettree=true;
        bool validassignment=false;
        if (convertsamplestospecies.size()==0) {
        	double splitwhenappropriateprob=1.0; //See where this is used below (to do initial assignment). This can be dropped if the initial assignments have too many species
			while (!validassignment) {
				convertsamplestospecies=intialconvertsamplestospeciesvector;
				if (assignmentbasedontriplettree) { //use the triplet tree to get assignments
					if (gsl_ran_flat(r,0,1)<tripletdistthreshold) {
					//New method, based on splitting on longest internal branches in starting nj tree. Basically, split on internal branches with longer than average lengths
						int nsamples=taxa->GetNumTaxonLabels();
						int numnontrivialclades=nsamples-2;
						convertsamplestospecies.assign(nsamples,1);
						if (debugmode) {
							cout<<"assembling initial convertsamplestospeciesvector, using NJ tree distances"<<endl;
							for (int k=0;k<convertsamplestospecies.size();k++) {
								cout<<convertsamplestospecies[k]<<" ";
							}
							cout<<endl;
						}
						int currentspecies=2;
						for (int currentpos=CladeVectorNJBrlen.size()-1;currentpos>=nsamples-1;currentpos--) { //start at the root, work up (based on other order on way down) 
							if (CladeVectorNJBrlen[currentpos]>meaninternalbrlen && (CladeVector[currentpos]).size()>=minsamplesperspecies && gsl_ran_bernoulli(r,splitwhenappropriateprob)==1) {
								for(int j=0;j<(CladeVector[currentpos]).size();j++) {
									convertsamplestospecies[(CladeVector[currentpos][j])]=currentspecies;
								}
								currentspecies++;
								if (debugmode) {
									for (int taxon=0; taxon<nsamples; taxon++) {
										cout<<convertsamplestospecies[taxon]<<" "<<taxa->GetTaxonLabel(taxon)<<endl;
									}
									cout<<endl;
								}
								
							}
						}
						bool goodshape=CheckConvertSamplesToSpeciesVector(true);
						goodshape=CombineSpeciesWithTooFewSamples(true);
						if (debugmode) {
							for (int taxon=0; taxon<nsamples; taxon++) {
								cout<<convertsamplestospecies[taxon]<<" "<<taxa->GetTaxonLabel(taxon)<<endl;
							}
							cout<<endl;
						}
						goodshape=CheckConvertSamplesToSpeciesTooManySpecies();
						if (!goodshape) {
							splitwhenappropriateprob=splitwhenappropriateprob*0.9; //we had too many splits, so drop the chance of splitting
						}
						if (goodshape) {
							goodshape=CheckConvertSamplesToSpeciesTooFewSpecies();
						}
						validassignment=goodshape;
						/* //Old method for getting starting assignments: tended to split good clades too often, not pay attention to relative support
							int nsamples=taxa->GetNumTaxonLabels();
						int numnontrivialclades=nsamples-2;
						convertsamplestospecies.assign(nsamples,1);
						int currentpos=nsamples-1;
						int currentspecies=2;
						while (currentpos<CladeVector.size()) {
							currentpos+=1+gsl_ran_binomial(r,.5,4);
							if (currentpos<CladeVector.size()) {
								for(int j=0;j<(CladeVector[currentpos]).size();j++) {
									convertsamplestospecies[(CladeVector[currentpos][j])]=currentspecies;
								}
							}
							currentspecies++;
						}
						bool goodshape=CheckConvertSamplesToSpeciesVector(true);
						goodshape=CombineSpeciesWithTooFewSamples(true);
						*/ //Old method for getting starting assignments
					}
					else { //Use triplet support distances
						int nsamples=taxa->GetNumTaxonLabels();
						int numnontrivialclades=nsamples-2;
						convertsamplestospecies.assign(nsamples,1);
						if (debugmode) {
							cout<<"assembling initial convertsamplestospeciesvector, using triplet support distances"<<endl;
							for (int k=0;k<convertsamplestospecies.size();k++) {
								cout<<convertsamplestospecies[k]<<" ";
							}
							cout<<endl;
						}
						int currentspecies=2;
						for (int currentpos=CladeVectorTripletSupport.size()-1;currentpos>=nsamples-1;currentpos--) { //start at the root, work up (based on other order on way down) 
							double supportvalue=CladeVectorTripletSupport[currentpos];
							if ((supportvalue>gsl_ran_flat(r,0.6,1)) && (CladeVector[currentpos]).size()>=minsamplesperspecies) { //only split on branches with 60% support or more
								for(int j=0;j<(CladeVector[currentpos]).size();j++) {
									convertsamplestospecies[(CladeVector[currentpos][j])]=currentspecies;
								}
								currentspecies++;
								if (debugmode) {
									for (int taxon=0; taxon<nsamples; taxon++) {
										cout<<convertsamplestospecies[taxon]<<" "<<taxa->GetTaxonLabel(taxon)<<endl;
									}
									cout<<endl;
								}
								
							}
						}
						bool goodshape=CheckConvertSamplesToSpeciesVector(true);
						goodshape=CombineSpeciesWithTooFewSamples(true);
						if (debugmode) {
							for (int taxon=0; taxon<nsamples; taxon++) {
								cout<<convertsamplestospecies[taxon]<<" "<<taxa->GetTaxonLabel(taxon)<<endl;
							}
							cout<<endl;
						}
						goodshape=CheckConvertSamplesToSpeciesTooManySpecies();
						if (goodshape) {
							goodshape=CheckConvertSamplesToSpeciesTooFewSpecies();
						}
						validassignment=goodshape;
					}
				}
				else {
					//message="You didn't do an intial assignment of taxa to species, so we'll do a random assignment";
					//PrintMessage();
					int ntax=taxa->GetNumTaxonLabels();
					//cout<<"ntax is "<<ntax<<endl;
					int samplesperspecies=GSL_MAX(minsamplesperspecies,1+gsl_ran_binomial(r,0.5,7)); // can change this; set now for on average 4.5 samples per species
					convertsamplestospecies.clear();
					vector <int> tempconvertsamplestospecies;
					int speciesid=1;
					int assignmentcount=0;
					for (int i=0;i<ntax;i++) {
						tempconvertsamplestospecies.push_back(speciesid);
						assignmentcount++;
						if (assignmentcount==samplesperspecies) {
							assignmentcount=0;
							speciesid++;
							samplesperspecies=GSL_MAX(minsamplesperspecies,1+gsl_ran_binomial(r,0.5,7)); //currently set for an average of 4.5 samples per species
																		   //cout<<"speciesid is "<<speciesid<<endl;
						}
					}
					if (assignmentcount<minsamplesperspecies && assignmentcount>0) { //Means the last species has too few samples in it, so we'll merge it with a random earlier species
						for (int i=0;i<tempconvertsamplestospecies.size();i++) {
							if(tempconvertsamplestospecies[i]==speciesid) {
								tempconvertsamplestospecies[i]=1+gsl_ran_binomial(r,0.5,speciesid-2);
							}
						}
					}
					gsl_permutation * c = gsl_permutation_alloc (tempconvertsamplestospecies.size());
					gsl_permutation_init (c);
					gsl_ran_shuffle (r, c->data, tempconvertsamplestospecies.size(), sizeof(size_t));
					for (int i=0; i<tempconvertsamplestospecies.size(); i++) {
						convertsamplestospecies.push_back(tempconvertsamplestospecies[gsl_permutation_get (c,i)]);
						//cout<<tempconvertsamplestospecies[gsl_permutation_get (c,i)]<<endl;
					}
					gsl_permutation_free(c);
					bool goodshape=CheckConvertSamplesToSpeciesTooManySpecies();
					if (goodshape) {
							goodshape=CheckConvertSamplesToSpeciesTooFewSpecies();
					}

					validassignment=goodshape;
				}
			}
		}
        int CurrentSppNum=0;
        for (int vectorpos=0;vectorpos<convertsamplestospecies.size();vectorpos++) {
            if (convertsamplestospecies[vectorpos]>CurrentSppNum) {
                CurrentSppNum=convertsamplestospecies[vectorpos];
            }
        }
        //cout<<"Starting vector = "<<endl;
        //for (int i=0;i<convertsamplestospecies.size();i++) {
        //	cout<<convertsamplestospecies[i]<<" ";
        //}
        //cout<<endl;
bestscorelocal=GSL_POSINF;
vector<ContainingTree> BestTreesThisRep;

//ContainingTree BestTree;
StartingTree.RandomTree(CurrentSppNum);
StartingTree.ConvertTaxonNamesToRandomTaxonNumbers();
//cout<<"Starting tree"<<endl;
//StartingTree.Write(cout);
//cout<<endl;
//StartingTree.Draw(cout);
//cout<<endl;
//cout<<"Root "<<StartingTree.GetRoot()<<endl;
//cout<<"Root child = "<<(StartingTree.GetRoot())->GetChild()<<endl;
//cout<<"Root child is leaf? "<<((StartingTree.GetRoot())->GetChild())->IsLeaf()<<endl;
//cout<<"\nTree health\n";
//StartingTree.ReportTreeHealth();
//cout<<"\n";
// char* inputforgtp=OutputForGTP(&CurrentTree);
// cout<<OutputForGTP(&StartingTree)<<endl;
// bestscorelocal=ReturnScore(OutputForGTP(&StartingTree));
//cout<<"currentsppnum = "<<CurrentSppNum<<endl;
//cout<<"starting tree health:\n";
//StartingTree.ReportTreeHealth();
if (useCOAL || useMS) {
		StartingTree.InitializeMissingBranchLengths();
}


//////////Copied from stuff below////////////////
			//brlen optimization
StartingTree.FindAndSetRoot();
StartingTree.Update();
if (useCOAL && StartingTree.GetNumLeaves()>1) {
				StartingTree.InitializeMissingBranchLengths();
				for (int brlenrep=0; brlenrep<20*(numbrlenadjustments-1); brlenrep++) {
					ContainingTree StartingTreeBrlenMod;
					StartingTreeBrlenMod.SetRoot(StartingTree.CopyOfSubtree(StartingTree.GetRoot()));
					StartingTreeBrlenMod.InitializeMissingBranchLengths();
					StartingTreeBrlenMod.RandomlyModifySingleBranchLength(markedmultiplier,brlensigma);
					vector<double> brlenscorevector=GetCombinedScore(&StartingTreeBrlenMod);
					if (brlenscorevector[0]<=bestscorelocal) {
						if (brlenscorevector[0]<bestscorelocal) {
							brlenrep=0; //so we restart from the new optimum
							bestscorelocal=brlenscorevector[0];
						}
						StartingTree.SetRoot(StartingTreeBrlenMod.CopyOfSubtree(StartingTreeBrlenMod.GetRoot()));
						StartingTree.FindAndSetRoot();
						StartingTree.Update();
						
					}
				}
				
}

			//Contour search
if (useCOAL  && StartingTree.GetNumLeaves()>1) {
				vector<double> speciestreebranchlengthvector;
				double startingwidth=contourstartingwidth; //start by looking at all brlen between pointestimate/startingwidth and startingwidth*pointestimated
				double startingnumbersteps=contourstartingnumbersteps; //works best if odd
				int maxrecursions=contourMaxRecursions;
				int recursions=0;
				int numberofedges=0;
				bool donecontour=false;
				ContainingTree StartingTreeBrlenMod;
				while (!donecontour) {
					recursions++;
					vector<double> midpointvector;
					vector<double> incrementwidths;
					vector<double> currentvector;
					speciestreebranchlengthvector.clear();
					ContourSearchVector.clear();
					StartingTree.FindAndSetRoot();
					StartingTree.Update();
					StartingTree.InitializeMissingBranchLengths();
					StartingTreeBrlenMod.SetRoot(StartingTree.CopyOfSubtree(StartingTree.GetRoot()));
					StartingTreeBrlenMod.Update();
					StartingTreeBrlenMod.InitializeMissingBranchLengths();
					vector<double> speciestreebranchlengthvector;
					NodeIterator <Node> n (StartingTreeBrlenMod.GetRoot());
					NodePtr currentnode = n.begin();
					NodePtr rootnode=StartingTreeBrlenMod.GetRoot();
					numberofedges=0;
					while (currentnode)
					{
						if (currentnode!=rootnode) {
							double edgelength=currentnode->GetEdgeLength();
							if (gsl_isnan(edgelength)) {
								edgelength=1.0;
							}
							speciestreebranchlengthvector.push_back(edgelength); //get midpoint edges
							double basestep=exp(2.0*log(startingwidth)/(startingnumbersteps-1));
							double smallestbrlen=edgelength*(pow(basestep,(0-startingnumbersteps+((startingnumbersteps+1.0)/2.0))));
							midpointvector.push_back(edgelength);
							currentvector.push_back(smallestbrlen); //start with minimum values, then  move up
							incrementwidths.push_back(basestep);
							numberofedges++;
						}
						currentnode = n.next();
					}
					bool donegrid=false;
					vector<int> increments;
					increments.assign(numberofedges,0);
					int origCOALaicmode=COALaicmode;
					COALaicmode=0;
					vector<double> startingscorevector=GetCombinedScore(&StartingTree);
					double currentscore=startingscorevector[0];
					while (!donegrid) {
						//cout<<"Looping over grid tries"<<endl;
						currentnode = n.begin();
						int nodenumber=0;
						while (currentnode)
						{
							if (currentnode!=rootnode) {
								currentnode->SetEdgeLength(currentvector[nodenumber]);
								nodenumber++;
							}
							currentnode = n.next();
						}
						vector<double> brlenscorevector=GetCombinedScore(&StartingTreeBrlenMod);
						double brlenscore=brlenscorevector[0]-currentscore;
						bool atmargin=false;
						for (int i=0; i<numberofedges; i++) {
							if ((increments[i]==0) || (increments[i]==startingnumbersteps-1)) {
								atmargin=true;
								
							}
						}
						if (atmargin) { //if we're bumping up against minimum branchlength, there's nowhere further to expand the grid
							for (int i=0; i<numberofedges; i++) {
								if (currentvector[i]==0) {
									atmargin=false;
									
								}
							}
						}
						//cout<<endl;
						if (atmargin) { //we're at a margin of the space; want to make sure that the region within two lnL is inside this region
							if (brlenscore<2 && recursions<maxrecursions) { //our region is too small, since points 2 lnL units away from the max are outside the region
								startingwidth*=1.5;
								donegrid=true; //break the while(!donegrid) loop; since donecontour isn't done, reinitialize everything
							}
						}
						vector<double> resultvector=currentvector;
						resultvector.push_back(brlenscore);
						ContourSearchVector.push_back(resultvector);
						if (brlenscore<0 && recursions<=maxrecursions) {
							currentscore=brlenscorevector[0];
							StartingTree.SetRoot(StartingTreeBrlenMod.CopyOfSubtree(StartingTreeBrlenMod.GetRoot()));
							StartingTree.FindAndSetRoot();
							StartingTree.Update();
							donegrid=true; //break the while(!donegrid) loop; since donecontour isn't done, reinitialize everything
							if (showtries) {
								cout<<"Better branch lengths found in grid search"<<endl;
							}
						}
						increments[0]++;
						if (!donegrid) {
							for (int itemtoexamine=0; itemtoexamine<numberofedges; itemtoexamine++) {
								if (increments[itemtoexamine]==startingnumbersteps) {
									if (itemtoexamine<numberofedges-1) { //means there's room to the left
										increments[itemtoexamine]=0;
										increments[itemtoexamine+1]++;
									}
									else {
										donegrid=true;
										donecontour=true;
									}
								}
							}
						}
						currentvector.clear();
						for (int i=0; i<numberofedges; i++) {
							double newbrlen=midpointvector[i]*(pow(incrementwidths[i],(increments[i]-startingnumbersteps+((startingnumbersteps+1)/2))));
							currentvector.push_back(newbrlen);
						}
					}
					
					COALaicmode=origCOALaicmode;
					
				}
				vector<double> totalbrlen;
				totalbrlen.assign(numberofedges-1,0);
				int numequaltrees=0;
				for (int i=0; i<ContourSearchVector.size(); i++) {
					if (ContourSearchVector[i][numberofedges]<=0) {
						for (int j=0; j<numberofedges; j++) {
							totalbrlen[j]+=ContourSearchVector[i][j];
							numequaltrees++;
						}
					}
				}
NodeIterator <Node> n (StartingTree.GetRoot());
NodePtr currentnode = n.begin();
NodePtr rootnode=StartingTree.GetRoot();
int edgenumber=0;
while (currentnode)
{
	if (currentnode!=rootnode) {
						//cout<<"new brlen = "<<totalbrlen[edgenumber]/(numequaltrees*1.0)<<endl;
		currentnode->SetEdgeLength(totalbrlen[edgenumber]/(numequaltrees*1.0));
		edgenumber++;
	}
	currentnode = n.next();
}
nextscorevector=GetCombinedScore(&StartingTree);
nextscore=nextscorevector[0];
}
//////////END Copied from stuff below////////////////

vector<double> bestscorelocalvector=GetCombinedScore(&StartingTree);
bestscorelocal=bestscorelocalvector[0];
if (bestscorelocal==bestscore) {
	//cout<<"Before RawBestTrees.push_back(StartingTree);"<<endl;
	//StartingTree.ReportTreeHealth();
    RawBestTrees.push_back(StartingTree);
	//TotalScores.push_back(bestscorelocalvector[0]);
	//GTPScores.push_back(bestscorelocalvector[1]);
	//StructScores.push_back(bestscorelocalvector[2]);	
    FormatAndStoreBestTree(&StartingTree,bestscorelocalvector);
    scoretype="*G\t";
}
else if (bestscorelocal<bestscore) {
    RawBestTrees.clear();
	//cout<<"Before RawBestTrees.push_back(StartingTree);"<<endl;
	//StartingTree.ReportTreeHealth();
    RawBestTrees.push_back(StartingTree);
    FormattedBestTrees.clear();
	TotalScores.clear();
	GTPScores.clear();
	StructScores.clear();
	BestConversions.clear();
	ContourSearchDescription.clear();
//	TotalScores.push_back(bestscorelocalvector[0]);
//	GTPScores.push_back(bestscorelocalvector[1]);
//	StructScores.push_back(bestscorelocalvector[2]);
    FormatAndStoreBestTree(&StartingTree,bestscorelocalvector);
    bestscore=bestscorelocal;
    scoretype="=G\t";
}
else {
    scoretype="*L\t";
}

//while (bestscorelocal<0) { //due to error in GTP
//    bestscorelocal=ReturnScore(OutputForGTP(&StartingTree));
//}
//cout<<"Before BestTreesThisRep.push_back(StartingTree);"<<endl;
//StartingTree.ReportTreeHealth();

BestTreesThisRep.push_back(StartingTree);

//BestTree=StartingTree;
int movecount=0;
if (status) {
    //cout<<"\n\n"<<OutputForGTP(&CurrentTree)<<"\n\n";
    // cout<<"Starting tree: \n\n"; //Rewrite the draw function to allow output to a file
    // CurrentTree.Draw(cout);
    // cout<<"\nScore is "<<bestscore<<"\n";
}
// if (BestTrees.size()==0) {
//    BestTrees.push_back(CurrentTree);
// }
bool improvement=true;
		bool moreswaps=true;
		bool morereassignments=true;
		bool moreincreases=true;
		bool moredecreases=true;
		bool morererootings=true;
		
while (improvement && (rearrlimit<0 || movecount<rearrlimit)) {
    //cout<<"\nimprovement, restarting\n";
    improvement=false;
	if (chosenmove!=6) { //only reset moves on topology change
		bool moreswaps=true;
		bool morereassignments=true;
		bool moreincreases=true;
		bool moredecreases=true;
		bool morererootings=true;
	}
    // cout<<"moreswaps = "<<moreswaps<<" morereassignments = "<<morereassignments<<" moreincreases = "<<moreincreases<<" moredecreases = "<<moredecreases<<" morererootings = "<<morererootings<<endl;
    assert(BestTreesThisRep.size()>0);
    //for(int i=0;i<BestTreesThisRep.size();i++) {
    //    BestTreesThisRep[i].Write(cout);
    //    cout<<endl;
    // }
    // (BestTreesThisRep.back()).Draw(cout);
    // (BestTreesThisRep.back()).ReportTreeHealth();
    (BestTreesThisRep.back()).Update();
    //(BestTreesThisRep.back()).ReportTreeHealth();
    // cout<<"GetRoot: "<<(BestTreesThisRep.back()).GetRoot()<<endl;
    assert(BestTreesThisRep.size()>0);
	//cout<<"Before  ContainingTree CurrentTree=BestTreesThisRep.back();"<<endl;
	//(BestTreesThisRep.back()).ReportTreeHealth();
    ContainingTree CurrentTree=BestTreesThisRep.back();
    CurrentTree.Update();
    CurrentTree.ResetBreakVector();
    CurrentTree.UpdateCherries();
    CurrentTree.SetLeafNumbers();
    if ((sppnumfixed==true) || CurrentTree.GetNumLeaves()<=minnumspecies) {
        moredecreases=false;
    }
    if ((sppnumfixed==true) || CurrentTree.GetNumLeaves()>=maxnumspecies) {
        moreincreases=false;
    }
    if (movefreqvector[0]==0 || CurrentTree.GetNumLeaves()<3) {
        moreswaps=false;
    }
    if (movefreqvector[1]==0 || CurrentTree.GetNumLeaves()==1) {
        morereassignments=false;
    }
    if (movefreqvector[2]==0) {
        moreincreases=false;
    }
    if (movefreqvector[3]==0 || CurrentTree.GetNumLeaves()==1) {
        moredecreases=false;
    }
    if (movefreqvector[4]==0 || CurrentTree.GetNumLeaves()<3) {
        morererootings=false;
    }
    //Get list of cherries to collapse; we do this at the start so that we try each cherry at random but only once.
    vector<int> TempCherriesToMash;
    vector<int> CherriesToMash;
    for (int i=0;i<CurrentTree.GetNumCherries(); i++) {
        TempCherriesToMash.push_back(i);
    }
    if (CurrentTree.GetNumCherries()>0) {
        gsl_permutation * p = gsl_permutation_alloc (TempCherriesToMash.size());
        gsl_permutation_init (p);
        gsl_ran_shuffle (r, p->data, TempCherriesToMash.size(), sizeof(size_t));
        for (int i=0; i<TempCherriesToMash.size(); i++) {
            CherriesToMash.push_back(TempCherriesToMash[gsl_permutation_get (p,i)]);
        }
        gsl_permutation_free(p);
    }
    else {
        moredecreases=false;
    }

    //Get list of nodes to reroot on
    vector<int> NodesToReRootOn=CurrentTree.GetPotentialNewRoots();
    if (NodesToReRootOn.size()==0) {
        morererootings=false;
    }
    //cout<<"There are potentially "<<NodesToReRootOn.size()<<" new roots\n";


    //Get list of leaves to split; we do this at the start so that we try each possible leaf (leaves with at least two samples) at random but only once.
    vector<int> TempLeavesToSplit;
    vector<int> LeavesToSplit;
    for (int i=0;i<CurrentTree.GetNumLeaves(); i++) {
        int numsamples=0;
        for (int j=0;j<convertsamplestospecies.size();j++) {
            if(convertsamplestospecies[j]==i+1) {
                numsamples++;
            }
        }
        if (numsamples>minsamplesperspecies) {
            TempLeavesToSplit.push_back(i+1);
        }
    }
	if (TempLeavesToSplit.size()>0) {
		gsl_permutation * q = gsl_permutation_alloc (TempLeavesToSplit.size());
		gsl_permutation_init (q);
		gsl_ran_shuffle (r, q->data, TempLeavesToSplit.size(), sizeof(size_t));
		for (int i=0; i<TempLeavesToSplit.size(); i++) {
			LeavesToSplit.push_back(TempLeavesToSplit[gsl_permutation_get (q,i)]);
		}
		gsl_permutation_free(q);
		if (showtries) {
			cout<<"Made LeavesToSplitVector of size "<<LeavesToSplit.size()<<endl<<"contents: ";
			for (int k=0;k<LeavesToSplit.size();k++) {
					cout<<LeavesToSplit[k]<<"\t";
			}
			cout<<endl;
		}
	}
	else {
		moreincreases=false;
	}

    SamplesToMove.clear();
    vector<int> TempSamplesToMove;
    int maxsamplesperspecies=0;
    vector<int> SamplesPerSpecies(1+CurrentTree.GetNumLeaves(),0); //so SamplesPerSpecies[0] is empty but then SamplesPerSpecies[X] is the number of samples for species X
    for (int i=0;i<convertsamplestospecies.size();i++) {
        SamplesPerSpecies[(convertsamplestospecies[i])]++;
        maxsamplesperspecies=GSL_MAX(maxsamplesperspecies,SamplesPerSpecies[(convertsamplestospecies[i])]);
    }
    if (maxsamplesperspecies<=minsamplesperspecies || movefreqvector[1]==0) {
        morereassignments=false;
    }
    else {
        vector<int> TempSampleDestinations;
        for (int i=0;i<CladeVector.size();i++) {
            TempSamplesToMove.push_back(i);
        }

        gsl_permutation * u = gsl_permutation_alloc (TempSamplesToMove.size());
        gsl_permutation_init (u);
        gsl_ran_shuffle (r, u->data, TempSamplesToMove.size(), sizeof(size_t));
        for (int i=0; i<TempSamplesToMove.size(); i++) {
            SamplesToMove.push_back(TempSamplesToMove[gsl_permutation_get (u,i)]);
        }
        gsl_permutation_free(u);

        SampleDestinations.clear();
        bool enoughdestinations=false;
        while (!enoughdestinations && SamplesToMove.size()>0) {
            SampleDestinations.clear();

            for (int i=0;i<CurrentTree.GetNumLeaves();i++) {
                if(TestMoveSamples(SamplesToMove.back(),i+1)) {
                    TempSampleDestinations.push_back(i+1);
                    enoughdestinations=true;
                }
            }
            if (!enoughdestinations) {
                SamplesToMove.pop_back(); //we're done trying to move that sample
            }
        }
        morereassignments=enoughdestinations;
        if(enoughdestinations) {
            gsl_permutation * v = gsl_permutation_alloc (TempSampleDestinations.size());
            gsl_permutation_init (v);
            gsl_ran_shuffle (r, v->data, TempSampleDestinations.size(), sizeof(size_t));
            for (int i=0; i<TempSampleDestinations.size(); i++) {
                SampleDestinations.push_back(TempSampleDestinations[gsl_permutation_get (v,i)]);
            }
            gsl_permutation_free(v);
        }
    }
    //cout<<"moreswaps = "<<moreswaps<<" morereassignments = "<<morereassignments<<" moreincreases = "<<moreincreases<<" moredecreases = "<<moredecreases<<" morererootings = "<<morererootings<<endl;

    if (movecount==0) {

        message="";
		if (jackknifesearch) {
			message+=jackrep;
			message+="\t";
		}
        message+=replicate;
        message+="\t";
        message+=movecount;
        message+="\t";
        message+=CurrentTree.GetNumLeaves();
        message+="\t\t";
        message+=scoretype;
        char outputstring[9];
        sprintf(outputstring,"%9.3f",bestscorelocal);
        message+=outputstring;
		if (!useCOAL && !useMS) {
			message+="\t";
			sprintf(outputstring,"%9.3f",bestscorelocalvector[1]);
			message+=outputstring;
			message+="\t";
			sprintf(outputstring,"%9.3f",bestscorelocalvector[2]);
			message+=outputstring;
		}
        message+="\t";
        sprintf(outputstring,"%9.3f",bestscorelocal);
        message+=outputstring;
        message+="\t";
        sprintf(outputstring,"%9.3f",GSL_MIN(bestscore,bestscorelocal));
        message+=outputstring;
        message+="\t";
        message+=int(FormattedBestTrees.size());
        if (moreswaps) {
            message+="\ts";
        }
        else {
            message+="\t_";
        }
        if (morereassignments) {
            message+="a";
        }
        else {
            message+="_";
        }
        if (moreincreases) {
            message+="i";
        }
        else {
            message+="_";
        }
        if (moredecreases) {
            message+="d";
        }
        else {
            message+="_";
        }
        if (morererootings) {
            message+="r";
        }
        else {
            message+="_";
        }
		if (movefreqvector[5]>0) {
			message+="b";
		}
        if (status) {
            PrintMessage();
        }
    }

    while ((moreswaps || morereassignments || moreincreases || moredecreases || morererootings) && (rearrlimit<0 || movecount<rearrlimit)) {
		//cout<<"while ((moreswaps || morereassignments || moreincreases || moredecreases || morererootings) && (rearrlimit<0 || movecount<rearrlimit)) {"<<endl;
        bool somethinghappened=true;
		//cout<<"just before ContainingTree NextTree=CurrentTree"<<endl;
		CurrentTree.FindAndSetRoot();
		CurrentTree.Update();
		//CurrentTree.ReportTreeHealth();
        ContainingTree NextTree=CurrentTree;
		//cout<<"just after ContainingTree NextTree=CurrentTree"<<endl;
		//cout<<"just before ContainingTree BrlenChangedTree"<<endl;
		ContainingTree BrlenChangedTree;
		//cout<<"just after ContainingTree BrlenChangedTree"<<endl;
        //   cout<<"\n\nOldTree\n"<<ReturnFinalSpeciesTree(CurrentTree)<<endl;
        // for (int i=0;i<convertsamplestospecies.size();i++) {
        //      cout<<convertsamplestospecies[i]<<" ";
        //  }
        // cout<<endl;
        NextTree.UpdateCherries();
        NextTree.Update();
        NextTree.GetNodeDepths();
        vector<int> Originalconvertsamplestospecies=convertsamplestospecies;

        //decide chosen move
        if (CurrentTree.GetNumLeaves()==1) {
            moreswaps=false;
            morereassignments=false;
            moredecreases=false;
            morererootings=false;
        }
        double randomvalue=double(gsl_ran_flat (r,0,1));
        chosenmove=0;
        nxsstring chosenmovestring="?";
        vector<double> possiblemovefreqvector;
        vector<int> possiblemovechoicevector;
        vector<nxsstring> possiblemoveabbrevvector;
        if (moreswaps) {
            possiblemovefreqvector.push_back(movefreqvector[0]);
            possiblemovechoicevector.push_back(1);
            possiblemoveabbrevvector.push_back("s");
        }
        if (morereassignments) {
            possiblemovefreqvector.push_back(movefreqvector[1]);
            possiblemovechoicevector.push_back(2);
            possiblemoveabbrevvector.push_back("a");
        }
        if (moreincreases) {
            possiblemovefreqvector.push_back(movefreqvector[2]);
            possiblemovechoicevector.push_back(3);
            possiblemoveabbrevvector.push_back("i");
        }
        if (moredecreases) {
            possiblemovefreqvector.push_back(movefreqvector[3]);
            possiblemovechoicevector.push_back(4);
            possiblemoveabbrevvector.push_back("d");
        }
        if (morererootings) {
            possiblemovefreqvector.push_back(movefreqvector[4]);
            possiblemovechoicevector.push_back(5);
            possiblemoveabbrevvector.push_back("r");
        }
		possiblemovefreqvector.push_back(movefreqvector[5]); //the brlen optimization
		if (movefreqvector[5]>0) {
			possiblemovechoicevector.push_back(6);
            possiblemoveabbrevvector.push_back("b");
		}
        double sumofpossiblemovefreqs=0;
        for (int k=0; k<possiblemovefreqvector.size(); k++) {
            sumofpossiblemovefreqs+=possiblemovefreqvector[k];
        }
        for (int k=0; k<possiblemovefreqvector.size(); k++) {
            possiblemovefreqvector[k]=(possiblemovefreqvector[k])/sumofpossiblemovefreqs;
            //  cout<<"possiblemovefreqvector["<<k<<"] = "<<possiblemovefreqvector[k]<<"\tpossiblemovechoicevector["<<k<<"] = "<<possiblemovechoicevector[k]<<endl;
        }
        double runningtotal=0;
        for (int k=0; k<possiblemovefreqvector.size(); k++) {
            runningtotal+=possiblemovefreqvector[k];
            // cout<<"randomvalue = "<<randomvalue<<" runningtotal = "<<runningtotal;
            if (randomvalue<=runningtotal) {
                chosenmove=possiblemovechoicevector[k];
                chosenmovestring=possiblemoveabbrevvector[k];
                //cout<<" chosenmove is "<<chosenmove;
                break;
            }
            //cout<<endl;
        }
        //cout<<"randomvalue is "<<randomvalue<<" chosenmove is "<<chosenmove<<" moreswaps = "<<moreswaps<<" morereassignments = "<<morereassignments<<" moreincreases = "<<moreincreases<<" moredecreases = "<<moredecreases<<" morererootings = "<<morererootings<<endl;
        movecount++;
        //  cout<<"convertsamplestospecies\n";
        //  for (int m=0;m<convertsamplestospecies.size();m++) {
        //      cout<<convertsamplestospecies[m]<<" ";
        //  }
        //  cout<<endl;
        //cout<<"chosenmove = "<<chosenmove<<endl;
        if (chosenmove==1) { //Try branch swap
            if (showtries) {
                cout<<"Trying branch swap"<<endl;
                cout<<"Start tree = \n";
                NextTree.Draw(cout);
            }
            NextTree.SetBreakVector(CurrentTree.GetBreakVector());
            NextTree.SetAttachVector(CurrentTree.GetAttachVector());
            NextTree.FindAndSetRoot();
            NextTree.Update();
            moreswaps=NextTree.NextSPR();
            if (showtries) {
                cout<<"Swap tree = \n";
                NextTree.Draw(cout);
            }
            CurrentTree.SetBreakVector(NextTree.GetBreakVector()); //due to how the vectors are updated during a swap.
            CurrentTree.SetAttachVector(NextTree.GetAttachVector());
			if (useCOAL || useMS) {
				NextTree.InitializeMissingBranchLengths();
			}
            nextscorevector=GetCombinedScore(&NextTree);
            nextscore=nextscorevector[0];
        }
        else if(chosenmove==2) {
            if (showtries) {
                cout<<"Moving a sample from one species to another\n";
            }
            if (showtries) {
                cout<<"Start assignment = (";
                for (int i=0;i<(CladeVector[SamplesToMove.back()]).size();i++) {
                    cout<<" "<<CladeVector[SamplesToMove.back()][i];
                }
                cout<<" )\n";
                for (int i=0; i<convertsamplestospecies.size();i++) {
                    cout<<convertsamplestospecies[i]<<" ";
                }
                cout<<"\n";
            }
            morereassignments=MoveSamples(Originalconvertsamplestospecies);
            if (showtries) {
                for (int i=0; i<convertsamplestospecies.size();i++) {
                    cout<<convertsamplestospecies[i]<<" ";
                }
                cout<<"\n";
            }
			if (useCOAL || useMS) {
				NextTree.InitializeMissingBranchLengths();
			}
            nextscorevector=GetCombinedScore(&NextTree);
            nextscore=nextscorevector[0];
        }
        else if(chosenmove==3) {
            //cout<<"LeavesToSplitVect\n";
            //for (int i=0;i<LeavesToSplit.size();i++) {
            //     cout<<" "<<LeavesToSplit[i];
            // }
            //increase the number of species
            if (showtries) {
				cout<<"LeavesToSplitVect\n";
				cout<<"vector size is "<<LeavesToSplit.size()<<endl;
				for (int i=0;i<LeavesToSplit.size();i++) {
					cout<<" "<<LeavesToSplit[i];
				}
                cout<<"Splitting a leaf\n";
            }
            int ChosenLeaf=LeavesToSplit.back();
            //cout<<"\nChosenLeaf = "<<ChosenLeaf<<endl;
            LeavesToSplit.pop_back();
            // cout<<"\nVector size now "<<LeavesToSplit.size();
            if (LeavesToSplit.size()==0) {
                moreincreases=false;
            }
            // cout<<"\nmoreincreases value = "<<moreincreases<<endl;
            vector<int>changevector;
            if (showtries) {
                cout<<"Split leaf start tree = \n";
                NextTree.Draw(cout);
            }
            changevector=NextTree.SplitLeaf(ChosenLeaf); //first element is split taxon, second element is new taxon
            if (showtries) {
                cout<<"Final tree = \n";
                NextTree.Draw(cout);
            }
			if (useCOAL || useMS) {
				NextTree.InitializeMissingBranchLengths();
			}
            //now try optimizing the new assignments (which descendant the samples go with) before actually getting the score
            vector<int> samplestomove;
            int sampletostay;
            for (int i=0;i<convertsamplestospecies.size();i++) {
                //cout<<convertsamplestospecies[i]<<" ";
                if(convertsamplestospecies[i]==changevector[0]) {
                    samplestomove.push_back(i);
                }
            }
            // cout<<endl;
            sampletostay=samplestomove.back();
            samplestomove.pop_back();
            //so idea here is to try all combinations, with one sample fixed in the old species and the others allowed to be in either species
            vector<int> Startingconvertsamplestospecies=convertsamplestospecies;
            vector<int> Bestconvertsamplestospecies=convertsamplestospecies;
            bool bestscorefound=false;
            double bestscoreforcombination=GSL_POSINF;
			vector<double> lastscorevector;
            //  cout<<"leaf split starting assignment"<<endl;
            //  for (int i=0;i<convertsamplestospecies.size();i++) {
            //      cout<<convertsamplestospecies[i]<<" ";
            //   }
            //   cout<<endl;
            size_t j;
            double numberofcomparisons=0;
            for(int l=1;l<=samplestomove.size();l++){
                numberofcomparisons+=gsl_sf_choose(samplestomove.size(),l);
				if (showtries) {
					cout<<"numberofcomparisons now "<<numberofcomparisons<<endl;
				}
				if (numberofcomparisons<1) { 
					cout<<"Number of comparisons was "<<numberofcomparisons<<" and samplestomove.size() was "<<samplestomove.size()<<endl;
				}
				assert(numberofcomparisons>0);
            }
            double probofacomb=(pow((1.0*numberofcomparisons),1.0/chosensubsampling))/(1.0*numberofcomparisons); //a way to reduce the search effort
			if (probofacomb*numberofcomparisons<10) { //so that if we choose a ridiculous number we expect to do at least ten swaps
				probofacomb=1;
			}
			else if (probofacomb!=probofacomb) {
				if (showtries) {
					cout<<"prob of a comb is "<<probofacomb<<" so we're adjusting it"<<endl;
				}
				probofacomb=GSL_MIN(100.0/numberofcomparisons,1);
			}
			if (showtries) {
				cout<<"Prob of a comb is "<<probofacomb<<" expected number of assignments to examine is "<<probofacomb*numberofcomparisons<<endl;
			}
            int combinationmoves=0;
            //while(bestscorefound==false) {
			gsl_combination *c;
			vector<int> LastTriedconvertsamplestospecies;
			for (j=1;j<=samplestomove.size();j++) { //always move
				c=gsl_combination_calloc(samplestomove.size(),j);
				do
				{
					combinationmoves++;
					if (gsl_ran_bernoulli(r,probofacomb)==1 || steepest || exhaustive) {
						if (showtries) {
							
							cout<<combinationmoves<<"/"<<numberofcomparisons<<" = ";
							cout<<(1.0*combinationmoves)/(1.0*numberofcomparisons)<<endl;
                                //ProgressBar(0);
							
						}
						convertsamplestospecies=Startingconvertsamplestospecies;
						for (int k=0;k<j;k++) {
							convertsamplestospecies[samplestomove[int(gsl_combination_get(c,k))]]=changevector[1]; //assign this taxon to the new species
						}
						LastTriedconvertsamplestospecies=convertsamplestospecies;
                            //  for (int m=0;m<convertsamplestospecies.size();m++) {
                            //      cout<<convertsamplestospecies[m]<<" ";
                            //   }
                            //for (int i=0;i<convertsamplestospecies.size();i++) {
                            //     cout<<convertsamplestospecies[i]<<" ";
                            // }
                            // cout<<endl;
                            //    for (int i=0;i<convertsamplestospecies.size();i++) {
                            //      cout<<convertsamplestospecies[i]<<" ";
                            //   }
                            //  cout<<endl;
						vector<double> newscoreforcombinationvector=GetCombinedScore(&NextTree);
						double newscoreforcombination=newscoreforcombinationvector[0];
						lastscorevector.swap(newscoreforcombinationvector);
						if (showtries) {
							cout<<"Score "<<newscoreforcombination<<" assign ( ";
							for (int i=0; i<convertsamplestospecies.size();i++) {
								cout<<convertsamplestospecies[i]<<" ";
							}
							cout<<" )"<<endl;
						}	
						if ((newscoreforcombination<bestscoreforcombination) && (1==gsl_finite(newscoreforcombination))) { //the second condition makes sure it isn't nan or inf
							if (showtries) {
								cout<<" Above is BETTER"<<endl;
							}
							bestscoreforcombination=newscoreforcombination;
							Bestconvertsamplestospecies=convertsamplestospecies;
							if (!steepest && !exhaustive) {
								bestscorefound=true;
							}
						}
                            // cout<<"\t"<<newscoreforcombination<<endl;
					}
				}
				while ((gsl_combination_next (c) == GSL_SUCCESS) && (bestscorefound==false));
				gsl_combination_free(c);
			}
			//gsl_combination_free(c); /moved up to stop leak
           // }
			
			
            convertsamplestospecies.swap(Bestconvertsamplestospecies);
            //  cout<<"Final assignment after split leaf"<<endl;
            //  for (int i=0;i<convertsamplestospecies.size();i++) {
            //      cout<<convertsamplestospecies[i]<<" ";
            //  }
            //  cout<<endl;
            //  cout<<"Bestscorefound = "<<bestscorefound<<endl;
            //convertsamplestospecies=Originalconvertsamplestospecies;
			if(isinf(bestscoreforcombination)==0) {//so the best score is NOT infinity
				if (useCOAL || useMS) {
					NextTree.InitializeMissingBranchLengths();
				}
				nextscorevector=GetCombinedScore(&NextTree);
				nextscore=nextscorevector[0];
			}
			else {
				convertsamplestospecies=LastTriedconvertsamplestospecies; //otherwise, an n-species tree is evaluated with the original convertsamplestospecies vector, which has n-1 species
				nextscorevector.swap(lastscorevector);
				nextscore=nextscorevector[0];
			}
        }
        else if(chosenmove==4) {  //reduce the number of species, if possible
                                  //cout<<"\nCherry vector"<<endl;
                                  //for (int i=0;i<CherriesToMash.size(); i++) {
                                  //    cout<<" "<<CherriesToMash[i];
                                  //}
                                  //cout<<endl;
            if (showtries) {
                cout<<"Collapsing a cherry\n";
            }
            int CherryToMash=CherriesToMash.back();
            //cout<<"Cherry to mash = "<<CherryToMash<<endl;
            CherriesToMash.pop_back(); //this is so we look at each cherry once
            if (CherriesToMash.size()==0) {
                moredecreases=false;
            }
            vector<int>changevector;
            if (showtries) {
                cout<<"Cherry collapse start tree = \n";
                NextTree.Draw(cout);
            }
            changevector=NextTree.CollapseCherry(CherryToMash);
            for (int i=0; i<convertsamplestospecies.size(); i++) {
                if(convertsamplestospecies[i]==changevector[1]) {
                    convertsamplestospecies[i]=changevector[0];
                }
                if(convertsamplestospecies[i]>changevector[1]) {
                    convertsamplestospecies[i]--; //so if we have taxa 1-8, and delete taxon 6, taxon 7 becomes the new 6 and taxon 8 becomes the new 7
                }
            }
            if (showtries) {
                cout<<"Next tree = \n";
                NextTree.Draw(cout);
            }
			if (useCOAL || useMS) {
				NextTree.InitializeMissingBranchLengths();
			}
            nextscorevector=GetCombinedScore(&NextTree);
            nextscore=nextscorevector[0];
        }
		
        else if (chosenmove==5) {
            if (showtries) {
                cout<<"Rerooting"<<endl;
            }
            int NodeToReRootOnNum=NodesToReRootOn.back();
            NodesToReRootOn.pop_back(); //this is so we look at each potential position once
            if (NodesToReRootOn.size()==0) {
                morererootings=false;
            }
            vector<int>changevector;
            if (showtries) {
                cout<<"Rerooting start tree = \n";
                NextTree.Draw(cout);
            }
            NextTree.ReRootTree(NextTree.SelectNodeToReRootOn(NodeToReRootOnNum));
            if (showtries) {
                cout<<"Next tree = \n";
                NextTree.Draw(cout);
            }
			if (useCOAL || useMS) {
				NextTree.InitializeMissingBranchLengths();
			}
            nextscorevector=GetCombinedScore(&NextTree);
            nextscore=nextscorevector[0];
        }
		else if (chosenmove==6) {
			if (showtries) {
				cout<<"Doing branch length optimization\n";
			}
			NextTree.FindAndSetRoot();
			NextTree.Update();
			NextTree.InitializeMissingBranchLengths();
			if (useCOAL) {
				NextTree.RandomlyModifySingleBranchLength(markedmultiplier,brlensigma);
			}
			else if (useMS) {
				if (0.2>gsl_ran_flat (r,0,1) || NextTree.GetNumLeaves()<3) {
					NextTree.ModifyTotalBranchLength(brlensigma);
				}
				else {
					NextTree.NodeSlideBranchLength(markedmultiplier);
				}
			}
			else {
				message+="Warning: attempting to do branch length estimation with a criterion that doesn't take it into account";
				PrintMessage();
			}
			nextscorevector=GetCombinedScore(&NextTree);
            nextscore=nextscorevector[0];
		}
        else {
            somethinghappened=false;
            movecount--;
        }
        if (somethinghappened) {
			//cout<<"Something happened"<<endl;
/*            if (NextTree.GetNumLeaves()==NextTree.GetNumInternals()) {
                errormsg="Error: num leaves = num internals\nLast move chosen was";
                errormsg+=chosenmove;
                NextTree.ReportTreeHealth();
                throw XNexus( errormsg);
				
            }*/
            assert(CheckConvertSamplesToSpeciesVector(false));

            scoretype="\t";
            bool modifiedscoretype=false;
			NextTree.FindAndSetRoot();
			NextTree.Update();

			
			//brlen optimization
			if (useCOAL && NextTree.GetNumLeaves()>1) { //only do this is there are at least two species (brlen doesn't matter for single species);
				NextTree.InitializeMissingBranchLengths();
				for (int brlenrep=0; brlenrep<numbrlenadjustments; brlenrep++) {
					BrlenChangedTree.SetRoot(NextTree.CopyOfSubtree(NextTree.GetRoot()));
					BrlenChangedTree.InitializeMissingBranchLengths();
					BrlenChangedTree.RandomlyModifySingleBranchLength(markedmultiplier,brlensigma);
					vector<double> brlenscorevector=GetCombinedScore(&BrlenChangedTree);
					if (brlenscorevector[0]<=nextscore) {
						if (brlenscorevector[0]<nextscore) {
							brlenrep=0; //so we restart from the new optimum
							nextscore=brlenscorevector[0];
							nextscorevector[0]=brlenscorevector[0]; //other elements are the same
						}
						NextTree.SetRoot(BrlenChangedTree.CopyOfSubtree(BrlenChangedTree.GetRoot()));
						NextTree.FindAndSetRoot();
						NextTree.Update();
						
					}
				}
				
			}
			 
			
			//Contour search
			if (useCOAL  && NextTree.GetNumLeaves()>1) {
				vector<double> speciestreebranchlengthvector;
				double startingwidth=contourstartingwidth; //start by looking at all brlen between pointestimate/startingwidth and startingwidth*pointestimated
				double startingnumbersteps=contourstartingnumbersteps; //works best if odd
				int maxrecursions=contourMaxRecursions;
				int recursions=0;
				int numberofedges=0;
				bool donecontour=false;
				while (!donecontour) {
					recursions++;
					//cout<<"donecontour"<<endl;
					vector<double> midpointvector;
					vector<double> incrementwidths;
					vector<double> currentvector;
					speciestreebranchlengthvector.clear();
					ContourSearchVector.clear();
					NextTree.FindAndSetRoot();
					NextTree.Update();
					NextTree.InitializeMissingBranchLengths();
					BrlenChangedTree.SetRoot(NextTree.CopyOfSubtree(NextTree.GetRoot()));
					BrlenChangedTree.Update();
					BrlenChangedTree.InitializeMissingBranchLengths();
					//BrlenChangedTree.ReportTreeHealth();
					vector<double> speciestreebranchlengthvector;
					NodeIterator <Node> n (BrlenChangedTree.GetRoot());
					NodePtr currentnode = n.begin();
					NodePtr rootnode=BrlenChangedTree.GetRoot();
					numberofedges=0;
					while (currentnode)
					{
						if (currentnode!=rootnode) {
							double edgelength=currentnode->GetEdgeLength();
							if (gsl_isnan(edgelength)) {
								edgelength=1.0;
							}
							speciestreebranchlengthvector.push_back(edgelength); //get midpoint edges
							double basestep=exp(2.0*log(startingwidth)/(startingnumbersteps-1));
							//cout<<"basestep = "<<basestep<<endl;
							//brlen if startingnumbersteps=5 and starting width is 4 is edgelength/4, edgelength/2, edgelength, edgelength*2, edgelength*4; aka 2^-2, 2^-1, 
							//double mindepth=GSL_MAX(edgelength*(1-startingwidth),0);
							//double maxdepth=edgelength*(1+startingwidth);
							//cout<<"mindepth = "<<mindepth<<" maxdepth = "<<maxdepth<<endl<<endl;
							double smallestbrlen=edgelength*(pow(basestep,(0-startingnumbersteps+((startingnumbersteps+1.0)/2.0))));
							midpointvector.push_back(edgelength);
							currentvector.push_back(smallestbrlen); //start with minimum values, then  move up
							incrementwidths.push_back(basestep);
							numberofedges++;
						}
						currentnode = n.next();
					}
					bool donegrid=false;
					vector<int> increments;
					increments.assign(numberofedges,0);
					int origCOALaicmode=COALaicmode;
					COALaicmode=0;
					vector<double> startingscorevector=GetCombinedScore(&NextTree);
					double currentscore=startingscorevector[0];
					while (!donegrid) {
						for (int i=0;i<currentvector.size();i++) {
							//cout<<currentvector[i]<<" ";
						}
						//cout<<endl;
						//cout<<"donegrid"<<endl;
						currentnode = n.begin();
						int nodenumber=0;
						while (currentnode)
						{
							if (currentnode!=rootnode) {
								currentnode->SetEdgeLength(currentvector[nodenumber]);
								nodenumber++;
							}
							currentnode = n.next();
						}
						vector<double> brlenscorevector=GetCombinedScore(&BrlenChangedTree);
						double brlenscore=brlenscorevector[0]-currentscore;
						//cout<<"score "<<brlenscore;
						bool atmargin=false;
						for (int i=0; i<numberofedges; i++) {
							if ((increments[i]==0) || (increments[i]==startingnumbersteps-1)) {
								atmargin=true;
								
							}
							//cout<<" "<<currentvector[i];
						}
						if (atmargin) { //if we're bumping up against minimum branchlength, there's nowhere further to expand the grid
							for (int i=0; i<numberofedges; i++) {
								if (currentvector[i]==0) {
									atmargin=false;
									
								}
							}
						}
						//cout<<endl;
						if (atmargin) { //we're at a margin of the space; want to make sure that the region within two lnL is inside this region
							if (brlenscore<2 && recursions<maxrecursions) { //our region is too small, since points 2 lnL units away from the max are outside the region
								//message="Starting width of ";
								//message+=startingwidth;
								//message+=" was too small, now increasing to ";
								startingwidth*=1.5;
								//message+=startingwidth;
								//PrintMessage();
								donegrid=true; //break the while(!donegrid) loop; since donecontour isn't done, reinitialize everything
							}
						}
						vector<double> resultvector=currentvector;
						resultvector.push_back(brlenscore);
						ContourSearchVector.push_back(resultvector);
						if (brlenscore<0 && recursions<=maxrecursions) {
							currentscore=brlenscorevector[0];
							//message="Better branch length found, restarting contour search";
							//PrintMessage();
							//recursions=0; //restart search
							NextTree.SetRoot(BrlenChangedTree.CopyOfSubtree(BrlenChangedTree.GetRoot()));
							NextTree.FindAndSetRoot();
							NextTree.Update();
							donegrid=true; //break the while(!donegrid) loop; since donecontour isn't done, reinitialize everything
							if (showtries) {
								cout<<"Better branch lengths found in grid search"<<endl;
							}
							
						}
						increments[0]++;
						if (!donegrid) {
							for (int itemtoexamine=0; itemtoexamine<numberofedges; itemtoexamine++) {
								if (increments[itemtoexamine]==startingnumbersteps) {
									if (itemtoexamine<numberofedges-1) { //means there's room to the left
										increments[itemtoexamine]=0;
										increments[itemtoexamine+1]++;
									}
									else {
										donegrid=true;
										donecontour=true;
									}
								}
							}
						}
						currentvector.clear();
						for (int i=0; i<numberofedges; i++) {
							//cout<<midpointvector[i]<<" * ("<<incrementwidths[i]<<"^"<<increments[i]-startingnumbersteps+((startingnumbersteps+1)/2)<<") = ";
							double newbrlen=midpointvector[i]*(pow(incrementwidths[i],(increments[i]-startingnumbersteps+((startingnumbersteps+1)/2))));
							currentvector.push_back(newbrlen);
							//cout<<newbrlen<<" ";
						}
						//cout<<endl;
						//cout<<"donegrid is "<<donegrid<<" and donecontour is "<<donecontour<<endl;
					}
					
					COALaicmode=origCOALaicmode;
					
				}
				vector<double> totalbrlen;
				totalbrlen.assign(numberofedges-1,0);
				int numequaltrees=0;
				for (int i=0; i<ContourSearchVector.size(); i++) {
					//if (logf_open) {
					//	for (int k=0; k<=numberofedges; k++) {
							//logf<<ContourSearchVector[i][k]<<"\t";
					//		cout<<ContourSearchVector[i][k]<<"\t";
					//	}
						//logf<<endl<<endl; 
					//cout<<endl;
					//}
					if (ContourSearchVector[i][numberofedges]<=0) {
						for (int j=0; j<numberofedges; j++) {
							totalbrlen[j]+=ContourSearchVector[i][j];
							numequaltrees++;
						}
					}
				}
				NodeIterator <Node> n (NextTree.GetRoot());
				NodePtr currentnode = n.begin();
				NodePtr rootnode=NextTree.GetRoot();
				int edgenumber=0;
				while (currentnode)
				{
					if (currentnode!=rootnode) {
						//cout<<"new brlen = "<<totalbrlen[edgenumber]/(numequaltrees*1.0)<<endl;
						currentnode->SetEdgeLength(totalbrlen[edgenumber]/(numequaltrees*1.0));
						edgenumber++;
					}
					currentnode = n.next();
				}
				nextscorevector=GetCombinedScore(&NextTree);
				nextscore=nextscorevector[0];
			}
			
			
			//cout<<"NextTree Health"<<endl;
			NextTree.FindAndSetRoot();
			NextTree.Update();
			//NextTree.ReportTreeHealth();
            if ((nextscore<bestscore) && (1==gsl_finite(nextscore))) {
                scoretype="*G\t";
                modifiedscoretype=true;
                improvement=true;
                RawBestTrees.clear();
                RawBestTrees.push_back(NextTree);
                FormattedBestTrees.clear();
				TotalScores.clear();
				GTPScores.clear();
				StructScores.clear();
				BestConversions.clear();
			//	TotalScores.push_back(nextscorevector[0]);
			//	GTPScores.push_back(nextscorevector[1]);
			//	StructScores.push_back(nextscorevector[2]);				
                FormatAndStoreBestTree(&NextTree,nextscorevector);
                //  (FormattedBestTrees.back()).Update();
                //  (FormattedBestTrees.back()).GetNodeDepths();
                //  (FormattedBestTrees.back()).Draw(cout);
                bestscore=nextscore;
                if (showtries) {
                    cout<<"GOT BETTER TREE with score of "<<nextscore<<endl<<endl;
                    NextTree.Draw(cout);
                }
            }
            else if ((nextscore==bestscore) && (1==gsl_finite(nextscore))) {
                scoretype="=G\t";
                modifiedscoretype=true;
                NextTree.Update();
                RawBestTrees.push_back(NextTree);
				//TotalScores.push_back(nextscorevector[0]);
				//GTPScores.push_back(nextscorevector[1]);
				//StructScores.push_back(nextscorevector[2]);				
                FormatAndStoreBestTree(&NextTree,nextscorevector);
                //    (FormattedBestTrees.back()).Update();
                //    (FormattedBestTrees.back()).GetNodeDepths();
                //    (FormattedBestTrees.back()).Draw(cout);
            }
			
            if ((nextscore<bestscorelocal) && (1==gsl_finite(nextscore))) {
                if (!modifiedscoretype) {
                    scoretype="*L\t";
                }
                bestscorelocal=nextscore;
                improvement=true;
                BestTreesThisRep.clear();
                NextTree.FindAndSetRoot();
                NextTree.Update();
                BestTreesThisRep.push_back(NextTree);
                CurrentTree.ResetBreakVector(); ///////figure out when to  reset this: any time you move to a new optimum
            }
            else if ((nextscore==bestscorelocal) && (1==gsl_finite(nextscore))) {
                if (!modifiedscoretype) {
                    scoretype="=L\t";
                }
                NextTree.FindAndSetRoot();
                NextTree.Update();
                BestTreesThisRep.push_back(NextTree);
                //cout<<"Swapping back the convertsamplestospecies vector\n";
                //convertsamplestospecies.swap(Originalconvertsamplestospecies); //need to reassign the original one
            }
            else {
                //cout<<"Swapping back the convertsamplestospecies vector\n";
                // convertsamplestospecies.swap(Originalconvertsamplestospecies); //need to reassign the original one
            }
            message="";
			if (jackknifesearch) {
				message+=jackrep;
				message+="\t";
			}
            message+=replicate;
            message+="\t";
            message+=movecount;
            message+="\t";
            if (chosenmove==3 || chosenmove==4) {
				message+=CurrentTree.GetNumLeaves();
				message+="->";
				message+=NextTree.GetNumLeaves();
            }
            else {
            	message+=NextTree.GetNumLeaves();
            }
            message+="\t";
            message+=chosenmovestring;
            message+="\t";
            message+=scoretype;
			if (gtptoohigh || triplettoohigh ) {
				message+=">";
			}			
			if (infinitescore) {
				message+="~";
			}
            char outputstring[9];
            sprintf(outputstring,"%9.3f",nextscore);
            message+=outputstring;
			if (!useCOAL && !useMS) {
				message+="\t";
				if (gtptoohigh) {
					message+=">";
				}						
				sprintf(outputstring,"%9.3f",nextscorevector[1]);
				message+=outputstring;
				message+="\t";
				if (gtptoohigh || triplettoohigh) { //since we abort gtp calculations if the triplet cost is already too high
					message+=">";
				}									
				sprintf(outputstring,"%9.3f",nextscorevector[2]);
				message+=outputstring;
			}
            message+="\t";
            sprintf(outputstring,"%9.3f",bestscorelocal);
            message+=outputstring;
            message+="\t";
            sprintf(outputstring,"%9.3f",GSL_MIN(bestscore,bestscorelocal));
            message+=outputstring;
            message+="\t";
            message+=int(FormattedBestTrees.size());
            if (moreswaps) {
                message+="\ts";
            }
            else {
                message+="\t_";
            }
            if (morereassignments) {
                message+="a";
            }
            else {
                message+="_";
            }
            if (moreincreases) {
                message+="i";
            }
            else {
                message+="_";
            }
            if (moredecreases) {
                message+="d";
            }
            else {
                message+="_";
            }
            if (morererootings) {
                message+="r";
            }
            else {
                message+="_";
            }
			if (movefreqvector[5]>0) {
				message+="b";
			}

			
            if (badgtpcount>0) {
                message+="\t!!!";
                message+=badgtpcount;
                message+="!!!";
            }
            if (status) {
                PrintMessage();
            }
            if ((steepest==false) && improvement) {
                if (showtries) {
                    cout<<"NOW BREAKING..."<<endl;
                }
                break;
            }
            else {
				convertsamplestospecies.assign( Originalconvertsamplestospecies.begin(), Originalconvertsamplestospecies.end() );
                //convertsamplestospecies.swap(Originalconvertsamplestospecies);
            }
			//cout<<"done outputting status line"<<endl;
        } //if something happened
		//cout<<"done if something happened loop"<<endl;
    } //while (moreswaps || morereassignments || moreincreases || moredecreases )
}//while improvement
 //DelDupes();
 //message="Replicate finished, now removing duplicate trees and saving best to file";
 //PrintMessage();
 //	ofstream outtreef;
 //	outtreef.open(treefilename.c_str());
 //outtreef<<"#nexus\nbegin trees;\n";
 //outtreef<<"[heuristic search, best results after replicate "<<replicate<<"\nSearch options: \n]\n";

//	for (int i=0; i<FormattedBestTrees.size(); i++) {
//		outtreef<<"tree sptre"<<i+1<<" = [&R] ";
//		outtreef<<ReturnFinalSpeciesTree(FormattedBestTrees[i]);
//	       outtreef<<endl;
//	}
//outtreef<<"end;";
//outtreef.close();
//  message="Best trees by the end of rep ";
//    message+=replicate;
//   message+="\n";
//    PrintMessage();
//    for (int i=0; i<FormattedBestTrees.size(); i++) {
//       (FormattedBestTrees[i]).Draw(cout);
//       cout<<endl;
//}
    }//nrep
     //DelDupes();
cout<<endl<<"Best trees overall"<<endl<<endl;
for (int i=0; i<FormattedBestTrees.size(); i++) {
    (FormattedBestTrees[i]).Update();
    (FormattedBestTrees[i]).GetNodeDepths();
    (FormattedBestTrees[i]).Draw(cout);
    cout<<endl;
}
message="\n\n#nexus\nbegin trees;\n";
for (int i=0; i<FormattedBestTrees.size(); i++) {
    message+="tree sptree";
    message+=i+1;
    message+=" = ";
    if (unrooted==1) {
        message+="[&U] ";
    }
    else {
        message+="[&R] ";
    }
    message+=ReturnFinalSpeciesTree(FormattedBestTrees[i]);
    message+="\n";
}
message+="end;\n\n";
PrintMessage();
gsl_matrix_free(TaxonDistance);
gsl_matrix_free(TaxonProportDistance);
}

void BROWNIE::HandleDettmanCollapse( NexusToken& token )
{
    for(;;)
    {
        token.GetNextToken();
		
        if( token.Equals(";") ) {
            DoDettmanCollapse();
            break;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading Dettman command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
    message="Now doing Dettman collapse\n";
    PrintMessage();
}



void BROWNIE::DoDettmanCollapse()
{
	//This just collapses branches on a partially-resolved input tree (the only branches are "independent evolutionary lineages") so that each sample is in a species. Assumes input tree has all branches but "independent evolutionary lineages" (and terminals) collapsed
	//basically, every taxon should be connecting to a branch with only taxa as descendants
	ContainingTree TreeToCollapse;
	TreeToCollapse.SetRoot((intrees.GetIthTree(chosentree-1)).CopyOfSubtree((intrees.GetIthTree(chosentree-1)).GetRoot()));
	bool isokay=false;
	while (!isokay) {
		TreeToCollapse.FindAndSetRoot();
		TreeToCollapse.Update();
		//TreeToCollapse.ReportTreeHealth();
		isokay=true;
		NodeIterator <Node> d (TreeToCollapse.GetRoot());
		NodePtr detnode = d.begin();
		while (detnode)
		{
			//cout<<"detnode is "<<detnode<<endl;
			if (isokay) {
				if (detnode->IsLeaf()) {
					NodePtr ancnode=detnode->GetAnc();
			//now make sure children of ancnode are all leaves
					if (ancnode) {
						NodePtr nodetotest=ancnode->GetChild();
						while (nodetotest!=NULL && isokay) {
							//cout<<"testing "<<nodetotest;
							if (!(nodetotest->IsLeaf())) {
								isokay=false;
								TreeToCollapse.SuppressInternalNode(nodetotest);
								//cout<<" delete below"<<endl;
							}
							else {
								nodetotest=nodetotest->GetSibling();
								//cout<<" okay"<<endl;
							}
						}
					}
				}
			}
			if (isokay) {
				detnode = d.next();
			}
			else {
				break;
			}
		}
	}
	TreeToCollapse.FindAndSetRoot();
	TreeToCollapse.Update();
	//TreeToCollapse.ReportTreeHealth();
	ofstream dettf;
	nxsstring dettfile="dettman.tre";
	dettf.open(dettfile.c_str());
	dettf<<"#nexus\nbegin trees;\n"; 
	dettf<<"tree collapsed = [&R]  ";
	TreeToCollapse.Write(dettf);
	dettf<<endl<<"end;";
	dettf.close();
	message="Collapsed tree has been saved to dettman.tre";
	PrintMessage();
}


//gives total score for structure within each putative species
double BROWNIE::GetTripletScore(ContainingTree *SpeciesTreePtr) { //Note that excess structure requires four or more individuals per species
                                                                  //cut each gene tree where it crosses a species tree boundary and reroot on the node connecting to the deleted edge
                                                                  //Watch out: don't compare two trees for the same gene (weighted trees)
                                                                  //Have matrix of expected random triplet scores (3 tax vs 4 tax, etc. ) so that you can subtract this from the observed overlap
    int oldchosentree=chosentree;
    int nspecies=SpeciesTreePtr->GetNumLeaves();
    double totalscore=0; //this is right now the default score if there are no triplets
    for (int i=1; i<=nspecies; i++) {
		if (triplettoohigh) {
			break;
		}
        int nsamplesinspecies=0;
        vector<nxsstring> taxatoexclude;
        for (int j=0; j<convertsamplestospecies.size();j++) {
            if (convertsamplestospecies[j]==i) {
                nsamplesinspecies++;
            }
            else {
                taxatoexclude.push_back(taxa->GetTaxonLabel(j));
                //	cout<<"taxatoexclude "<<taxa->GetTaxonLabel(j)<<endl;
            }
        }
        if (nsamplesinspecies>=3) { //so there's actually a triplet; otherwise, score is ____default____.
                                    //cout<<"at least three species"<<endl;
            double totalweight=0;
            vector<vector <ContainingTree> > GeneTreesVector;
            vector<vector <double> > GeneTreesWeights;
            int numberofgenes=0;
            vector<int> TreesPerGene;
            vector<ContainingTree> OneGeneTreeVector;
            vector<double> OneGeneTreeWeights;
            for (int chosentreenum=0; chosentreenum<trees->GetNumTrees(); chosentreenum++) {
				if (triplettoohigh) {
					break;
				}				
				double newweight=trees->GetTreeWeight(chosentreenum);
				if (!jackknifesearch) {
										//cout<<"Tree weight is "<<newweight<<endl;
					totalweight+=newweight;
					if (totalweight>=1) { //new gene (remember, if we have several bootstrap trees from one gene, we don't want triplet scores against them. Weight of this means that we're starting a new gene
						if (OneGeneTreeVector.size()>0) {
							GeneTreesVector.push_back(OneGeneTreeVector);
							GeneTreesWeights.push_back(OneGeneTreeWeights);
							numberofgenes++;
							TreesPerGene.push_back(OneGeneTreeVector.size());
							OneGeneTreeVector.clear();
						}
						totalweight=newweight;
					}
				}
				else {
					bool usethistree=true;
					if (jackknifevector[i]==0) {
						usethistree=false;
					}
					if (usethistree) {
						if (geneidvector[chosentreenum]>numberofgenes+1) {
							if (OneGeneTreeVector.size()>0) {
								GeneTreesVector.push_back(OneGeneTreeVector);
								GeneTreesWeights.push_back(OneGeneTreeWeights);
								numberofgenes++;
								TreesPerGene.push_back(OneGeneTreeVector.size());
								OneGeneTreeVector.clear();
							}
						}
					}
				}
				bool usethistree=true;
				if (jackknifesearch) {
				if (jackknifevector[i]==0) {
					usethistree=false;
				}
				}
				if (usethistree) {
					//Now we can deal with the gene tree, having re-initialized the vector if need be
					Tree CurrentGeneTreeTreeFmt=intrees.GetIthTree(chosentreenum);
					//	cout<<"Gene tree is "<<endl;
					//	CurrentGeneTreeTreeFmt.Draw(cout);
					ContainingTree CurrentGeneTree;
					CurrentGeneTree.SetRoot(CurrentGeneTreeTreeFmt.CopyOfSubtree(CurrentGeneTreeTreeFmt.GetRoot()));
					vector<ContainingTree> SplitTreeVector;
					SplitTreeVector=CurrentGeneTree.SplitOnTaxon(taxatoexclude);
					for (int k=0; k<SplitTreeVector.size(); k++) {
						OneGeneTreeVector.push_back(SplitTreeVector[k]);
						OneGeneTreeWeights.push_back(newweight);
					}
				}
            }
            if (OneGeneTreeVector.size()>0) { //add the last gene
                GeneTreesVector.push_back(OneGeneTreeVector);
                GeneTreesWeights.push_back(OneGeneTreeWeights);
                numberofgenes++;
                TreesPerGene.push_back(OneGeneTreeVector.size());
                OneGeneTreeVector.clear();
            }
            //now compare all genes to all other genes
            for (int chosengene1=0;chosengene1<numberofgenes-1; chosengene1++) {
				if (triplettoohigh) {
					break;
				}				
                for (int chosengene2=chosengene1+1;chosengene2<numberofgenes; chosengene2++) {
					if (triplettoohigh) {
						break;
					}					
                    for (int chosentreenum1=0; chosentreenum1<TreesPerGene[chosengene1]; chosentreenum1++) {
						if (triplettoohigh) {
							break;
						}						
                        //ContainingTree Tree1=(GeneTreesVector[chosengene1][chosentreenum1]);
                        assert(GeneTreesVector.size()>chosengene1);
                        assert((GeneTreesVector[chosengene1]).size()>chosentreenum1);
                        double Tree1Wt=GeneTreesWeights[chosengene1][chosentreenum1];
                        int Tree1Ntax=(GeneTreesVector[chosengene1][chosentreenum1]).GetNumLeaves();
                        for (int chosentreenum2=0; chosentreenum2<TreesPerGene[chosengene2]; chosentreenum2++) {
							if (triplettoohigh) {
								break;
							}							
                            double Tree2Wt=GeneTreesWeights[chosengene2][chosentreenum2];
                            assert(GeneTreesVector.size()>chosengene2);
                            assert((GeneTreesVector[chosengene2]).size()>chosentreenum2);
                            //	(GeneTreesVector[chosengene2][chosentreenum2]).Write(cout);
                            //	cout<<endl;
                            (GeneTreesVector[chosengene2][chosentreenum2]).Update();
                            int Tree2Ntax=(GeneTreesVector[chosengene2][chosentreenum2]).GetNumLeaves();
                            if ((Tree1Ntax>2) && (Tree2Ntax>2)) { //so, at least a triplet in each
                                ContainingTree t1=GeneTreesVector[chosengene1][chosentreenum1]; //Need to store copies of trees because they're modified (leaves deleted)
                                ContainingTree t2=GeneTreesVector[chosengene2][chosentreenum2];
                                int taxaincommon=PrepareTreesForTriplet(&t1,&t2);
                                if (taxaincommon>2) {
									//cout<<"t1=\n";
									//t1.Draw(cout);
									//cout<<"\nt2=\n";
									//t2.Draw(cout);
                                    vector<int> tripletoverlapoutput=GetTripletOverlap(&t1,&t2,taxaincommon);
                                    int maxnumber=tripletoverlapoutput[0];
                                    int numberdisagree=tripletoverlapoutput[1];
                                    int numberunresolved=tripletoverlapoutput[2];
                                    int numberagree=maxnumber-numberdisagree-numberunresolved; //Note that this is the number of triplets resolved IN BOTH TREES that agree
                                    double newscore=ComputeTripletCost(numberagree,maxnumber,taxaincommon,Tree1Wt,Tree1Ntax,Tree2Wt,Tree2Ntax,numberofgenes);
                                    // cout<<"maxnum="<<maxnumber<<" dis="<<numberdisagree<<" un="<<numberunresolved<<" agr="<<numberagree<<" newscore="<<newscore<<endl;
                                    totalscore+=newscore;
									if ((totalscore*structwt)>bestscorelocal) {
										totalscore=(0.0001+bestscorelocal)/structwt;
										triplettoohigh=true;
										break;
									}
                                }
                            }
                        }
                    }
                }
            }
            
        }
        
    }
    chosentree=oldchosentree;
    return totalscore;
}

//computes the TaxonDistance matrix. each entry (i,j) is the number of times taxon i and taxon j are each others' closest relatives in a triplet
void BROWNIE::GetTaxonTaxonTripletDistances() {
	TripletCounts.clear(); //Clear the triplet counts map
	TripleCounts.clear(); //Clears the triples counts map (number of time each set of three taxa co-occur on a tree)
    int nsamples=taxa->GetNumTaxonLabels();
    //cout<<"nsamples = "<<nsamples<<endl;
    if (nsamples==0) {
             errormsg = "ERROR: There are no samples recorded. This can happen if the taxa block is not formatted correctly.";
            throw XNexus( errormsg);
    }
    TaxonDistance=gsl_matrix_calloc(nsamples,nsamples); //on diagonal is number of triplets containing the taxon, other elements are as above.
                                                        //	gsl_matrix *TaxonTripletsNotOnSameSide=gsl_matrix_calloc(,);
                                                        //	gsl_matrix *TaxonTripletsTotal=gsl_matrix_calloc(intrees.GetNumTrees(),intrees.GetNumTrees()); //at end, do element by element division, watch for division by zero
    for (int i=0;i<intrees.GetNumTrees();i++) {
			message="now getting triplets on tree ";
			message+=i+1;
			message+=" of ";
			message+=intrees.GetNumTrees();
			PrintMessage();
		bool usethistree=true;
		if (jackknifesearch) {
			if (jackknifevector[i]==0) {
				usethistree=false;
			}
		}
		if (usethistree) {
			Tree t1=intrees.GetIthTree(i);
			gsl_combination *combo;
			combo=gsl_combination_calloc(t1.GetNumLeaves(),3);
			vector<nxsstring> LeafLabelVect; //this section taken from containingtree.cpp (no easy way to convert tree to containing tree)
			NodeIterator <Node> n (t1.GetRoot());
			NodePtr currentnode = n.begin();
			t1.MakeNodeList(); //moved from below
			while (currentnode)
			{
				if(currentnode->IsLeaf()) {
					LeafLabelVect.push_back(currentnode->GetLabel());
				}
				currentnode = n.next();
			}
			do
			{
				NodePtr p;
				NodePtr q;
				nxsstring a=LeafLabelVect[gsl_combination_get(combo,0)];
				nxsstring b=LeafLabelVect[gsl_combination_get(combo,1)];
				nxsstring c=LeafLabelVect[gsl_combination_get(combo,2)];
				std::map<Node *, int, std::less<Node *> > depth;
				PreorderIterator <Node> n (t1.GetRoot());
				int count = 0;
				q = n.begin();
				while (q)
				{
					depth[q] = count++;
					q = n.next();
				}
				
            //Now get depths for each pair
				vector<int> LCADepthVectorT1;
				//t1.MakeNodeList(); //Move this up a level -- calculate it once, see if it goes faster
				p = t1.GetLeafWithLabel(a.c_str());
				q = t1.GetLeafWithLabel(b.c_str());
				while (depth[p] != depth[q])
				{
					if (depth[p] < depth[q])
						q = q->GetAnc();
					else
						p = p->GetAnc();
				}
				LCADepthVectorT1.push_back(depth[p]);
				p = t1.GetLeafWithLabel(b.c_str());
				q = t1.GetLeafWithLabel(c.c_str());
				while (depth[p] != depth[q])
				{
					if (depth[p] < depth[q])
						q = q->GetAnc();
					else
						p = p->GetAnc();
				}
				LCADepthVectorT1.push_back(depth[p]);
				p = t1.GetLeafWithLabel(a.c_str());
				q = t1.GetLeafWithLabel(c.c_str());
				while (depth[p] != depth[q])
				{
					if (depth[p] < depth[q])
						q = q->GetAnc();
					else
						p = p->GetAnc();
				}
				LCADepthVectorT1.push_back(depth[p]);
				int T1DepthMax=0; //root has depth 0, others have higher depths.
				int T1DepthMaxIndex=0;
				int anum=taxa->FindTaxon(a);
				int bnum=taxa->FindTaxon(b);
				int cnum=taxa->FindTaxon(c);
            //cout<<anum<<": "<<a<<" "<<bnum<<": "<<b<<" "<<cnum<<": "<<c<<endl;
				for (int j=0;j<3;j++) {
					if(LCADepthVectorT1[j]>T1DepthMax) {
						T1DepthMax=LCADepthVectorT1[j];
						T1DepthMaxIndex=j;
					}
				}
            //j=0->ab, j=1->bc j=2->ac
				gsl_matrix_set(TaxonDistance,anum,anum,1+gsl_matrix_get(TaxonDistance,anum,anum));
				gsl_matrix_set(TaxonDistance,bnum,bnum,1+gsl_matrix_get(TaxonDistance,bnum,bnum));
				gsl_matrix_set(TaxonDistance,cnum,cnum,1+gsl_matrix_get(TaxonDistance,cnum,cnum));
				nxsstring tripletlabel="";
				if (T1DepthMaxIndex==0) {
					gsl_matrix_set(TaxonDistance,anum,bnum,1+gsl_matrix_get(TaxonDistance,anum,bnum));
					gsl_matrix_set(TaxonDistance,bnum,anum,gsl_matrix_get(TaxonDistance,anum,bnum));
					tripletlabel+=GSL_MIN(anum,bnum);
					tripletlabel+="_";
					tripletlabel+=GSL_MAX(anum,bnum);
					tripletlabel+="_";
					tripletlabel+=cnum;
					TripletCounts[tripletlabel]++;
				}
				else if (T1DepthMaxIndex==1) {
					gsl_matrix_set(TaxonDistance,bnum,cnum,1+gsl_matrix_get(TaxonDistance,bnum,cnum));
					gsl_matrix_set(TaxonDistance,cnum,bnum,gsl_matrix_get(TaxonDistance,bnum,cnum));
					tripletlabel+=GSL_MIN(cnum,bnum);
					tripletlabel+="_";
					tripletlabel+=GSL_MAX(cnum,bnum);
					tripletlabel+="_";
					tripletlabel+=anum;
					TripletCounts[tripletlabel]++;
				}
				else if (T1DepthMaxIndex==2) {
					gsl_matrix_set(TaxonDistance,anum,cnum,1+gsl_matrix_get(TaxonDistance,anum,cnum));
					gsl_matrix_set(TaxonDistance,cnum,anum,gsl_matrix_get(TaxonDistance,anum,cnum));
					tripletlabel+=GSL_MIN(cnum,anum);
					tripletlabel+="_";
					tripletlabel+=GSL_MAX(cnum,anum);
					tripletlabel+="_";
					tripletlabel+=bnum;
					TripletCounts[tripletlabel]++;
				}
				int taxonarray[3]={anum, bnum, cnum};
				sort(taxonarray,taxonarray+3);
				nxsstring triplelabel="";
				triplelabel+=taxonarray[0];
				triplelabel+="_";
				triplelabel+=taxonarray[1];
				triplelabel+="_";
				triplelabel+=taxonarray[2];
				TripleCounts[triplelabel]++;
				//if (debugmode) {
				//	cout<<"just got triplet for taxa "<<anum<<", "<<bnum<<", and "<<cnum<<endl;
				//}
			}
			while (gsl_combination_next(combo) == GSL_SUCCESS);
			gsl_combination_free(combo);
		}
    }
    TaxonProportDistance=gsl_matrix_calloc(nsamples,nsamples);
	double maxnumbertriplets=0;
for (int i=0;i<nsamples;i++) {
	for (int j=0;j<nsamples;j++) {
		if (i!=j) {
			maxnumbertriplets=GSL_MAX(gsl_matrix_get(TaxonDistance,i,j),maxnumbertriplets);
		}
	}
}
maxnumbertriplets++; //This is to ensure that taxa are closer to themselves than to other taxa
    //double maxnumbertriplets=gsl_matrix_max(TaxonDistance);
if (debugmode) {
	cout<<"Original similarity matrix\n";
}
for (int i=0;i<nsamples;i++) {
	if (debugmode) {
		cout<<endl<<i<<"\t";
	}
	for (int j=0;j<nsamples;j++) {
		if (debugmode) {
			cout<<gsl_matrix_get(TaxonDistance,i,j)<<"\t";
		}
	}
}
if (debugmode) {
	cout<<endl<<"\nNew distance matrix\n";
}



    for (int i=0;i<nsamples;i++) {
		if (debugmode) {
				cout<<endl<<i<<"\t";
		}
        for (int j=0;j<nsamples;j++) {
			if (i==j) {
				gsl_matrix_set(TaxonProportDistance,i,j,0.0);
			}
			else {
				gsl_matrix_set(TaxonProportDistance,i,j,1-((1.0*gsl_matrix_get(TaxonDistance,i,j))/maxnumbertriplets));
			}
            //gsl_matrix_set(TaxonProportDistance,i,j,1-((1.0*gsl_matrix_get(TaxonDistance,i,j))/maxnumbertriplets));
			if (debugmode) {
					cout<<gsl_matrix_get(TaxonProportDistance,i,j)<<"\t";
			}
        }
    }
	if (debugmode) {
		cout<<endl;
	}
    //for (int i=0;i<nsamples;i++) {
    //		cout<<endl;
    //		cout<<i<<": "<<taxa->GetTaxonLabel(i)<<"\t";
    //		for (int j=0;j<nsamples;j++) {
    //			cout<<"\t"<<gsl_matrix_get(TaxonProportDistance,i,j);
    //		}
    //	}
    //for (int i=0;i<nsamples;i++) {
    //	cout<<endl;
    //	cout<<i<<": "<<taxa->GetTaxonLabel(i)<<"\t";
    //	for (int j=0;j<nsamples;j++) {
    //		cout<<"\t"<<gsl_matrix_get(TaxonDistance,i,j);
    //	}
    //}
    //	cout<<"\nTaxondistance dimensions = "<<TaxonDistance->size1<<","<<TaxonDistance->size2<<endl;
    //	cout<<endl<<"last element = ";
    //	cout<<gsl_matrix_get(TaxonDistance,-1+(TaxonDistance->size1),-1+(TaxonDistance->size1));
}

ContainingTree BROWNIE::ComputeTripletNJTree() {
    ContainingTree NJresult;
    vector<NodePtr> NewNodeVector;
    vector<bool> UnattachedNodeVector;
    vector<int> MatrixToNodeVector;
    int nsamples=taxa->GetNumTaxonLabels();
    int nsamplesremaining=nsamples;
    gsl_matrix *dmatrix=gsl_matrix_calloc(nsamples,nsamples);
    for (int i=0;i<nsamples;i++) {
        NodePtr p=NJresult.NewNode();
        p->SetLabel(taxa->GetTaxonLabel(i));
        p->SetLeaf(true);
        p->SetLeafNumber(i+1);
        p->SetWeight(0);
        NJresult.SetEdgeLengths(true);
        NJresult.SetRooted(true);
        NewNodeVector.push_back(p);
        UnattachedNodeVector.push_back(true); //true means unattached
        MatrixToNodeVector.push_back(i);
        for (int j=0;j<nsamples;j++) {
            gsl_matrix_set(dmatrix,i,j,gsl_matrix_get(TaxonProportDistance,i,j));
        }
    }
    while (nsamplesremaining>2) {
        //cout<<"\n\nMatrixToNodeVector=";
        ////for (int i=0;i<nsamplesremaining;i++) {
        //	cout<<" "<<MatrixToNodeVector[i];
        //}
        //cout<<"\n\nNewNodeVector=";
        //for (int i=0;i<NewNodeVector.size();i++) {
        //	cout<<" "<<NewNodeVector[i]<<"("<<(NewNodeVector[i])->GetLabel()<<")";
        //}
        //for (int i=0;i<nsamplesremaining;i++) {
        //	cout<<endl<<i<<": ";
        //	for (int j=0;j<nsamplesremaining;j++) {
        //		cout<<"\t"<<gsl_matrix_get(dmatrix,i,j);
        //	}
        //}
        //cout<<endl;
gsl_vector *rvector=gsl_vector_calloc(nsamplesremaining);
for (int i=0;i<nsamplesremaining;i++) {
    for (int j=0;j<nsamplesremaining;j++) {
        if(i!=j) {
            gsl_vector_set(rvector,i,(gsl_vector_get(rvector,i))+gsl_matrix_get(dmatrix,i,j));
        }
    }
}
int imin;
int jmin;
double Mijmin=GSL_POSINF;
for (int i=0;i<nsamplesremaining;i++) {
    for (int j=i+1;j<nsamplesremaining;j++) {
        double Mij=gsl_matrix_get(dmatrix,i,j)-((gsl_vector_get(rvector,i)+gsl_vector_get(rvector,j))/(nsamplesremaining-2.0));
        if (Mij<Mijmin) {
            Mijmin=Mij;
            imin=i;
            jmin=j;
        }
    }
}
NodePtr u=NJresult.NewNode();
u->SetWeight(0);
NodePtr i=NewNodeVector[MatrixToNodeVector[imin]];
NodePtr j=NewNodeVector[MatrixToNodeVector[jmin]];
//cout<<"Now joining imin: "<<imin<<": "<<i<<"("<<i->GetLabel()<<") and jmin: "<<jmin<<": "<<j<<"("<<j->GetLabel()<<") to make "<<u<<endl;
u->SetLeaf(false);
u->SetChild(i);
i->SetAnc(u);
i->SetSibling(j);
j->SetAnc(u);
gsl_matrix *dmatrixnew=gsl_matrix_calloc(nsamplesremaining-1,nsamplesremaining-1);
gsl_matrix *dmatrixorig=gsl_matrix_calloc(nsamplesremaining,nsamplesremaining);
gsl_matrix_memcpy(dmatrixorig,dmatrix);
i->SetEdgeLength(gsl_matrix_get(dmatrix,imin,jmin)/2.0+(1.0*gsl_vector_get(rvector,imin)-1.0*gsl_vector_get(rvector,jmin))/(2.0*(1.0*nsamplesremaining-2.0)));
j->SetEdgeLength(gsl_matrix_get(dmatrix,imin,jmin)-(i->GetEdgeLength()));
NewNodeVector.push_back(u);
UnattachedNodeVector.push_back(true);
UnattachedNodeVector[MatrixToNodeVector[imin]]=false;
UnattachedNodeVector[MatrixToNodeVector[jmin]]=false;
bool seenj=false;
for (int i=0;i<MatrixToNodeVector.size();i++) {
    if (seenj==false) {
        if (i==jmin) {
            seenj=true;
        }
        else if (i==imin) {
            MatrixToNodeVector[i]=-1+NewNodeVector.size();
        }
    }
    else {
        if (i==imin) {
            MatrixToNodeVector[i-1]=-1+NewNodeVector.size();
        }
        else {
            MatrixToNodeVector[i-1]=MatrixToNodeVector[i];
        }
    }
}
MatrixToNodeVector.pop_back();
for (int i=0;i<nsamplesremaining;i++) {
    for (int j=0;j<nsamplesremaining;j++) {
        if(i==imin) {
            double newvalue=0.0;
            newvalue+=gsl_matrix_get(dmatrixorig,j,imin);
            newvalue+=gsl_matrix_get(dmatrixorig,j,jmin);
            newvalue-=gsl_matrix_get(dmatrixorig,imin,jmin);
            newvalue*=0.5;
            //			cout<<"changed pos ("<<i<<","<<j<<") from value "<<gsl_matrix_get(dmatrix,i,j)<<" or "<<gsl_matrix_get(dmatrixorig,i,j);
            gsl_matrix_set(dmatrix,i,j,newvalue);
            //			cout<<" to have value "<<gsl_matrix_get(dmatrixorig,j,imin)<<"+"<<gsl_matrix_get(dmatrixorig,j,jmin)<<"-"<<gsl_matrix_get(dmatrixorig,imin,jmin)<<"="<<gsl_matrix_get(dmatrix,i,j)<<endl;
        }
        if(j==imin) {
            double newvalue=0.0;
            newvalue+=gsl_matrix_get(dmatrixorig,imin,i);
            newvalue+=gsl_matrix_get(dmatrixorig,jmin,i);
            newvalue-=gsl_matrix_get(dmatrixorig,imin,jmin);
            newvalue*=0.5;
            //			cout<<"changed pos ("<<i<<","<<j<<") from value "<<gsl_matrix_get(dmatrix,i,j)<<" or "<<gsl_matrix_get(dmatrixorig,i,j);
            gsl_matrix_set(dmatrix,i,j,newvalue);
            //			cout<<" to have value "<<gsl_matrix_get(dmatrixorig,imin,i)<<"+"<<gsl_matrix_get(dmatrixorig,jmin,i)<<"-"<<gsl_matrix_get(dmatrixorig,imin,jmin)<<"="<<gsl_matrix_get(dmatrix,i,j)<<endl;
        }
    }
}
gsl_matrix_set(dmatrix,imin,imin,0);
gsl_matrix_set(dmatrix,jmin,jmin,0);
gsl_matrix_set(dmatrix,imin,jmin,0);
gsl_matrix_set(dmatrix,jmin,imin,0);
int rowoffset=0;
gsl_matrix *dmatrixintermed=gsl_matrix_calloc(nsamplesremaining-1,nsamplesremaining);
for (int i=0;i<nsamplesremaining;i++) {
    if (i==jmin) {
        rowoffset=1;
    }
    else {
        for (int j=0;j<nsamplesremaining;j++) {
            gsl_matrix_set(dmatrixintermed,i-rowoffset,j,gsl_matrix_get(dmatrix,i,j));
        }
    }
}
//cout<<"\ndintermed\n";
//for (int i=0;i<nsamplesremaining-1;i++) {
//	cout<<endl<<i<<":\t";
//	for (int j=0;j<nsamplesremaining;j++) {
//		cout<<"\t"<<gsl_matrix_get(dmatrixintermed,i,j);
//	}
//}
for (int i=0;i<nsamplesremaining-1;i++) {
    int coloffset=0;
    for (int j=0;j<nsamplesremaining;j++) {
        if (j==jmin) {
            coloffset=1;
        }
        else {
            gsl_matrix_set(dmatrixnew,i,j-coloffset,gsl_matrix_get(dmatrixintermed,i,j));
        }
    }
}

nsamplesremaining--;
gsl_vector_free(rvector);
gsl_matrix_free(dmatrix);
dmatrix=gsl_matrix_calloc(nsamplesremaining,nsamplesremaining);
gsl_matrix_memcpy(dmatrix,dmatrixnew);
gsl_matrix_free(dmatrixintermed);
gsl_matrix_free(dmatrixnew);
gsl_matrix_free(dmatrixorig);
    }
//now two nodes left to attach
NodePtr u=NJresult.NewNode();
u->SetWeight(0);
vector<NodePtr> NodesToConnect;
for(int i=0;i<UnattachedNodeVector.size();i++) {
    if(UnattachedNodeVector[i]) {
        NodesToConnect.push_back(NewNodeVector[i]);
    }
}
assert(NodesToConnect.size()==2);
NodePtr i=NodesToConnect[0];
NodePtr j=NodesToConnect[1];
u->SetLeaf(false);
u->SetChild(i);
i->SetAnc(u);
i->SetSibling(j);
j->SetAnc(u);
i->SetEdgeLength(0.5*gsl_matrix_get(dmatrix,0,0));
j->SetEdgeLength(0.5*gsl_matrix_get(dmatrix,0,0));
for (int i=0;i<NewNodeVector.size();i++) {
    NodePtr p=NewNodeVector[i];
    //cout<<"p = "<<p;
    //cout<<" "<<p->GetLabel();
    //cout<<" Anc "<<p->GetAnc()<<endl;
}
NJresult.SetRoot(u);
//cout<<"Root is "<<NJresult.GetRoot()<<endl;
//NJresult.ReportTreeHealth();
//NJresult.GetNodeDepths();
NJresult.Update();
//	NJresult.Draw(cout);
//	cout<<endl;
//	NJresult.Write(cout);
//	NJresult.SetPathLengths();
//	double minpathlength=GSL_POSINF;
//	double maxpathlength=GSL_NEGINF;
//	NodePtr MaxLengthNode;
//	NodePtr MinLengthNode;
//	NodeIterator <Node> n (NJresult.GetRoot());
//	NodePtr currentnode = n.begin();
//	while (currentnode) {
//		if (currentnode->IsLeaf()) {
//			if (currentnode->GetPathLength()<minpathlength) {
//				minpathlength=currentnode->GetPathLength();
//				MinLengthNode=currentnode;
//			}
//			else if (currentnode->GetPathLength()>maxpathlength) {
//				maxpathlength=currentnode->GetPathLength();
//				MaxLengthNode=currentnode;
//			}
//		}
//		currentnode = n.next();
//	}
//	//now, midpoint root the tree
//	double desiredpathlength=(maxpathlength+minpathlength)/2.0;
//	currentnode=MaxLengthNode;
//	while (currentnode->GetPathLength()>desiredpathlength) {
//		if (((currentnode->GetAnc())->GetPathLength())<=desiredpathlength) {
//			NodePtr newroot=NJresult.ReRootTree(currentnode);
//			break;
//		}
//		currentnode=currentnode->GetAnc();
//	}
vector<int> origconvertsamplestospecies;
origconvertsamplestospecies.swap(convertsamplestospecies);
convertsamplestospecies.clear();
for (int i=0;i<nsamples;i++) {
    convertsamplestospecies.push_back(i+1);
}
ContainingTree NJrenamed=NJresult;
nxsstring NewLeafLabel;
NodeIterator <Node> n (NJrenamed.GetRoot());
NodePtr currentnode = n.begin();
while (currentnode)
{
    if (currentnode->IsLeaf())
    {
        NewLeafLabel="taxon";
        NewLeafLabel+=currentnode->GetLeafNumber();
        currentnode->SetLabel(NewLeafLabel);
    }
    currentnode = n.next();
}
//cout<<"Renamed tree"<<endl;
NJrenamed.Update();
//NJrenamed.Draw(cout);
//cout<<endl;
//cout<<"Now geting initial GTP score, it =";
NodePtr OldRoot=(NJrenamed.GetRoot())->GetChild(); //since the previous root itself will be deleted, need to reroot on its descendant.
NodePtr tmp=NJrenamed.ReRootTree(OldRoot); //so that the child its sib don't change when we reset
NJrenamed.Update();
double bestGTPscore=GetGTPScoreNew(&NJrenamed);
//cout<<bestGTPscore;
//cout<<endl;
int NodeToReRootOn=0;
//cout<<"Now starting nodeIterator m\n";
NodeIterator <Node> m (NJrenamed.GetRoot());
//cout<<"NOw setting currentnode=m.begin()\n";
NodePtr mcurrentnode = m.begin();
//cout<<"Now in while(currentnode)\n";
vector<NodePtr> nodestorerootNJrenamedon;
while (mcurrentnode) {
    if (mcurrentnode!=(NJrenamed.GetRoot())) {
        nodestorerootNJrenamedon.push_back(mcurrentnode);
    }
    mcurrentnode=m.next();
}
for (int i=0;i<nodestorerootNJrenamedon.size();i++) {
    NodePtr NewRoot=NJrenamed.ReRootTree(nodestorerootNJrenamedon[i]);
    NJrenamed.Update();
    double newGTPscore=GetGTPScoreNew(&NJrenamed);
    if (newGTPscore<bestGTPscore) {
        NodeToReRootOn=i+1;
        bestGTPscore=newGTPscore;
    }
    //  cout<<"FOr rerooting "<<i+1<<endl;
    //  NJrenamed.Draw(cout);
    // cout<<"\nscore is "<<newGTPscore<<endl;
    NewRoot=NJrenamed.ReRootTree(OldRoot); //Reset it to initial tree
}

//cout<<"OPtimal rerooting on node "<<NodeToReRootOn<<endl;

int currentnodenum=0;
convertsamplestospecies.swap(origconvertsamplestospecies);
if (NodeToReRootOn!=0) {
    NodeIterator <Node> o (NJresult.GetRoot());
    currentnode=o.begin();
    while (currentnode) {
        if (currentnode!=(NJresult.GetRoot())) {
            currentnodenum++;
            if (currentnodenum==NodeToReRootOn) {
                NodePtr NewRoot=NJresult.ReRootTree(currentnode);
            }
        }
        currentnode=o.next();
    }
}



NJresult.Update();
message="Tree of samples: ";
PrintMessage();
NJresult.Draw(cout);
message="";
PrintMessage();

//now modify CladeVector;
CladeVector.clear();
CladeVectorNJBrlen.clear(); //This is used to record internal brlen for the triplet tree, which is used to help do the initial assignments
CladeVectorTripletSupport.clear(); //stores the mean internal brlen of the starting nj tree, used for deciding on initial splits
for (int i=0;i<(2*nsamples)-2;i++) {
    vector<int> tempvector;
    CladeVector.push_back(tempvector);
	CladeVectorNJBrlen.push_back(0.0); //Just initializing
	CladeVectorTripletSupport.push_back(0.0); //Just initializing
}
if (debugmode) {
NJresult.ReportTreeHealth();
}
NodeIterator <Node> l (NJresult.GetRoot()); //goes from tips down
currentnode=l.begin();
int nodecount=0;
double brlensum=0.0;
while (currentnode) {
    vector<int> tempvector; //stores the list of terminal taxa as we go down the tree
    if (currentnode->IsLeaf()) {
        int leafnumtrans=(currentnode->GetLeafNumber())-1;
		if (debugmode) {
			cout<<"looking at leaf with node number "<<leafnumtrans<<" and address "<<currentnode<<endl;
		}
        currentnode->SetIndex(leafnumtrans);
        tempvector.push_back(leafnumtrans);
        CladeVector[leafnumtrans]=tempvector;
		CladeVectorNJBrlen[leafnumtrans]=-1.0*(currentnode->GetEdgeLength()); //terminal branches are stored as negative lengths so we can distinguish them later
    }
    else if(currentnode!=NJresult.GetRoot()) {
        int nodenumber=nodecount+NJresult.GetNumLeaves();
		if (debugmode) {
			cout<<"looking at internal with node number "<<nodenumber<<" and address "<<currentnode<<endl;
		}		
        currentnode->SetIndex(nodenumber);
        int desc1=(currentnode->GetChild())->GetIndex();
        int desc2=((currentnode->GetChild())->GetSibling())->GetIndex();
        for (int i=0;i<(CladeVector[desc1]).size();i++) {
            tempvector.push_back(CladeVector[desc1][i]);
        }
        for (int i=0;i<(CladeVector[desc2]).size();i++) {
            tempvector.push_back(CladeVector[desc2][i]);
        }
        CladeVector[nodenumber]=tempvector;
		CladeVectorNJBrlen[nodenumber]=currentnode->GetEdgeLength();
		brlensum+=currentnode->GetEdgeLength();
		

		
        nodecount++;
    }
	
    currentnode=l.next();
}
meaninternalbrlen=brlensum/(1.0*nodecount);
if (debugmode) {
	cout<<endl<<"average = "<<meaninternalbrlen<<endl;
}

//Get triplet support on edges: have to loop through the tree again (first had to assign clade vectors to tree)
NodeIterator <Node> q (NJresult.GetRoot()); //goes from tips down
currentnode=q.begin();
nodecount=0;
while (currentnode) {
    if (currentnode->IsLeaf()) {
        int leafnumtrans=(currentnode->GetLeafNumber())-1;
		if (debugmode) {
			cout<<"looking at leaf with node number "<<leafnumtrans<<" and address "<<currentnode<<endl;
		}
		CladeVectorTripletSupport[leafnumtrans]=0.0; //triplet support for an OTU is meaningless
    }
    else if(currentnode!=NJresult.GetRoot()) {
        int nodenumber=nodecount+NJresult.GetNumLeaves();
		if (debugmode) {
			cout<<"looking at internal with node number "<<nodenumber<<" and address "<<currentnode<<endl;
		}		
		
		//Get triplet support on edges
		/*
		 If this is the tree, and we want to get support for the edge below the MRCA of taxon p1samp2 and p1samp3, 
		 we need the proportion of relevant triplets with (samp2 or samp6) forming a clade with samp3 relative to (samp1 or samp7)
		 (we want to know the support ONLY at that edge, not using triplets relevant for branches higher or lower in the tree).
		 
		 +- p1samp1
							  +--|
		 |  +- p1samp7
		 +-|
		 | |  +- p1samp2
		 | | +|
		 | | |+- p1samp6
		 | +-|
		 |   +-- p1samp3
		 ---------|
		 |    +- p1samp4
		 +----|
		 +- p1samp5
		 
		 basically, for node with two descendant clades A and B and with sister node (clade) C, want the proportion of triplets matching (A,B),C, where the 
		 taxa are pulled from those clades
		 */		 
		
		//Remember, this tree must be binary (thank goodness!)
		vector<int> avector=CladeVector[(currentnode->GetChild())->GetIndex()]; //Get the child of the current node, get its node number, get the list of taxa at that node
		vector<int> bvector=CladeVector[((currentnode->GetChild())->GetSibling())->GetIndex()]; //Get the sibling of the child of the current node, get its node number, get the list of taxa at that node
		vector<int> cvector;
		if ((currentnode->GetSibling())!=NULL) {
			cvector=CladeVector[(currentnode->GetSibling())->GetIndex()]; //get list of taxa in sister node
			if (debugmode) {
				cout<<"Sister node = "<<currentnode->GetSibling()<<endl;
			}
			
		}
		else { //This node is the sibling, so have to go down one step and back up to get its actual sib
			cvector=CladeVector[((currentnode->GetAnc())->GetChild())->GetIndex()];
			if (debugmode) {
				cout<<"Sister node = "<<(currentnode->GetAnc())->GetChild()<<endl;
			}
		}
		int numberofmatchingtriplets=0;
		int numberoftriplespossible=0;
		if (debugmode) {
			cout<<"vector sizes (a, b, c): "<<avector.size()<<" "<<bvector.size()<<" "<<cvector.size()<<endl;
		}
		for (int aindex=0;aindex<avector.size();aindex++) {
			for (int bindex=0;bindex<bvector.size();bindex++) {
				nxsstring tripletstringroot="";
				int anum=avector[aindex];
				int bnum=bvector[bindex];
				tripletstringroot+=GSL_MIN(anum,bnum);
				tripletstringroot+="_";
				tripletstringroot+=GSL_MAX(anum,bnum);
				tripletstringroot+="_";
				for (int cindex=0;cindex<cvector.size();cindex++) {
					nxsstring tripletstring=tripletstringroot;
					int cnum=cvector[cindex];
					tripletstring+=cnum;
					int taxonarray[3]={anum, bnum, cnum};
					sort(taxonarray,taxonarray+3);
					nxsstring triplelabel="";
					triplelabel+=taxonarray[0];
					triplelabel+="_";
					triplelabel+=taxonarray[1];
					triplelabel+="_";
					triplelabel+=taxonarray[2];
					map<nxsstring,int>::iterator iter = TripletCounts.find(tripletstring);
					if( iter != TripletCounts.end() ) {
						numberofmatchingtriplets+=iter->second;
					}
					map<nxsstring,int>::iterator iter2 = TripleCounts.find(triplelabel);
					if( iter != TripleCounts.end() ) {
						numberoftriplespossible+=iter2->second;
					}					
					if (debugmode) {
						cout<<"Triplet = "<<tripletstring<<" triple = "<<triplelabel<<endl;
					}
				}
			}
		}
		if (debugmode) {
			cout<<"number of matching triplets = "<<numberofmatchingtriplets<<" number of triplets possible = "<<numberoftriplespossible;
		}
		CladeVectorTripletSupport[nodenumber]=(1.0*numberofmatchingtriplets)/(1.0*numberoftriplespossible);
		assert (CladeVectorTripletSupport[nodenumber]<=1 && CladeVectorTripletSupport[nodenumber]>=0);
		if (debugmode) {
			cout<<" proport = "<<CladeVectorTripletSupport[nodenumber]<<endl;
		}
        nodecount++;
    }
	
    currentnode=q.next();
}






if (debugmode) {
	
	for( map<nxsstring, int>::iterator iter = TripletCounts.begin(); iter != TripletCounts.end(); iter++ ) {
		cout << "Triplet "<<(*iter).first << " found " << (*iter).second << " times" << endl;
	}
	
	cout<<endl<<endl<<endl;
	for( map<nxsstring, int>::iterator iter = TripleCounts.begin(); iter != TripleCounts.end(); iter++ ) {
		cout << "Triple "<<(*iter).first << " found " << (*iter).second << " times" << endl;
	}
}
//for (int i=0;i<CladeVector.size();i++) {
//    cout<<endl;
//    for (int j=0;j<(CladeVector[i]).size();j++) {
//        cout<<" "<<CladeVector[i][j];
//   }
//}
//cout<<endl;
gsl_matrix_free(dmatrix);
return NJresult;
}


vector<nxsstring> BROWNIE::ReturnClade(Tree *T, nxsstring a, nxsstring b, nxsstring c) {
	NodePtr p;
	NodePtr q;
	std::map<Node *, int, std::less<Node *> > depth;
	PreorderIterator <Node> n ((*T).GetRoot());
	int count = 0;
	q = n.begin();
	while (q)
	{
		depth[q] = count++;
		q = n.next();
	}
	
	//Now get depths for each pair
	vector<int> LCADepthVectorT1;
	(*T).MakeNodeList();
	p = (*T).GetLeafWithLabel(a.c_str());
	q = (*T).GetLeafWithLabel(b.c_str());
	while (depth[p] != depth[q])
	{
		if (depth[p] < depth[q])
			q = q->GetAnc();
		else
			p = p->GetAnc();
	}
	LCADepthVectorT1.push_back(depth[p]);
	p = (*T).GetLeafWithLabel(b.c_str());
	q = (*T).GetLeafWithLabel(c.c_str());
	while (depth[p] != depth[q])
	{
		if (depth[p] < depth[q])
			q = q->GetAnc();
		else
			p = p->GetAnc();
	}
	LCADepthVectorT1.push_back(depth[p]);
	p = (*T).GetLeafWithLabel(a.c_str());
	q = (*T).GetLeafWithLabel(c.c_str());
	while (depth[p] != depth[q])
	{
		if (depth[p] < depth[q])
			q = q->GetAnc();
		else
			p = p->GetAnc();
	}
	LCADepthVectorT1.push_back(depth[p]);
	int T1DepthMax=0; //root has depth 0, others have higher depths.
	int T1DepthMaxIndex=-1;
	//int anum=taxa->FindTaxon(a);
	//int bnum=taxa->FindTaxon(b);
	//int cnum=taxa->FindTaxon(c);
	//cout<<anum<<": "<<a<<" "<<bnum<<": "<<b<<" "<<cnum<<": "<<c<<endl;
	for (int j=0;j<3;j++) {
		if(LCADepthVectorT1[j]>T1DepthMax) {
			T1DepthMax=LCADepthVectorT1[j];
			T1DepthMaxIndex=j;
		}
	}
	if (LCADepthVectorT1[0]==LCADepthVectorT1[1] && LCADepthVectorT1[0]==LCADepthVectorT1[2]) {
		T1DepthMaxIndex=-1;
	}
	//j=0->ab, j=1->bc j=2->ac
	//gsl_matrix_set(TaxonDistance,anum,anum,1+gsl_matrix_get(TaxonDistance,anum,anum));
	//gsl_matrix_set(TaxonDistance,bnum,bnum,1+gsl_matrix_get(TaxonDistance,bnum,bnum));
	//gsl_matrix_set(TaxonDistance,cnum,cnum,1+gsl_matrix_get(TaxonDistance,cnum,cnum));
	vector<nxsstring> returnvector;
	if (T1DepthMaxIndex==0) {
		returnvector.push_back(a);
		returnvector.push_back(b);
	}
	else if (T1DepthMaxIndex==1) {
		returnvector.push_back(b);
		returnvector.push_back(c);
	}
	else if (T1DepthMaxIndex==2) {
		returnvector.push_back(a);
		returnvector.push_back(c);
	}
	else if (T1DepthMaxIndex==-1) {
		returnvector.push_back("nope");
	}
	return returnvector;
}


void BROWNIE::InitializeQuartetCounts() {
	ProgressBar(trees->GetNumTrees());
	for (int treenum=0;treenum<trees->GetNumTrees(); treenum++) {
		Tree T=intrees.GetIthTree(treenum);
		vector<nxsstring> LeafLabelVector;
		NodeIterator <Node> n (T.GetRoot());
		NodePtr currentnode = n.begin();
		while (currentnode)
		{
			if(currentnode->IsLeaf()) {
				LeafLabelVector.push_back(currentnode->GetLabel());
			}
			currentnode = n.next();
		}
		gsl_combination *comb;
		comb=gsl_combination_calloc(T.GetNumLeaves(),4);
		do {
			nxsstring a=LeafLabelVector[gsl_combination_get(comb,0)];
			nxsstring b=LeafLabelVector[gsl_combination_get(comb,1)];
			nxsstring c=LeafLabelVector[gsl_combination_get(comb,2)];
			nxsstring d=LeafLabelVector[gsl_combination_get(comb,3)];
			map< nxsstring, int> LabelToNum;
			map< int, nxsstring > NumToLabel;
			vector<int> LabelNums;
			vector<int> OutputVector;
			LabelToNum[a]=taxa->FindTaxon(a);
			NumToLabel[LabelToNum[a]]=a;
			LabelToNum[b]=taxa->FindTaxon(b);
			NumToLabel[LabelToNum[b]]=b;
			LabelToNum[c]=taxa->FindTaxon(c);
			NumToLabel[LabelToNum[c]]=c;
			LabelToNum[d]=taxa->FindTaxon(d);
			NumToLabel[LabelToNum[d]]=d;
			int minlabel=-1;
			int maxlabel=100000000;
			map<nxsstring, int>::const_iterator itr;
			
			for(itr = LabelToNum.begin(); itr != LabelToNum.end(); ++itr){
				minlabel=GSL_MIN((*itr).second,minlabel);
				maxlabel=GSL_MAX((*itr).second,maxlabel);
				LabelNums.push_back((*itr).second);
			}
			
			vector<nxsstring> taxainclade;
			bool validquartet=false;
			taxainclade=ReturnClade(&T,LeafLabelVector[gsl_combination_get(comb,0)],LeafLabelVector[gsl_combination_get(comb,1)],LeafLabelVector[gsl_combination_get(comb,2)]);
			if(taxainclade.size()==2) {
				taxainclade=ReturnClade(&T,taxainclade[0],taxainclade[1],LeafLabelVector[gsl_combination_get(comb,3)]);	
			}
			if (taxainclade.size()==2) {
				OutputVector.push_back(minlabel);
				if (LabelToNum[taxainclade[0]]==minlabel) {
					
					OutputVector.push_back(LabelToNum[taxainclade[1]]);
					vector<int>remainingoutput;
					for (int i=0;i<4;i++) {
						if(LabelNums[i]!=OutputVector[0] && LabelNums[i]!=OutputVector[1]) {
							remainingoutput.push_back(LabelNums[i]);
						}
					}
					OutputVector.push_back(GSL_MIN(remainingoutput[0],remainingoutput[1]));
					OutputVector.push_back(GSL_MAX(remainingoutput[0],remainingoutput[1]));
				}
				else if (LabelToNum[taxainclade[1]]==minlabel) {
					OutputVector.push_back(LabelToNum[taxainclade[0]]);
					vector<int> remainingoutput;
					for (int i=0;i<4;i++) {
						if(LabelNums[i]!=OutputVector[0] && LabelNums[i]!=OutputVector[1]) {
							remainingoutput.push_back(LabelNums[i]);
						}
					}
					OutputVector.push_back(GSL_MIN(remainingoutput[0],remainingoutput[1]));
					OutputVector.push_back(GSL_MAX(remainingoutput[0],remainingoutput[1]));
				}
				else {
					vector<int> remainingoutput;
					for (int i=0;i<4;i++) {
						if(LabelNums[i]!=OutputVector[0] && LabelNums[i]!=OutputVector[1] && LabelNums[i]!=minlabel) {
							OutputVector.push_back(LabelNums[i]);
						}
					}
					OutputVector.push_back(GSL_MIN(LabelToNum[taxainclade[0]],LabelToNum[taxainclade[1]]));
					OutputVector.push_back(GSL_MAX(LabelToNum[taxainclade[0]],LabelToNum[taxainclade[1]]));
				}
				
			}
			else {
				OutputVector.push_back(maxlabel);
				vector<int> remainingoutput;
				int newmaxlabel=-1;
				for (int i=0;i<4;i++) {
					if(LabelNums[i]!=maxlabel) {
						newmaxlabel=GSL_MAX(LabelNums[i],newmaxlabel);
					}
				}
				OutputVector.push_back(newmaxlabel);
				int newmaxlabel2=-1;
				int newminlabel2=100000000;
				for (int i=0;i<4;i++) {
					if(LabelNums[i]!=maxlabel && LabelNums[i]!=newmaxlabel) {
						newmaxlabel2=GSL_MAX(LabelNums[i],newmaxlabel2);
						newminlabel2=GSL_MAX(LabelNums[i],newminlabel2);
					}
				}
				OutputVector.push_back(newmaxlabel2);
				OutputVector.push_back(newminlabel2);
			}
			quartetcounts[OutputVector]++;
		}
		while (gsl_combination_next(comb) == GSL_SUCCESS);
		gsl_combination_free(comb);
		ProgressBar(0);
	}
	
}

void BROWNIE::HandleNast( NexusToken& token ) {
	bool finishexecuting=true;
	for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (finishexecuting) {
                DoNast();
            }
            break;
        }
        else if( token.Abbreviation("?") ) {
			finishexecuting=false;
		}
        else if(token.Abbreviation("NPercent")) {
            nxsstring numbernexus;
			double initialnumbernexus;
            numbernexus = GetNumber(token);
            initialnumbernexus=atof( numbernexus.c_str() );
			if (initialnumbernexus>1 || initialnumbernexus<0) {
				message="Error: npercent should be between 0 and 1, inclusive";
				PrintMessage();
				finishexecuting=false;
			}
			else {
				npercent=initialnumbernexus;
			}
        }
	}
	
}

void BROWNIE::DoNast() {
	if (quartetcounts.size()==0) {
		message="Now initializing quartet counts";
		PrintMessage();
		InitializeQuartetCounts();
	}
	
}



vector<int> BROWNIE::GetTripletOverlap(ContainingTree *t1, ContainingTree *t2,int taxaincommon) {
    int maxnumber=0;
    int numberdisagree=0;
    int numunresolved=0; //number unresolved in one or the other trees
	int numunresolvedinboth=0; //number unresolved in both trees
	int numunresolvedinT1only=0;
	int numunresolvedinT2only=0;
	(*t1).MakeNodeList(); //added here, took out from containingtree:GetLCADepthVectorCommand
	(*t2).MakeNodeList();
    gsl_combination *c;
    c=gsl_combination_calloc(taxaincommon,3); //this is for triplets, change 3 to 4 (and make some other mods below) for quartets
    vector<nxsstring> LeafLabelVect=(*t1).GetLeafLabelVector();
    do
    {
        maxnumber++;
        vector<int> LCADepthVectorT1=(*t1).GetLCADepthVector(LeafLabelVect[gsl_combination_get(c,0)],LeafLabelVect[gsl_combination_get(c,1)],LeafLabelVect[gsl_combination_get(c,2)]);
        vector<int> LCADepthVectorT2=(*t2).GetLCADepthVector(LeafLabelVect[gsl_combination_get(c,0)],LeafLabelVect[gsl_combination_get(c,1)],LeafLabelVect[gsl_combination_get(c,2)]);
        int T1DepthMax=0; //root has depth 0, others have higher depths.
        int T2DepthMax=0;
        int T1DepthMaxIndex=0;
        int T2DepthMaxIndex=0;
        if ((LCADepthVectorT1[0]==LCADepthVectorT1[1] && LCADepthVectorT1[0]==LCADepthVectorT1[2] && LCADepthVectorT1[1]==LCADepthVectorT1[2]) || (LCADepthVectorT2[0]==LCADepthVectorT2[1] && LCADepthVectorT2[0]==LCADepthVectorT2[2] && LCADepthVectorT2[1]==LCADepthVectorT2[2])) {
            numunresolved++;
			if ((LCADepthVectorT1[0]==LCADepthVectorT1[1] && LCADepthVectorT1[0]==LCADepthVectorT1[2] && LCADepthVectorT1[1]==LCADepthVectorT1[2]) && (LCADepthVectorT2[0]==LCADepthVectorT2[1] && LCADepthVectorT2[0]==LCADepthVectorT2[2] && LCADepthVectorT2[1]==LCADepthVectorT2[2])) {
				numunresolvedinboth++;
			}
			else if ((LCADepthVectorT1[0]==LCADepthVectorT1[1] && LCADepthVectorT1[0]==LCADepthVectorT1[2] && LCADepthVectorT1[1]==LCADepthVectorT1[2])) {
				numunresolvedinT1only++;
			}
			else if ((LCADepthVectorT2[0]==LCADepthVectorT2[1] && LCADepthVectorT2[0]==LCADepthVectorT2[2] && LCADepthVectorT2[1]==LCADepthVectorT2[2])) {
				numunresolvedinT2only++;
			}
        }
        else {
			if(debugmode) {
				cout<<"LCADepthVectorT1 = ("<<LCADepthVectorT1[0]<<" "<<LCADepthVectorT1[1]<<" "<<LCADepthVectorT1[2]<<")"<<endl;
				cout<<"LCADepthVectorT2 = ("<<LCADepthVectorT2[0]<<" "<<LCADepthVectorT2[1]<<" "<<LCADepthVectorT2[2]<<")"<<endl;
			}			
            for (int j=0;j<3;j++) {
                if(LCADepthVectorT1[j]>T1DepthMax) {
                    T1DepthMax=LCADepthVectorT1[j];
                    T1DepthMaxIndex=j;
                }
                if(LCADepthVectorT2[j]>T2DepthMax) {
                    T2DepthMax=LCADepthVectorT2[j];
                    T2DepthMaxIndex=j;
                }
            }
			if(debugmode) {
				cout<<"T1DepthMaxIndex = "<<T1DepthMaxIndex<<endl;
				cout<<"T2DepthMaxIndex = "<<T2DepthMaxIndex<<endl;
			}						
            if (T2DepthMaxIndex!=T1DepthMaxIndex) {
                numberdisagree++;
            }
        }
    }
    while (gsl_combination_next(c) == GSL_SUCCESS);
    gsl_combination_free(c);
    vector<int> tripletoverlapoutput;
    tripletoverlapoutput.push_back(maxnumber);
    tripletoverlapoutput.push_back(numberdisagree);
    tripletoverlapoutput.push_back(numunresolved);
	tripletoverlapoutput.push_back(numunresolvedinboth);
	tripletoverlapoutput.push_back(numunresolvedinT1only);
	tripletoverlapoutput.push_back(numunresolvedinT2only);
	if(debugmode) {
		cout<<"tripletoverlapoutput = (";
		for(int vectorposition=0;vectorposition<6;vectorposition++) {
			cout<<" "<<tripletoverlapoutput[vectorposition];
		}
		cout<<" )"<<endl;
	}						
    return tripletoverlapoutput;
}


void BROWNIE::HandleAccuracy( NexusToken& token )
{
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
			ComputeAccuracy();
            break;
        }
        else if( token.Abbreviation("?") ) {
            message="Usage: Accuracy\n\n";
            message+="Outputs accuracy stats, assuming the first tree is the true tree.\n";
            PrintMessage();
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading ACCURACY command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}


void BROWNIE::ComputeAccuracy() {
	Tree TrueTreeGeneTreeFmt=intrees.GetIthTree(0); //we assume the first tree in the file is the accurate one
	ContainingTree TrueTree;
	TrueTree.SetRoot(TrueTreeGeneTreeFmt.CopyOfSubtree(TrueTreeGeneTreeFmt.GetRoot()));
	//message="TreeNum\tTreeName\tNumTripletsProperlyResolved\tNumTripletsProperlyUnresolved\tNumTripletsImproperlyResolved\tNumTripletsImproperlyUnresolved";
	message="\t\tCorrect\t\tIncorrect\tNumber\nNum\tName\tRes\tUn\tRes\tUn\tSpecies";
	PrintMessage();
	for (int selectedtree=1;selectedtree<trees->GetNumTrees();selectedtree++) {
		ContainingTree ModifiedTrueTree=TrueTree;
		Tree CurrentGeneTreeTreeFmt=intrees.GetIthTree(selectedtree);
		ContainingTree CurrentGeneTree;
		CurrentGeneTree.SetRoot(CurrentGeneTreeTreeFmt.CopyOfSubtree(CurrentGeneTreeTreeFmt.GetRoot()));
		CurrentGeneTree.FindAndSetRoot();
		CurrentGeneTree.Update();
		
		//Get number of species. The idea here is that each node immediately ancestral to one species is a species. This covers the case of a species tree with two clades, representing 2 species, as well as when there's a paraphyletic species tree (some tips connecting to the root and the others to some descendant node)
		typedef map<NodePtr, int> AncestorCountType;
        AncestorCountType AncestorCount;
		NodeIterator <Node> n (CurrentGeneTree.GetRoot());
		NodePtr currentnode = n.begin();
		while (currentnode)
		{
			if (currentnode->IsLeaf()) {
				++AncestorCount[currentnode->GetAnc()];
			}
			currentnode=n.next();
			
		}
		int speciescount=AncestorCount.size();
		if (debugmode) {
			CurrentGeneTree.ReportTreeHealth();
		}
		
		int ntaxincommon=PrepareTreesForTriplet(&ModifiedTrueTree,&CurrentGeneTree);
		vector<int> tripletoverlapoutput=GetTripletOverlap(&ModifiedTrueTree,&CurrentGeneTree,ntaxincommon);
		message="";
		message+=selectedtree+1;
		message+="\t";
		//if (selectedtree==0) {
		//	message+="TRUTH";
		//}
		//else {
			message+=trees->GetTreeName(selectedtree);
		//}
		message+="\t";
		
		double maxnumber=tripletoverlapoutput[0];
		double numberdisagree=tripletoverlapoutput[1];
		double numunresolved=tripletoverlapoutput[2];
		double numunresolvedinboth=tripletoverlapoutput[3];
		double numunresolvedinTTrueonly=tripletoverlapoutput[4];
		double numunresolvedinTOtheronly=tripletoverlapoutput[5];
		//message+=double((tripletoverlapoutput[0]-tripletoverlapoutput[1]-tripletoverlapoutput[4]-tripletoverlapoutput[5])/tripletoverlapoutput[0]);
		message+=(maxnumber-numunresolved-numberdisagree)/maxnumber;
		message+="\t";
		//message+=double(tripletoverlapoutput[3]/tripletoverlapoutput[0]);
		message+=numunresolvedinboth/maxnumber;
		message+="\t";
		//message+=double(tripletoverlapoutput[4]/tripletoverlapoutput[0]);
		message+=(numunresolvedinTTrueonly+numberdisagree)/maxnumber;
		message+="\t";
		//message+=double(tripletoverlapoutput[5]/tripletoverlapoutput[0]);
		message+=numunresolvedinTOtheronly/maxnumber;
		message+="\t";
		message+=GSL_MAX(1,speciescount); //If there's one species, there will only be the root node
		PrintMessage();
	}
}

//Creates batch file for partitioned edge support for a given tree;
void BROWNIE::HandlePartitionedEdgeSupport ( NexusToken& token) {
	int numberofpartitions=1;
	bool donehelp=false;
	for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") && !donehelp) {
			BatchPartitionedEdgeSupport(numberofpartitions);
            break;
        }
        else if( token.Abbreviation("?") ) {
            message="Usage: PartitionedEdgeSupport partitions=<int>\n\n";
            message+="This creates a file (partitionededgesupport.tre) containing the chosen tree and then all trees one NNI move away from this tree, grouped by edges.\n";
            PrintMessage();
			donehelp=true;
        }
		else if( token.Abbreviation("Partitions") ) {
            nxsstring numbernexus = GetNumber(token);
            numberofpartitions=atoi( numbernexus.c_str() ); //convert to int
            message="You have chosen to use ";
            message+=numberofpartitions;
            PrintMessage();
            if (numberofpartitions<1) {
                errormsg = "Error: must select a number greater than zero";
                numberofpartitions=1;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading PARTITIONEDEDGESUPPORT command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
	
}

void BROWNIE::BatchPartitionedEdgeSupport (int numberofpartitions) {
	int totaledges=0;
	ContainingTree OriginalTree;
	Tree OriginalTreeTreeFmt=intrees.GetIthTree(chosentree-1);
	OriginalTree.SetRoot(OriginalTreeTreeFmt.CopyOfSubtree(OriginalTreeTreeFmt.GetRoot()));
	ofstream partedgef;
	nxsstring partedgefile="partitionededgesupport.tre";
	partedgef.open(partedgefile.c_str());
	partedgef<<"#nexus\nbegin trees;\n"; 
	OriginalTree.Update();
	OriginalTree.GetNodeDepths();
	int NodesTouched=0;
	NodeIterator <Node> npe (OriginalTree.GetRoot());
	//OriginalTree.ReportTreeHealth();
    NodePtr currentnodepe = npe.begin();
	//cout<<"starting currentnode is "<<currentnodepe<<endl;
    while (currentnodepe)
    {
		NodesTouched++;
        if (currentnodepe->IsLeaf() || currentnodepe==OriginalTree.GetRoot() ) { // we don't want to do NNI on edges that don't matter
        }
		else {
			nxsstring newlabel="edge";
			newlabel+=NodesTouched;
			currentnodepe->SetLabel(newlabel);
			ContainingTree NewTree=OriginalTree;
			NewTree.NonRandomNNIAtNode(NodesTouched,1);
			partedgef<<"tree edge"<<NodesTouched<<"_res1 = [&R] ";
			NewTree.Write(partedgef);
			partedgef<<endl;
			ContainingTree NewTree2=OriginalTree;
			NewTree2.NonRandomNNIAtNode(NodesTouched,2);
			partedgef<<"tree edge"<<NodesTouched<<"_res2 = [&R] ";
			NewTree2.Write(partedgef);
			partedgef<<endl;
			//OriginalTree.ReportTreeHealth();
			totaledges++;
		}
//cout<<"currentnode before next is "<<currentnodepe<<endl;
        currentnodepe=npe.next();
		//cout<<"currentnode after next is "<<currentnodepe<<endl;

    }
	partedgef<<"[tree original = [&R] ";
	OriginalTree.Write(partedgef);
	partedgef<<"]"<<endl;
	partedgef<<"tree originalwithlabels = [&R] ";
	OriginalTree.SetInternalLabels(true);
	OriginalTree.Write(partedgef);
	partedgef<<endl;
	partedgef<<"end;";
	partedgef.close();
	message="There were two alternate resolutions for each of ";
	message+=totaledges;
	message+=" internal edges, plus the original tree, for ";
	message+=1+(2*totaledges);
	message+=" trees total in partitionededgesupport.tre.";
	PrintMessage();
}


double BROWNIE::ComputeTripletCost(int numberagree,int maxnumber,int ntaxincommon, double Tree1Wt,int Tree1Ntax,double Tree2Wt,int Tree2Ntax, int numberofgenes) {
    double initialscorealgorithm;
    double cdfundermodel;
    if (ntaxincommon<=50) { //change this to use exact values
                            //cdfundermodel=gsl_cdf_binomial_P(numberagree, 2.0/3.0, maxnumber);
        cdfundermodel=CDFvector[ntaxincommon-3][maxnumber-numberagree]; //assumes gene trees are resolved
                                                                        //	cout<<"ntax="<<ntaxincommon<<" cdfundermodelVector="<<cdfundermodel<<" ";
                                                                        //	if(ntaxincommon>5) {
                                                                        //		double x=(1.0*(maxnumber-numberagree))-(2.0*maxnumber/3.0);
                                                                        //		double calculatedvar=sqrt(((1.0/30.0)*double(gsl_sf_choose(ntaxincommon,5)))+((10.0/27.0)*double(gsl_sf_choose(ntaxincommon,4)))+((2.0/9.0)*double(gsl_sf_choose(ntaxincommon,3))));
                                                                        //		assert(calculatedvar>0);
                                                                        //		cout<<"cdfundermodelApprox="<<gsl_cdf_gaussian_P(x,calculatedvar)<<" ";
                                                                        //	}
    }
    else { //use approximation from D.E. Critchlow, D.K. Pearl, and C. Qian 1996. "The triples distance for rooted bifurcating phylogenetic trees" Syst Biol 45(3): 323-334.
        double x=(1.0*(maxnumber-numberagree))-(2.0*maxnumber/3.0);
        double calculatedvar=sqrt(((1.0/30.0)*double(gsl_sf_choose(ntaxincommon,5)))+((10.0/27.0)*double(gsl_sf_choose(ntaxincommon,4)))+((2.0/9.0)*double(gsl_sf_choose(ntaxincommon,3))));
        assert(calculatedvar>0);
        cdfundermodel=gsl_cdf_gaussian_P(x,calculatedvar);
        //	cout<<"ntax="<<ntaxincommon<<" cdfundermodelApprox="<<cdfundermodel<<" ";
        //	for (int i=0;i<=35;i++) {
        //		 x=(1.0*i)-(2.0*35/3.0);
        //		 calculatedvar=sqrt(((1.0/30.0)*double(gsl_sf_choose(7,5)))+((10.0/27.0)*double(gsl_sf_choose(7,4)))+((2.0/9.0)*double(gsl_sf_choose(7,3))));
        //		cout<<"i="<<i<<" cdf="<<gsl_cdf_gaussian_P(x,calculatedvar)<<endl;
        //	}
    }
    if (cdfundermodel>pthreshold) {
        initialscorealgorithm=0; //so we don't penalize random structure
    }
    else {
        initialscorealgorithm=Tree1Wt*Tree2Wt*((1.0/cdfundermodel)-1.0)/(1.0*numberofgenes); //Motivation: treeweighting obvious; use number of identical triplets as score; divide by number of genes since we do O(numberofgenes*numberofgenes) comparisons for structure but only O(numberofgenes) costs for GTP and want the costs to be similar regardless of gene number
		if (debugmode) {
			if (isinf(initialscorealgorithm)==1) {
				cout<<"Error: initialscorealgorithm = "<<initialscorealgorithm<<" cdfundermodel = "<<cdfundermodel;
				if (ntaxincommon<=32) {
					cout<<" calculated using cdfundermodelVector"<<endl;
				}
				else {
					cout<<" calculated using Critchlow et al approximation"<<endl;
				}
			}
		}
    }
                                                                                         //initialscorealgorithm=Tree1Wt*Tree2Wt*(maxnumber-numberdisagree)/numberofgenes;

    return initialscorealgorithm;
}

//Gets trees ready to compare using triplet or quartet scores. After this, just use
//QTValues Q;
//CompareTriplets(t1,t2, Q);
//SummaryStats(Q);
//and get the results you want from Q.
int BROWNIE::PrepareTreesForTriplet(ContainingTree *t1, ContainingTree *t2) {
    int taxaincommon=0;
    (*t1).Update();
    (*t2).Update();
    (*t1).SetLeafNumbers();
    (*t2).SetLeafNumbers();
    //cout<<"starting to compare "<<endl;
    //(*t1).Write(cout);
    //cout<<endl;
    //(*t2).Write(cout);
    //cout<<endl;
    taxaincommon=PruneToOverlappingLeaves(t1, t2);
    if (taxaincommon>2) {
        //   (*t1).BuildLeafClusters();
        //   (*t1).BuildLabelClusters();
        //   (*t1).MakeNodeList();
        //   (*t2).BuildLeafClusters();
        //   (*t2).BuildLabelClusters();
        //   (*t2).MakeNodeList();
        (*t1).SetLeafNumbers();
        (*t2).SetLeafNumbers();
		assert((*t1).GetNumLeaves()==(*t2).GetNumLeaves());

    }
    //cout<<"continuing to compare "<<endl;
    //(*t1).Write(cout);
    //cout<<endl;
    //(*t2).Write(cout);
    //cout<<endl<<taxaincommon<<endl<<endl;
    return taxaincommon;
}

//Uses some code modified from Rod Page's supertrees program to prune leaves from two trees so they have the same leaf set
int BROWNIE::PruneToOverlappingLeaves(ContainingTree *t1, ContainingTree *t2)
{
	//cout<<"\n\n\nStarting Prune\n\n";
	//(*t1).ReportTreeHealth();
	//(*t2).ReportTreeHealth();
    (*t1).MakeNodeList();
    (*t2).MakeNodeList();
    int maxoverlap=(*t1).GetNumLeaves();
    int leavesincommon=0;
    vector<NodePtr> NodesToDelete;
    NodeIterator <Node> n ((*t1).GetRoot());
    NodePtr currentnode = n.begin();
    while (currentnode)
    {
        if (currentnode->IsLeaf()) {
            NodePtr matchingLeaf = (*t2).GetLeafWithLabel(currentnode->GetLabel());
            //cout<<"Current node is "<<currentnode<<" "<<currentnode->GetLabel()<<" and matching leaf is "<<matchingLeaf<<endl;
            if (matchingLeaf != NULL)
            {
                // This leaf is in t1 AND t2. Set the LeafNumber of this leaf
                // in t1 to match that in t2
                leavesincommon++;
                currentnode->SetLeafNumber(leavesincommon); //We need to do this & the following step so with MakeNodeList, the leaf number is 1:#Leaves
                matchingLeaf->SetLeafNumber(leavesincommon);

            }
            else
            {
                // This leaf is not in t2 so we will prune it from the tree
                NodesToDelete.push_back(currentnode);
                //cout<<"Going to delete "<<currentnode<<" with label "<<currentnode->GetLabel()<<" from t1"<<endl;
                maxoverlap--;

            }
        }
        currentnode=n.next();

    }
    if (maxoverlap>0) {
        for (int i=0;i<NodesToDelete.size();i++) {
            (*t1).DeleteSubtendingEdge(NodesToDelete[i]);
            delete (NodesToDelete[i]);
        }
        NodesToDelete.clear();
        (*t1).SuppressInternalNodesWithOneDescendant();
        (*t1).Update();
        (*t1).GetNodeDepths();
		(*t1).SetLeafNumbers();
        (*t2).SetLeafNumbers();
        //now do the same for the reverse comparison

        (*t1).MakeNodeList();
        (*t2).MakeNodeList();
		//(*t1).ReportTreeHealth();
		//(*t2).ReportTreeHealth();
        leavesincommon=0;
        NodeIterator <Node> m ((*t2).GetRoot());
        currentnode = m.begin();
        while (currentnode)
        {
            if (currentnode->IsLeaf()) {
                NodePtr matchingLeaf = (*t1).GetLeafWithLabel (currentnode->GetLabel());
                if (matchingLeaf != NULL)
                {
                    leavesincommon++;
                    currentnode->SetLeafNumber(leavesincommon); //We need to do this & the following step so with MakeNodeList, the leaf number is 1:#Leaves
                    matchingLeaf->SetLeafNumber(leavesincommon);
                }
                else
                {
                    NodesToDelete.push_back(currentnode);
                    //cout<<"Going to delete "<<currentnode<<" with label "<<currentnode->GetLabel()<<" from t2"<<endl;
                }
            }
            currentnode=m.next();
        }
        for (int i=0;i<NodesToDelete.size();i++) {
			//cout<<"\n------------------- \nGoing to delete "<<NodesToDelete[i]<<endl<<endl;
            (*t2).DeleteSubtendingEdge(NodesToDelete[i]);
            delete (NodesToDelete[i]);
			
			//(*t2).ReportTreeHealth();
			//cout<<"\nJust deleted "<<NodesToDelete[i]<<endl<<endl;
        }
        if (leavesincommon>2) {
            (*t2).SuppressInternalNodesWithOneDescendant();
            (*t2).Update();
            (*t2).GetNodeDepths();
        }
    }
			//(*t1).ReportTreeHealth();
			//(*t2).ReportTreeHealth();
    return (GSL_MIN(leavesincommon,maxoverlap));
}



//Eliminate duplicate trees in the FormattedBestTrees vector
void BROWNIE::DelDupes()
{
    vector<ContainingTree> FormattedBestTreesBkup;
    FormattedBestTreesBkup.swap(FormattedBestTrees); //empties besttrees into best treesbkup
    assert(FormattedBestTreesBkup.size()>0);
    FormattedBestTrees.push_back(FormattedBestTreesBkup.back());
    (FormattedBestTrees.back()).FindAndSetRoot();
    (FormattedBestTrees.back()).Update();
    FormattedBestTreesBkup.pop_back(); // transfer a tree back into best trees bkup
    for (int i=0; i<FormattedBestTreesBkup.size(); i++) {
        (FormattedBestTreesBkup[i]).FindAndSetRoot();
        (FormattedBestTreesBkup[i]).Update();
        bool newtree=true;
        for (int j=0; j<FormattedBestTrees.size(); j++) {
            ContainingTree t1=FormattedBestTreesBkup[i];
            ContainingTree t2=FormattedBestTrees[j];
            // cout<<"Compare \n\t";
            // (FormattedBestTreesBkup[i]).Write(cout);
            //cout<<"\n\t";
            //  (FormattedBestTrees[j]).Write(cout);
            int ntaxt1=t1.GetNumLeaves();
            int ntaxt2=t2.GetNumLeaves();
            int ntaxincommon=PrepareTreesForTriplet(&t1,&t2);
            // cout<<"\n\tT1 leaves before = "<<ntaxt1<<" after "<<t1.GetNumLeaves()<<"\n\tT2 leaves before = "<<ntaxt2<<" after "<<t2.GetNumLeaves()<<"\n\tntaxincommon = "<<ntaxincommon<<endl;
            if (ntaxincommon==ntaxt1 && ntaxincommon==ntaxt2) {
                if (ntaxincommon>=3) { //so, all leaves are in common, but is the topology the same? If there are only two or one species, yes, so newtree=false. If there are at least 3 taxa, there are different possible topologies, so maybe, thus check below
                    vector<int> tripletoverlapoutput=GetTripletOverlap(&t1,&t2,ntaxincommon);
                    int numberdisagree=tripletoverlapoutput[1];
                    if (numberdisagree==0) {
                        newtree=false;
                    }
                    //      cout<<"\tNumDisagree = "<<numberdisagree<<endl;
                }
                else {
                    newtree=false; //all leaves in common, must have same topology, so same tree
                }
            }
            //ShowHeader(cout);
            //ShowQTRecord(cout,Q);
            //cout<<endl<<endl;
        }
        if (newtree) {
            FormattedBestTrees.push_back(FormattedBestTreesBkup[i]); //BestTreesBkup is going out of scope soon, anyway, so don't need to bother pop_backing
        }
    }
    for (int k=0; k<FormattedBestTrees.size(); k++) {
        (FormattedBestTrees[k]).Update();
        (FormattedBestTrees[k]).GetNodeDepths();
        (FormattedBestTrees[k]).ClearInternalLabels();
    }
}

/**
* @method HandleCitation [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the CITATION command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleCitation( NexusToken& token )
{
    nxsstring numbernexus;
    bool donenothing=true;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (donenothing) {
				PrintCitations();
            }
            break;
        }
        else if( token.Abbreviation("?") ) {
            donenothing=false;
            message="Usage: Citation\n\n";
            message+="Prints out list of papers that you should read and cite for the analyses\nrun this session.";
            PrintMessage();
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading CHOOSE command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}

void BROWNIE::PrintCitations()
{
	message="Papers whose methods you have used so far in this session:\n[you should certainly read and probably cite them]";
	if (citationarray[0]) {
		message+="\n\nCitation for this program and for rate comparison methods:\n   O'Meara, B.C., C. Ane, M.J. Sanderson, and P.C. Wainwright. 2006. \"Testing for different rates of continuous trait evolution using likelihood.\" Evolution 60(5): 922-933.";
	}
	if (citationarray[1]) {
		message+="\n\nCitation for constant mean, constant pull OU and ACDC transformations (d and g parameters, respectively):\n   Blomberg, S.P., T. Garland, Jr., and A.R. Ives. 2003. \"Testing for phylogenetic signal in comparative data: Behavioral traits are more labile.\" Evolution 57(4) 717-745.";
	}
	if (citationarray[2]) {
		message+="\n\nCitations for Ornstein-Uhlenbeck model with multiple means but one attraction and rate parameter:";
		message+="\n\n   Butler, M.A., King, A.A. 2004. \"Phylogenetic comparative analysis: a modeling approach for adaptive evolution.\" American Naturalist. 164(6):683-695.";
		message+="\n\n   Hansen, T.F., 1997. \"Stabilizing selection and the comparative analysis of adaptation.\" Evolution, 51:1341-1351.";
		message+="\n\n   O'Meara, B.C. Brownie v2.1.1. Distributed by the author at http://www.brianomeara.info/brownie";
				}
	if (citationarray[3]) {
		message+="\n\nCitation for species delimitation approach:";
		message+="\n\n   O'Meara, B.C. 2009. \"New Heuristic Methods for Joint Species Delimitation and Species Tree Inference\" Systematic Biology. doi: 10.1093/sysbio/syp077";
	}
	if (citationarray[4]) {
		message+="\n\nCitation for loss only model:";
		message+="\n\n   McBride, C.S., J.R. Arguello, B.C. O'Meara 2007 \"Five Drosophila Genomes Reveal Nonneutral Evolution and the Signature of Host Specialization in the Chemoreceptor Superfamily\" Genetics: 177(3): 1395-1395.";
	}
	PrintMessage();
}

/**
* @method HandleChoose [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the CHOOSE command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleChoose( NexusToken& token )
{
    nxsstring numbernexus;
    bool donenothing=true;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (donenothing) {
                message="Usage: Choose [tree=<integer>] [char=<integer>]\n\n";
                PrintMessage();
            }
            break;
        }
        else if( token.Abbreviation("Tree") ) {
            donenothing=false;
            numbernexus = GetNumber(token);
            chosentree=atoi( numbernexus.c_str() ); //convert to int
            message="You have chosen tree number ";
            message+=chosentree;
            PrintMessage();
            if (chosentree<1) {
                errormsg = "Error: must select a number greater than zero";
                chosentree=1;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            if (chosentree>trees->GetNumTrees()) {
                errormsg = "Error: you chose tree number ";
                errormsg += chosentree;
                errormsg += " but there are only ";
                errormsg += trees->GetNumTrees();
                errormsg += " trees loaded.\n";
                errormsg += "Tree 1 has been selected by default.";
                chosentree=1;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        else if( token.Abbreviation("Char") ) {
            donenothing=false;
            numbernexus = GetNumber(token);
            chosenchar=atoi( numbernexus.c_str() ); //convert to int
            message="You have chosen char number ";
            message+=chosenchar;
            PrintMessage();
            if (chosenchar<1) {
                errormsg = "Error: must select a number greater than zero";
                chosenchar=1;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            if (chosenchar>continuouscharacters->GetNChar()) {
                errormsg = "Error: you chose char number ";
                errormsg += chosenchar;
                errormsg += " but there are only ";
                errormsg += continuouscharacters->GetNChar();
                errormsg += " characters loaded.\n";
                errormsg += "Character 1 has been selected by default.";
                chosenchar=1;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
		else if( token.Abbreviation("Discretechar") ) {
            donenothing=false;
            numbernexus = GetNumber(token);
            discretechosenchar=-1+atoi( numbernexus.c_str() ); //convert to int
            message="You have chosen discrete character number ";
            message+=discretechosenchar+1;
            PrintMessage();
            if (discretechosenchar<0) {
                errormsg = "Error: must select a number greater than zero";
                discretechosenchar=0;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            if (discretechosenchar>-1+discretecharacters->GetNChar()) {
                errormsg = "Error: you chose char number ";
                errormsg += discretechosenchar;
                errormsg += " but there are only ";
                errormsg += discretecharacters->GetNChar();
                errormsg += " discrete characters loaded.\n";
                errormsg += "Character 1 has been selected by default.";
                discretechosenchar=0;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
		
        else if( token.Abbreviation("?") ) {
            donenothing=false;
            message="Usage: Choose [tree=<integer>] [char=<integer>] [discrete=<integer>]\n\n";
            message+="Selects one tree and/or character to use for subsequent analyses. Continuous characters\n";
			message+="are chosen by default ('char=' or 'c='), discrete characters must be chosen using\n";
			message+="'discretechar=' or 'd='.\n\n";
            message+="The remaining trees and characters are still stored in memory and can be chosen later.\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --\n";
            message+="Tree         <integer-value>                      ";
            if (trees->GetNumTrees()>0) {
                message+=chosentree;
            }
            else {
                message+="[no trees loaded]";
            }
            message+="\nChar         <integer-value>                      ";
            if (continuouscharacters->GetNChar()>0) {
                message+=chosenchar;
            }
            else {
                message+="[no continuous characters loaded]";
            }
			message+="\nDiscrete     <integer-value>                      ";
            if (!discretecharloaded) {
                message+="[no discrete characters loaded]";
                
            }
            else {
                message+=discretechosenchar+1;
            }
			
            message+="\n\n";
            PrintMessage();
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading CHOOSE command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}

void BROWNIE::HandleCompareRandomTrees( NexusToken& token)
{
    bool donenothing=false;
    int ltax=3;
    int rtax=3;
    int nreps=1000;
    bool loopmode=false;
    bool automode=false;
    int requiredminimum=5;
    nxsstring numbernexus;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (donenothing) {
                message="Usage: Compare [ltax=<integer> rtax=<integer> nreps=<integer> [loop] [auto]]\n\nIf you use loop, ltax=taxmax and rtax=taxmin";
                PrintMessage();
            }
            else {
                if (loopmode && !automode) {
                    message="loopmode=yes automode=no\nltax\trtax\tnreps\tdiffs\tobs\tprobability\tcum prob";
                    PrintMessage();
                    int taxmax=GSL_MAX(rtax,ltax);
					int taxmin=GSL_MIN(rtax,ltax);
					nxsstring filename="comparisoncodeMax";
					filename+=taxmax;
					filename+="Min";
					filename+=taxmin;
					filename+=".cpp";
					ofstream comparisonfile;
					comparisonfile.open( filename.c_str() );
					if (debugmode) {
						cout<<"taxmax="<<taxmax<<endl<<"taxmin="<<taxmin<<endl;
					}	
                    for (int ntax=taxmin;ntax<=taxmax;ntax++) {
                        ltax=ntax;
                        rtax=ntax;
						if (debugmode) {
							cout<<"ntax="<<ntax<<endl;
						}	
                        int minnumleaves=GSL_MIN(ltax,rtax);
                        int maxnumsametriplets=(minnumleaves*(minnumleaves-1)*(minnumleaves-2))/6; //minnumleaves choose 3
						vector<int> differencevector(maxnumsametriplets+1,0); // the +1 is so we can have zero in the vector, too
						for (int i=0;i<nreps;i++) {
							ContainingTree t1;
							ContainingTree t2;
							t1.RandomTree(ltax);
							t2.RandomTree(rtax);
							t1.ConvertTaxonNamesToRandomTaxonNumbers();
							t2.ConvertTaxonNamesToRandomTaxonNumbers();
							int taxaincommon=PrepareTreesForTriplet(&t1,&t2);
							// QTValues Q;
							// CompareTriplets(t1,t2, Q);
							// SummaryStats(Q);
							// int numberdisagree=Q.d;
							vector<int> tripletoverlapoutput=GetTripletOverlap(&t1,&t2,taxaincommon);
							int numberdisagree=tripletoverlapoutput[1];
							if (debugmode) {
								t1.Draw(cout);
								t2.Draw(cout);
								cout<<"rep "<<i+1<<" "<<numberdisagree<<endl;
							}
							differencevector[numberdisagree]++;
						}
						
						double cumprob=0;
						comparisonfile<<"contentsofrow.clear();\n";
						for (int i=0; i<differencevector.size(); i++) {
							message="";
							message+=ltax;
							message+="\t";
							message+=rtax;
							message+="\t";
							message+=nreps;
							message+="\t";
							message+=i;
							message+="\t";
							message+=differencevector[i];
							message+="\t";
							char outputstring[11];
							double ratio=(1.0*differencevector[i])/nreps;
							cumprob+=ratio;
							sprintf(outputstring,"%11.10f",ratio);
							message+=outputstring;
							message+="\t";
							sprintf(outputstring,"%11.10f",cumprob);
							message+=outputstring;
							PrintMessage();
							comparisonfile<<"contentsofrow.push_back("<<cumprob<<");\n";
						}
						comparisonfile<<"CDFvector.push_back(contentsofrow); //comparison where ntax="<<ntax<<endl;
                    }
					comparisonfile.close();
                }
                else if (automode) {
                    int minnumleaves=GSL_MIN(ltax,rtax);
                    int maxnumsametriplets=(minnumleaves*(minnumleaves-1)*(minnumleaves-2))/6; //minnumleaves choose 3
					vector<int> differencevector(maxnumsametriplets+1,0); // the +1 is so we can have zero in the vector, too
					int vectorminimum=0;
					int loopcounttotal=0;
					int loopcount=0;
					message="Loops\tNumber of samples in the 4 most similar & 4 least similar scores\tTime remaining";
					PrintMessage();
					nxsstring filename="comparisonof";
					filename+=GSL_MIN(ltax,rtax);
					filename+=".txt";
					ofstream comparisonfile;
					comparisonfile.open( filename.c_str() );
					time_t start,end;
					time (&start);
					int oldminimum=0;
					double secondstotal=-1;
					while (vectorminimum<requiredminimum) {
						ContainingTree t1;
						ContainingTree t2;
						t1.RandomTree(ltax);
						t2.RandomTree(rtax);
						t1.ConvertTaxonNamesToRandomTaxonNumbers();
						t2.ConvertTaxonNamesToRandomTaxonNumbers();
						int taxaincommon=PrepareTreesForTriplet(&t1,&t2);
							// QTValues Q;
							// CompareTriplets(t1,t2, Q);
							// SummaryStats(Q);
							// int numberdisagree=Q.d;
						vector<int> tripletoverlapoutput=GetTripletOverlap(&t1,&t2,taxaincommon);
						int numberdisagree=tripletoverlapoutput[1];
						comparisonfile<<1.0*numberdisagree/maxnumsametriplets<<endl;
						differencevector[numberdisagree]++;
						vectorminimum=requiredminimum*10;
						for (int i=0; i<differencevector.size(); i++) {
							if (differencevector[i]<vectorminimum) {
								vectorminimum=differencevector[i];
							}
						}
						if(vectorminimum>oldminimum) {
							oldminimum=vectorminimum;
							time(&end);
							double timedifference=difftime(end,start);
							secondstotal=(1.0*timedifference)*(1.0*requiredminimum)/(1.0*vectorminimum); // all the 1.0 to make them treated as doubles.
																										 // cout<<"timedif is "<<timedifference<<" sectotal is "<<secondstotal<<" vectorminimum is "<<vectorminimum<<" requiredminimum is "<<requiredminimum<<endl;
						}
						loopcount++;
						loopcounttotal++;
						if (loopcount==5000) {
							message="";
							message+=loopcounttotal;
							message+="\t";
							for (int k=0;k<GSL_MIN(4,differencevector.size());k++) {
								message+="\t";
								message+=differencevector[k];
							}
							message+="\t...";
							for (int k=GSL_MAX(0,differencevector.size()-4);k<differencevector.size();k++) {
								message+="\t";
								message+=differencevector[k];
							}
							message+="\t";
							time(&end);
							if (secondstotal>3600) {
								int hours=int(floor((secondstotal-difftime(end,start))/3600.0));
								int minutes=int(floor(((secondstotal-difftime(end,start))/60.0)-(hours*60)));
								message+=hours;
								message+="h ";
								message+=minutes;
								message+=" min";
							}
							else if (secondstotal>120) {
								int minutes=int(ceil((secondstotal-difftime(end,start))/60.0));
								if (minutes<1) {
									message+="almost done";
								}
								else {
									message+=minutes;
									message+=" min";
								}
							}
							else if (secondstotal>0) {
								int seconds=int(ceil((secondstotal-difftime(end,start))));
								if (seconds<1) {
									message+="almost done";
								}
								else {
									message+=seconds;
									message+=" sec";
								}
							}
							else {
								message+="Too long";
							}
							PrintMessage();
							loopcount=0;
						}
					}
					comparisonfile.close();
					message="";
					message+=loopcounttotal;
					message+="\t";
					for (int k=0;k<GSL_MIN(4,differencevector.size());k++) {
						message+="\t";
						message+=differencevector[k];
					}
					message+="\t...";
					for (int k=GSL_MAX(0,differencevector.size()-4);k<differencevector.size();k++) {
						message+="\t";
						message+=differencevector[k];
					}
					message+="\n\t\t\tElapsed time: ";
					time(&end);
					message+=floor(difftime(end,start)/60);
					message+=" min";
					message+="\n";
					PrintMessage();
					message="ltax\trtax\tnreps\tdiffs\tobs\tproportion\tprobability\tcum prob";
					PrintMessage();
					double cumprob=0;
					int maxcount=0;
					for (int i=0; i<differencevector.size(); i++) {
						if (differencevector[i]>maxcount) {
							maxcount=differencevector[i];
						}
					}
					for (int i=0; i<differencevector.size(); i++) {
						message="";
						message+=ltax;
						message+="\t";
						message+=rtax;
						message+="\t";
						message+=loopcounttotal;
						message+="\t";
						message+=i;
						message+="\t";
						message+=differencevector[i];
						message+="\t";
						char outputstring[11];
						double proportion=(1.0*i/maxnumsametriplets);
						sprintf(outputstring,"%11.10f",maxnumsametriplets);
						message+=outputstring;
						message+="\t";
						double ratio=(1.0*differencevector[i])/loopcounttotal;
						cumprob+=ratio;
						sprintf(outputstring,"%11.10f",ratio);
						message+=outputstring;
						message+="\t";
						sprintf(outputstring,"%11.10f",cumprob);
						message+=outputstring;
						message+="\t";
						for (int k=0;k<floor(10.0*(differencevector[i])/maxcount);k++) {
							message+="*";
						}
						PrintMessage();
					}
					
					gsl_vector * rawdata = gsl_vector_calloc(loopcounttotal);
					int insertionposition=0;
					for (int i=0; i<differencevector.size(); i++) {
						for (int j=0;j<differencevector[i];j++) {
							gsl_vector_set(rawdata,insertionposition,1.0*i/maxnumsametriplets);
							insertionposition++;
						}
					}
					
					message="\n";
					gsl_sort_vector(rawdata);
					int cdfresolution=100;
					for (int k=1;k<=cdfresolution;k++) {
						double proportiondesired=1.0*k/cdfresolution;
						message+="CDF percentile = ";
						message+=proportiondesired;
						message+=" proportion mismatch = ";
						message+=gsl_stats_quantile_from_sorted_data(rawdata->data,1,loopcounttotal,proportiondesired);
						message+="\n";
					}
					
					message+="\nmean is ";
					message+=gsl_stats_mean(rawdata->data, 1, loopcounttotal);
					message+="\nvariance is ";
					message+=gsl_stats_variance(rawdata->data, 1, loopcounttotal);
					PrintMessage();
                }
				
                else {
                    int minnumleaves=GSL_MIN(ltax,rtax);
                    int maxnumsametriplets=(minnumleaves*(minnumleaves-1)*(minnumleaves-2))/6; //minnumleaves choose 3
					vector<int> differencevector(maxnumsametriplets+1,0); // the +1 is so we can have zero in the vector, too
					ProgressBar(nreps);
					for (int i=0;i<nreps;i++) {
						ProgressBar(0);
						ContainingTree t1;
						ContainingTree t2;
						t1.RandomTree(ltax);
						t2.RandomTree(rtax);
						t1.ConvertTaxonNamesToRandomTaxonNumbers();
						t2.ConvertTaxonNamesToRandomTaxonNumbers();
						int taxaincommon=PrepareTreesForTriplet(&t1,&t2);
							// QTValues Q;
							// CompareTriplets(t1,t2, Q);
							// SummaryStats(Q);
							// int numberdisagree=Q.d;
						vector<int> tripletoverlapoutput=GetTripletOverlap(&t1,&t2,taxaincommon);
						int numberdisagree=tripletoverlapoutput[1];
						differencevector[numberdisagree]++;
					}
					message="ltax\trtax\tnreps\tdiffs\tobs\tprobability\tcum prob";
					PrintMessage();
					double cumprob=0;
					for (int i=0; i<differencevector.size(); i++) {
						message="";
						message+=ltax;
						message+="\t";
						message+=rtax;
						message+="\t";
						message+=nreps;
						message+="\t";
						message+=i;
						message+="\t";
						message+=differencevector[i];
						message+="\t";
						char outputstring[11];
						double ratio=(1.0*differencevector[i])/nreps;
						cumprob+=ratio;
						sprintf(outputstring,"%11.10f",ratio);
						message+=outputstring;
						message+="\t";
						sprintf(outputstring,"%11.10f",cumprob);
						message+=outputstring;
						PrintMessage();
					}
					nxsstring filename="comparisonof";
					filename+=GSL_MIN(ltax,rtax);
					filename+=".txt";
					ofstream comparisonfile;
					comparisonfile.open( filename.c_str() );
					gsl_vector * rawdata = gsl_vector_calloc(nreps);
					int insertionposition=0;
					for (int i=0; i<differencevector.size(); i++) {
						for (int j=0;j<differencevector[i];j++) {
							comparisonfile<<i<<endl;
							gsl_vector_set(rawdata,insertionposition,1.0*i/maxnumsametriplets);
							insertionposition++;
						}
					}
					comparisonfile.close();
					gsl_sort_vector(rawdata);
					for (int k=1;k<=10;k++) {
						double proportiondesired=k/10.0;
						cout<<"CDF percentile = "<<proportiondesired<<" proportion mismatch = "<<gsl_stats_quantile_from_sorted_data(rawdata->data,1,nreps,proportiondesired)<<endl;
					}
					cout<<"mean is "<<gsl_stats_mean(rawdata->data, 1, nreps)<<endl;
					cout<<"variance is "<<gsl_stats_variance(rawdata->data, 1, nreps)<<endl;
                }
            }
            break;
        }
        else if( token.Abbreviation("Ltax") ) {
            donenothing=false;
            numbernexus = GetNumber(token);
            ltax=atoi( numbernexus.c_str() ); //convert to int
            if (ltax<3) {
                errormsg = "Error: must select a number greater than 2";
                ltax=3;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        else if( token.Abbreviation("Rtax") ) {
            donenothing=false;
            numbernexus = GetNumber(token);
            rtax=atoi( numbernexus.c_str() ); //convert to int
            if (rtax<3) {
                errormsg = "Error: must select a number greater than 2";
                rtax=3;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        else if( token.Abbreviation("LOop") ) {
            loopmode=true;
        }
        else if ( token.Abbreviation("AUto") ) {
            automode=true;
        }
        else if( token.Abbreviation("Nreps") ) {
            donenothing=false;
            numbernexus = GetNumber(token);
            nreps=atoi( numbernexus.c_str() ); //convert to int
            if (nreps<1) {
                errormsg = "Error: must select a number greater than 0";
                nreps=1000;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
		else if( token.Abbreviation("?") ) {
            donenothing=false;
            message="Usage: Compare [ltax=<integer> rtax=<integer> nreps=<integer> <loop> <auto>]\n\n";
            PrintMessage();
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading COMPARE command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}



/**
* @method HandleSet [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the SET command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleSet( NexusToken& token )
{
    nxsstring numbernexus;
    bool donenothing=true;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (donenothing) {
                message="Usage: Set [maxspecies=<integer>]\n\n";
                PrintMessage();
            }
            break;
        }
        else if( token.Abbreviation("Maxspecies") ) {
            donenothing=false;
            numbernexus = GetNumber(token);
            maxnumspecies=atoi( numbernexus.c_str() ); //convert to int
            message="You have chosen to limit the search to no more than ";
            message+=maxnumspecies;
            PrintMessage();
            if (maxnumspecies<1) {
                errormsg = "Error: must select a number greater than zero";
                maxnumspecies=taxa->GetNumTaxonLabels();
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            if (maxnumspecies>taxa->GetNumTaxonLabels()) {
                errormsg = "Error: you chose a maximum of ";
                errormsg += maxnumspecies;
                errormsg += " but there are only ";
                errormsg += taxa->GetNumTaxonLabels();
                errormsg += " samples loaded.\n";
                maxnumspecies=taxa->GetNumTaxonLabels();
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        else if( token.Abbreviation("Compare")) {
            ContainingTree t1;
            ContainingTree t2;
            t1.RandomTree(6);
            t2.RandomTree(4);
            t1.FindAndSetRoot();
            t1.Update();
            t1.GetNodeDepths();
            t1.ConvertTaxonNamesToRandomTaxonNumbers();
            t2.ConvertTaxonNamesToRandomTaxonNumbers();
            t1.SetLeafNumbers();
            t2.SetLeafNumbers();
            int taxaincommon=PruneToOverlappingLeaves(&t1,&t2);
            cout<<"Taxa in common are "<<taxaincommon<<endl;
            t1.BuildLeafClusters();
            t1.BuildLabelClusters();
            t1.MakeNodeList();
            t2.FindAndSetRoot();
            t2.Update();
            t2.GetNodeDepths();
            t2.BuildLeafClusters();
            t2.BuildLabelClusters();
            t2.MakeNodeList();

            cout<<"clusters tree 1\n";
            t1.ShowClusters();
            cout<<"clusters tree 2\n";
            t2.ShowClusters();
            cout<<"Now comparing\n";
            t1.Draw(cout);
            cout<<endl;
            t2.Draw(cout);
            cout<<endl;
            QTValues Q;
            CompareTriplets(t1,t2, Q);
            SummaryStats(Q);
			int taxaincommon2=PrepareTreesForTriplet(&t1,&t2);
							// QTValues Q;
							// CompareTriplets(t1,t2, Q);
							// SummaryStats(Q);
							// int numberdisagree=Q.d;
			vector<int> tripletoverlapoutput=GetTripletOverlap(&t1,&t2,taxaincommon2);
			int numberdisagree=tripletoverlapoutput[1];
// int numberdisagree=Q.d;
            int maxnumber=tripletoverlapoutput[0];
            cout<<"Out of "<<maxnumber<<" only "<<numberdisagree<<" disagreed\n";
        }
        else if( token.Abbreviation("Randtree")) {
            donenothing=false;
            ContainingTree testtree;
            testtree.RandomTree(6);
            cout<<endl;
            // testtree.Write(cout);
            testtree.ConvertTaxonNamesToOrderedTaxonNumbers();
            cout<<endl;
            //testtree.Write(cout);
            cout<<endl;
            testtree.Update();
            //testtree.Draw(cout);
            cout<<endl<<"Tree has "<<testtree.GetNumLeaves()<<" leaves and "<<testtree.GetNumInternals()<<" internals"<<endl;
            testtree.TestRerooting();
            // testtree.Update();
            //  testtree.Draw(cout);
            cout<<endl<<"Tree has "<<testtree.GetNumLeaves()<<" leaves and "<<testtree.GetNumInternals()<<" internals"<<endl;
        }
        else if( token.Abbreviation("Write")) {
            donenothing=false;
            ContainingTree testtree;
            int nspecies=0;
            for (int vectorpos=0;vectorpos<convertsamplestospecies.size();vectorpos++) {
                if (convertsamplestospecies[vectorpos]>nspecies) {
                    nspecies=convertsamplestospecies[vectorpos];
                }
            }
            testtree.RandomTree(nspecies);
            testtree.Draw(cout);
            testtree.ConvertTaxonNamesToRandomTaxonNumbers();
            testtree.Draw(cout);
            cout<<OutputForGTP(&testtree);
            //printUsage();
            //cout<<"ReturnScore = "<<ReturnScore(OutputForGTP(&testtree),unrooted);
        }
        else if( token.Abbreviation("?") ) {
            donenothing=false;
            message="Usage: Set [maxspecies=<integer>]\n\n";
            message+="Sets the maximum number of species to test.\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --\n";
            message+="MaxSpecies   <integer-value>                      ";
            message+=maxnumspecies;
            PrintMessage();
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading CHOOSE command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}



void BROWNIE::HandleModel( NexusToken& token )
{
	// Retrieve all tokens for this command, stopping only in the event
	// of a semicolon or an unrecognized keyword
	//
    for(;;)
    {
        token.GetNextToken();
		
        if( token.Equals(";") ) {
            break;
        }
        else if( token.Abbreviation("?") ) {
			//message="Usage: Model type= [BM | OU | ACDC] assign=[1 | CS | CC | TI | TC] ou=[M | MV]\n";
        //    message="Usage: Model type = [ BM1 | OU1 | ACDC1 | BMS | BMC | BMAN | BMAO | OUSMVA | OUCMVA | OUSMA | OUCMA | OUSM | OUCM | TSBMI |TSBMC]\n\n"; //add this once everything's been tested
            message="Usage: Model type = [ BM1 | BMS | BMC | OUSM | ... ] states = ( vector )\n";
            message+="";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------------------------------- Current setting -----";
            message+="\n Type        BM1 | BMS | BMC | OUSM | OUCM                                  ";

		//	message+="\n Type        BM | OU | ACDC                                                           ";
		//	message+="\n Assign      1 | CS | CC | TI | TC | RBS | RBO                                                    ";
		//	message+="\n OU          M | MV                                                                   ";
         //   message+="\n Type        BM1 | OU1 | ACDC1 | BM[S|C] | OU[SMVA|CMVA|SMA|CMA|CM|SM] | TS[BMS|BMC]  ";
            if (chosenmodel==1 || chosenmodel==2) {
                message+="BM1";
            }
            else if(chosenmodel==3) {
                message+="OU1";
            }
            else if(chosenmodel==4) {
                message+="ACDC";
            }
            else if(chosenmodel==5) {
                message+="BMS";
            }
            else if(chosenmodel==6) {
                message+="BMC";
            }
            else if(chosenmodel==12) {
                message+="OUSM";
            }            
			else if(chosenmodel==20) {
				message+="Pagel's Kappa";
			}
			else if(chosenmodel==21) {
				message+="Pagel's Delta";
			}
			else if(chosenmodel==22) {
				message+="Pagel's Lambda";
			}
            else {
                message+=chosenmodel;
            }
			message+="\n             OUSM | OUCM";
            message+="\n States      vector                                                        (";
            for (int i=0; i<staterestrictionvector.size();i++) {
                message+=" ";
                message+=staterestrictionvector[i];
            }
            message+=" )";
			message+="\n Changes     int                                                           ";
			message+=maxstartstops;
			/*if((optimalvalueslabels.size())>1) {
				message+="\n ManualOptima";
				cout<<"optimalvaluescontinuouschar size="<<optimalvaluescontinuouschar->size<<endl;
				cout<<"optimaldiscretecharstatefreq size="<<optimaldiscretecharstatefreq->size<<endl;
				for (int i=0; i<optimalvaluescontinuouschar->size;i++) {
					cout<<" "<<gsl_vector_get(optimalvaluescontinuouschar,i);
				}
				cout<<endl;
			}*/
           // message+="\n\nType\nBM1    = Brownian motion, one rate parameter";
			//ROse suggests OUCMA approach would be better
			
           // message+="\nOU: Ornstein-Uhlenbeck";
		// message+="\nACDC: Acceleration/deceleration model of Blomberg et al (2003)";
           // message+="\nAssign=1: Assign one model over the whole tree";
          //  message+="\nAssign=CS: Assign one model for each discrete Character State mapped on the tree";
           // message+="\nAssign=CC: Assign one model on branches where the discrete Character Changes, and a different model to the other branches";
         //   message+="\nAssign=TI: Assign one model for each Time Interval";
         //   message+="\nAssign=TC: Assign one model on branches that Cross a Time interval, and a different model to the other branches";
         //   message+="\nOU=M: Allow means (optimal values) to vary across OU models, but keep variance and attraction parameters the same";
          //  message+="\nOU-MV: Allow different means and variances across OU models, but keep attraction parameter the same";
			message+="\n\nBM1     = Brownian motion, one rate parameter";
			message+="\nBMS     = Brownian motion, with different rate parameters for each state on a tree";
			message+="\nBMC     = Brownian motion, with one rate parameter for branches with state changes and another for branches without changes";
			//message+="\nOU1     = Blomberg et al. Ornstein-Uhlenbeck (one attraction parameter (d), one mean)";
			//message+="\nACDC    = Blomberg et al. Acceleration/Deceleration (g parameter)";
			message+="\nOUSM    = Ornstein-Uhlenbeck with one mean per discrete state (attraction and rate parameters constant across tree)";
			message+="\nOUCM    = Ornstein-Uhlenbeck with independent means on branches with and without changes in a discrete character";
			message+="\n         (attraction and rate parameters constant across tree)";
			//message+="\nKappa   = Pagel's Kappa parameter (continuous char)";
			//message+="\nDelta   = Pagel's Delta parameter (continuous char)"; //NOTE: Calculation of this seems off
			//message+="\nLambda  = Pagel's Lambda parameter (continuous char)"; //NOTE: Calculation of this seems off

		//	message+="\nOUSMVA = Ornstein-Uhlenbeck, with different mean, variance, and attraction parameters for each state on a tree";
		//	message+="\nOUCMVA = Ornstein-Uhlenbeck, with different mean, variance, and attraction parameters on branches with and without state changes";
		//	message+="\nOUSMA  = Ornstein-Uhlenbeck, with different mean and attraction parameters for each state, but one variance parameter";
		//	message+="\nOUCMA  = Ornstein-Uhlenbeck, with different mean and attraction parameters on branches with and without state changes, but one variance parameter";
		//	message+="\nTSMBI  = Brownian motion, with different rate parameters for each time interval (set using timeslice)";
		//	message+="\nTSMBC  = Brownian motion, with different rate parameters for branches that cross and do not cross time interval boundaries";
			
			message+="\n\nState vector allows restrictions, so that character states 0 and 2, for example, may be\n     viewed by the program as identical. To do this, you'd enter:\n\n       states=(0 1 0 2 3 4 5 6 7 8)";
			/*if(optimalvalueslabels.size()>1) {
				message+="\n\nManualOptima allows you to change the values recorded from the last optimization\n  and use this in character simulations.";
				message+="\n  Current parameter labels:";
				for (int i=0;i<optimalvalueslabels.size();i++) {
					message+=" ";
					message+=optimalvalueslabels[i];
				}
			}*/
			message+="\n\nChanges is the maximum number of times a particular character state can be present on a root to tip lineage.\n  For example, if a taxon sister to all other taxa starts in state 0, changes to state 1,\n  and then changes to state 0, state 0 has beeen present on that branch twice."; 
            PrintMessage();
        }
		else if( token.Abbreviation("Changes") ) {
            nxsstring numbernexus = GetNumber(token);
            int newmaxchanges=atoi( numbernexus.c_str() ); //convert to int
            if (newmaxchanges<1) {
                errormsg = "Error: must select a number greater than zero";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
			else {
				maxstartstops=2*newmaxchanges;
			}
        }
		
		else if (token.Abbreviation("Manualoptima")) {
            int numparams=(optimalvaluescontinuouschar->size)-1;
            for (int currentparam=0; currentparam<numparams; currentparam++) {
                double newvalue;
                if ((cout<<"Parameter "<<currentparam+1<<" (Old value: "<<gsl_vector_get(optimalvaluescontinuouschar,currentparam)<<") New value: ") && (!(cin >> newvalue))) {
                    cout << "Using old value:";
                    cin.clear();
					newvalue=gsl_vector_get(optimalvaluescontinuouschar,currentparam);
                    //cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }
                gsl_vector_set(optimalvaluescontinuouschar,currentparam,newvalue);
            }
			GetOptimalVCVAndTraitsContinuous();
			int ntax=optimalVCV->size1;
			gsl_vector *tipsresid=gsl_vector_calloc(ntax);
			tipsresid=GetTipValues(globalchosentaxset,chosenchar);
			for (int taxon=0;taxon<ntax;taxon++) {
				gsl_vector_set(tipsresid,taxon,gsl_vector_get(tipsresid,taxon)-gsl_vector_get(optimalTraitMeans,taxon));
			}
			gsl_vector_set(optimalvaluescontinuouschar,numparams,GetLScore(optimalVCV,tipsresid,1));
			message="New lscore with current model settings: ";
			message+=gsl_vector_get(optimalvaluescontinuouschar,numparams);
			PrintMessage();
        }
        else if( token.Abbreviation("Type") ) {
			citationarray[0]=true;
            nxsstring chosenmodelinput=GetFileName(token);
            if (token.Abbreviation("BM1")) {
                if (tipvariancetype==1) {
                    chosenmodel=2;
                }
                else {
                    chosenmodel=1;
                }
                message="You have chosen a single rate Brownian motion model";
                PrintMessage();
            }
            else if(token.Abbreviation("OU1")) {
                chosenmodel=3;
                citationarray[1]=true;
                message="You have chosen the OU1 model";
                PrintMessage();
            }
            else if(token.Abbreviation("ACDC1")) {
                chosenmodel=4;
                citationarray[1]=true;
                message="You have chosen the ACDC model";
                PrintMessage();
            }
            else if(token.Abbreviation("BMS")) {
                chosenmodel=5;
				citationarray[0]=true;
                message="You have chosen Brownian motion with one rate parameter per state";
                PrintMessage();
            }
            else if(token.Abbreviation("BMC")) {
                chosenmodel=6;
				citationarray[0]=true;
                message="You have chosen Brownian motion with different rates on branches with and without state changes.";
                PrintMessage();
            }
            else if(token.Abbreviation("OUSMVA")) {
                chosenmodel=7;
                message="You have chosen the OU model with different means, brownian rate parameters, and attraction parameters for each state.";
                PrintMessage();
            }
            else if(token.Abbreviation("OUCMVA")) {
                chosenmodel=8;
                message="You have chosen the OU model with different means, brownian rate parameters, and attraction parameters on branches with and without state changes.";
                PrintMessage();
            }
            else if(token.Abbreviation("OUSMA")) {
                chosenmodel=10;
                message="You have chosen the OU model with different means and attraction parameters for each state.";
                PrintMessage();
            }
            else if(token.Abbreviation("OUCMA")) {
                chosenmodel=11;
                message="You have chosen the OU model with different means and attraction parameters on branches with and without state changes.";
                PrintMessage();
            }
            else if(token.Abbreviation("OUSM")) {
                chosenmodel=12;
                message="You have chosen the OU model with different means for each state.";
				citationarray[2]=true;
                PrintMessage();
            }
            else if(token.Abbreviation("OUCM")) {
                chosenmodel=13;
                message="You have chosen the OU model with different means on branches with and without state changes.";
				citationarray[2]=true;
                PrintMessage();
            }
            else if(token.Abbreviation("TSBMI")) {
                chosenmodel=14;
                message="You have chosen to have a different rate parameter in each time interval.";
                PrintMessage();
            }
            else if(token.Abbreviation("TSBMC")) {
                chosenmodel=15;
                message="You have chosen to have a different rate parameter on branches that do and do not cross a time interval.";
                PrintMessage();
            }
            else if (token.Abbreviation("BMAN")) {
                chosenmodel=16;
                message="You have chosen to try all possible assignments of zero rates to branches, leaving the other brlen unchanged";
                PrintMessage();
            }
            else if (token.Abbreviation("BMAO")) {
                chosenmodel=17;
                message="You have chosen to try all possible assignments of zero rates to branches, assigning a branch length of one to the other branches";
                PrintMessage();
            }
            else if (token.Abbreviation("BMPN")) {
                chosenmodel=18;
                message="You have chosen to try, at each internal node, each assignment of length 0 to one descendant branch and the original brlen to the other";
                PrintMessage();
            }
            else if (token.Abbreviation("BMPO")) {
                chosenmodel=19;
                message="You have chosen to try, at each internal node, each assignment of length 0 to one descendant branch and a length of 1 to the other";
                PrintMessage();
            }
			else if (token.Abbreviation("KAppa")) {
                chosenmodel=20;
                message="You have chosen to try Pagel's Kappa parameter model";
                PrintMessage();
            }
			else if (token.Abbreviation("DElta")) {
                chosenmodel=21;
                message="You have chosen to try Pagel's Delta parameter model";
                PrintMessage();
            }
			else if (token.Abbreviation("LAmbda") || token.Abbreviation("LAmda")) {
                chosenmodel=22;
                message="You have chosen to try Pagel's Lambda parameter model";
                PrintMessage();
            }
            else {
                errormsg = "Unexpected option (";
                errormsg += chosenmodelinput;
                errormsg += ") encountered reading Model command";
                throw XNexus( errormsg);
            }
        }
        else if( token.Abbreviation("States") ) {
            token.GetNextToken();
            token.GetNextToken(); //eat the equals sign
            vector<int> temporarystatevector;
            if (!token.Equals("(")) {
                errormsg="Expecting next token to be a \'(\' but instead got ";
                errormsg+=token.GetToken();
                throw XNexus( errormsg);
            }
            int inputcount=0;
            while (!token.Equals(")")) {
                nxsstring numbernexus;
                numbernexus=GetNumberOnly(token);
                if (numbernexus!=")") {
                    temporarystatevector.push_back(atoi( numbernexus.c_str() ));
                    inputcount++;
                }
                else {
                    break;
                }
            }
            if (inputcount<1) {
                errormsg="You should have entered at least one number";
                throw XNexus( errormsg);
            }
            else {
                int maxstate=0;
                for (int i=0;i<GSL_MIN(temporarystatevector.size(),staterestrictionvector.size());i++) {
                    staterestrictionvector[i]=temporarystatevector[i];
                    
                }
                for (int i=0;i<staterestrictionvector.size();i++) {
                    maxstate=GSL_MAX(staterestrictionvector[i],maxstate);
                }
                bool complete=false;
                while (!complete) {
                    for (int j=0;j<=maxstate;j++) {
                        int count=0;
                        for (int i=0;i<staterestrictionvector.size();i++) {
                            if(staterestrictionvector[i]==j) {
                                count++;
                            }
                        }
                        if (count==0) {
                            for (int i=0;i<staterestrictionvector.size();i++) {
                                if(staterestrictionvector[i]>=j) {
                                    staterestrictionvector[i]=staterestrictionvector[i]-1;
                                }
                            }
                            maxstate--;
                            complete=false;
                            break;
                        }
                        else {
                            complete=true;
                        }
                    }
                }
                //int sizediff=staterestrictionvector.size()-temporarystatevector.size();
               //// if (sizediff>0) {
               //     for (int i=temporarystatevector.size();i<staterestrictionvector.size();i++) {
               //         staterestrictionvector[i]=i-1-(-1+temporarystatevector.size()-maxstate);
               //     }
              //  }
            }
			
            ///below here is old way of reading vector
          //  nxsstring numbernexus;
           // numbernexus=GetNumber(token);
          //  staterestrictionvector[0]=atoi( numbernexus.c_str() );
         //   for (int position=1;position<=9;position++) {
         //       numbernexus=GetNumberOnly(token);
         //       if (debugmode) {
         //           message="Now reading for position ";
          //          message+=position;
          //          message+=" state is ";
          //          message+=numbernexus;
          //          PrintMessage();
          //      }
           //     staterestrictionvector[position]=atoi( numbernexus.c_str() );
           // }
        }
		
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading Model command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}

void BROWNIE::HandleLoss ( NexusToken& token )
{
	citationarray[4]=true;
	// Retrieve all tokens for this command, stopping only in the event
	// of a semicolon or an unrecognized keyword
	//
	bool donesomething=false;
    for(;;)
    {
        token.GetNextToken();
		
        if( token.Equals(";") ) {
			if (donesomething==false) {
			//double rateA=0.0000000000000000001;
			//while (rateA<100) {
			//	double neglnL=CalculateDiscreteLindy1(rateA);
			//	message="rateA=";
			//	message+=rateA;
			//	message+=" -lnL=";
			//	message+=neglnL;
			//	PrintMessage();
			//	rateA*=100.0;
			//}
			//double neglnL=CalculateDiscreteLindy1(0.022322);
			//message="rateA=";
			//message+=0.022322;
			//message+=" -lnL=";
			//message+=neglnL;
			//PrintMessage();
				
			//BROWNIE my_fnA(maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput,trees, taxa, assumptions, characters);
				gsl_vector *optimalrateA=gsl_vector_calloc(3);
				//gsl_vector_memcpy(optimalrate,my_fn.OptimizeRateWithGivenTipVariance());
				gsl_vector_memcpy(optimalrateA,LindyGeneralOptimization(1));
			//for(int i=0;i<optimalrateA->size;i++){
			//	message="i=";
			//	message+=i;
			//	message+=" val is ";
			//	message+=gsl_vector_get(optimalrateA,i);
			//	PrintMessage();
			//}
			//BROWNIE my_fnB(maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput,trees, taxa, assumptions, characters);
				gsl_vector *optimalrateB=gsl_vector_calloc(5);
				gsl_vector_memcpy(optimalrateB,LindyGeneralOptimization(2));
		//	for(int i=0;i<optimalrateB->size;i++){
		//		message="i=";
		//		message+=i;
		//		message+=" val is ";
		//		message+=gsl_vector_get(optimalrateB,i);
		//		PrintMessage();
		//	}
				double modelAaicc=(2.0*gsl_vector_get(optimalrateA,2))+2+4.0/((discretecharacters->GetNChar())-2); //AICc, n=1;
				double modelBaicc=(2.0*gsl_vector_get(optimalrateB,4))+4+12.0/((discretecharacters->GetNChar())-3);
				double modelAaic=(2.0*gsl_vector_get(optimalrateA,2))+2*1;
				double modelBaic=(2.0*gsl_vector_get(optimalrateB,4))+2*2;				
				message="ModelA\n\trate = ";
				message+=gsl_vector_get(optimalrateA,0);
				message+=" +/- ";
				message+=gsl_vector_get(optimalrateA,1);
				message+="\n\t-lnL = ";
				message+=gsl_vector_get(optimalrateA,2);
				message+="\n\tAIC = ";
				message+=modelAaic-GSL_MIN(modelAaic,modelBaic);
				message+="\n\tAICc = ";
				message+=modelAaicc-GSL_MIN(modelAaicc,modelBaicc);				
				message+="\nModelB\n\trate0 = ";
				message+=gsl_vector_get(optimalrateB,0);
				message+=" +/- ";
				message+=gsl_vector_get(optimalrateB,2);
				message+="\n\trate1 = ";
				message+=gsl_vector_get(optimalrateB,1);
				message+=" +/- ";
				message+=gsl_vector_get(optimalrateB,3);
				message+="\n\t-lnL = ";
				message+=gsl_vector_get(optimalrateB,4);
				message+="\n\tAIC = ";
				message+=modelBaic-GSL_MIN(modelAaic,modelBaic);
				message+="\n\tAICc = ";
				message+=modelBaicc-GSL_MIN(modelAaicc,modelBaicc);	
				PrintMessage();
			}
            break;
        }
        else if( token.Abbreviation("?") ) {
			donesomething=true;
            message="Usage: Loss \n\n";
            message+="";
            message+="Available options: none\n\n";
		message+="This evaluates two models. For each model, the ancestral state is assumed to be 1 (presence) and the model is applied across all characters. \n  ModelA: Same loss rate on all branches\n  ModelB: Independent loss rates on branches with label 0 and label 1.";
            PrintMessage();
        }
		
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading Loss command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}

void BROWNIE::HandlePagelDiscrete ( NexusToken& token)
{
	//This does, basically Pagel 1994's discrete model for correlation of two traits, but does it on both binary and multistate traits.
	bool donesomething=false;
    bool donenothing=true;
	bool treeloop=false;
 	nxsstring tmessage;
    ofstream tablef;
    nxsstring tablefname;
    bool tablef_open=false;
    bool name_provided=false;	
	bool appending=true;
	bool replacing=false;
	globalstates=false;
	int char1=0; //is 0 offset
	int char2=1;
	int pagelchosenmodel=1;
	discretechosenmodel=4;
	tmessage="";
	int vectorposition=0;
	int position=-1; //used only in user-set model
	for(;;)
    {
        token.GetNextToken();
		
        if( token.Equals(";") ) {
			if (donesomething==false) {
				if( appending && replacing ) {
					errormsg = "Cannot specify APPEND and REPLACE at the same time";
					throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
				}		
				bool exists = FileExists( tablefname.c_str() );
				bool userok = true;
				if (appending && name_provided) {
					tablef_open = true;
					tablef.open( tablefname.c_str(), ios::out | ios::app );
					message = "\nAppending to discrete model output file (creating it if need be) ";
					message += tablefname;
					PrintMessage();
				}
				else if (name_provided) {
					if( exists && !replacing && !UserSaysOk( "Ok to replace?", "Discrete model output file specified already exists" ) )
						userok = false;
					if( userok && !tablef_open) {
						tablef_open = true;
						tablef.open( tablefname.c_str() );
					}
					if( exists && userok ) {
						message = "\nReplacing discrete model output file ";
						message += tablefname;
					}
					else if( userok ) {
						message = "\nDiscrete model output file ";
						message += tablefname;
						message += " opened";
					}
					else {
						errormsg = "Aborting the discrete optimization so as not to overwrite the file.\n";
						throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
					}
					PrintMessage();
				}
				
				assert(char1>=0);
				assert(char2>=0);
				int originaldiscretechosenchar=discretechosenchar;
				discretechosenchar=discretecharacters->GetNChar(); //due to 0 offset, this will be number of the new char 
				int char1numstates=discretecharacters->GetObsNumStates(char1); //0 offset in discretematrix
				int char2numstates=discretecharacters->GetObsNumStates(char2);
				int stateConversionMatrix[char1numstates][char2numstates];
				vector<double> temporaryratematfixedvector; //Will be a vector containing JUST the fixed values
				vector<int> temporaryratematassignvector; //Will be a vector containing ints corresponding to "pointers" to either  fixed or variable values. If entries are non-negative,
													//they point to entries in temporaryratematfixedvector (i.e., value of 2 means the rate is whatever is stored at temporaryratematfixedvector[2])
				string tempfreerateletterstring;		//Allows mapping of letters on input to negative values in temporaryratematassignvector vector. New letters are appended, old ones are looked up
				usermatrix=""; //just a string to store the description
				
				int newstate=0;
				for (int i=0; i<char1numstates; i++) {
					for (int j=0; j<char2numstates; j++) {
						stateConversionMatrix[i][j]=newstate;
						//cout<<"i="<<i<<" j="<<j<<" newstate="<<newstate<<endl;
						newstate++;
					}
				}
				discretecharacters->AddCharacters(1);
                int ntax=taxa->GetNumTaxonLabels();
				//for (int taxonid=0; taxonid<ntax; taxonid++) {
				//	cout<<endl<<"taxonid="<<taxonid;
				//	for (int charid=0; charid<discretecharacters->GetNChar(); charid++) {
				//		cout<<" "<<discretecharacters->GetInternalRepresentation(taxonid,charid);
				//	}
				//}
				//cout<<endl;
				discretecharacters->MaximizeSymbols();
				for (int taxonid=0; taxonid<ntax; taxonid++) {
					int index1=discretecharacters->GetInternalRepresentation(taxonid,char1);
					int index2=discretecharacters->GetInternalRepresentation(taxonid,char2);
					//cout<<"index1="<<index1<<" index2="<<index2<<endl;
					//cout<<"stateConversionMatrix[discretecharacters->GetState(taxonid,char1)][atoi(discretecharacters->GetState(taxonid,char2)]))="<<stateConversionMatrix[index1][index2]<<endl;
					discretecharacters->SetState(taxonid,-1+discretecharacters->GetNChar(),stateConversionMatrix[index1][index2]);
				}
				localnumbercharstates=(discretecharacters->GetObsNumStates(char1))*(discretecharacters->GetObsNumStates(char2));
				
				//Now create the rate matrix
				string assigmentMatrix[localnumbercharstates][localnumbercharstates]; //Stores the entries
				string assigmentMatrixChar1[char1numstates][char1numstates]; //Stores variables for rates
				string assigmentMatrixChar2[char2numstates][char2numstates]; //Stores variables for rates
				int freeparameterLabelIndex=0;
				string freeparameterLabels="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
				
				
				/* this creates a matrix like:
					
					0	1
					0	-	a
					1	b	0
					
					if the model is non-time reversible within a character (i.e.,possibly different gain loss rates, as specified with
																			"INDNon" or "DEPNon") or
					
					0	1
					0	-	a
					1	a	0
					
					otherwise. Same for char2, though its parameters start with c (if two params have been used up with char1) or b (if one param has been used up with char1)
					
					Basic idea is that the optimization function assigns one rate parameter for every unique letter. So two rates with the same letter are assigned to the same parameter
					*/
				for (int i=0; i<char1numstates; i++) {
					for (int j=0; j<char1numstates; j++) {
						if (i<j) { //upper part of matrix
							if (pagelchosenmodel==1 || pagelchosenmodel==3) { //a reversible model within each char (gain and loss rate within a char is equal)
								assigmentMatrixChar1[i][j]=freeparameterLabels[freeparameterLabelIndex];
								assigmentMatrixChar1[j][i]=freeparameterLabels[freeparameterLabelIndex]; //so fill upper and lower parts with same rate
								freeparameterLabelIndex++;
							}
							else if (pagelchosenmodel==2 || pagelchosenmodel==4) { //a nonreversible model within each char
								assigmentMatrixChar1[i][j]=freeparameterLabels[freeparameterLabelIndex];
								freeparameterLabelIndex++; //so, we use the next rate parameter for the other rate
								assigmentMatrixChar1[j][i]=freeparameterLabels[freeparameterLabelIndex]; //so fill upper and lower parts with same rate
								freeparameterLabelIndex++;
							}
						}
						else if (i==j) {
							assigmentMatrixChar1[i][j]="0";
						}
					}
				}
				for (int i=0; i<char2numstates; i++) {
					for (int j=0; j<char2numstates; j++) {
						if (i<j) { //upper part of matrix
							if (pagelchosenmodel==1 || pagelchosenmodel==3) { //a reversible model within each char (gain and loss rate within a char is equal)
								assigmentMatrixChar2[i][j]=freeparameterLabels[freeparameterLabelIndex];
								assigmentMatrixChar2[j][i]=freeparameterLabels[freeparameterLabelIndex]; //so fill upper and lower parts with same rate
								freeparameterLabelIndex++;
							}
							else if (pagelchosenmodel==2 || pagelchosenmodel==4) { //a nonreversible model within each char
								assigmentMatrixChar2[i][j]=freeparameterLabels[freeparameterLabelIndex];
								freeparameterLabelIndex++; //so, we use the next rate parameter for the other rate
								assigmentMatrixChar2[j][i]=freeparameterLabels[freeparameterLabelIndex]; //so fill upper and lower parts with same rate
								freeparameterLabelIndex++;
							}
						}
						else if (i==j) {
							assigmentMatrixChar2[i][j]="0";
						}
					}
				}
				
				//Fill in Q matrix by row, omitting the diagonal elements. So, first 00 by [00], 01, 02, etc., then 01 by 00, [01], 02, etc.
				/*idea here is to create a matrix like
					00	01	10	11
					00	-	c1	a1	0
					01	d1	-	0	a2
					10	b1	0	-	c2
					11	0	b2	d2	-
					
					when the original char matrices were
Char1:				Char2
					0	1				0	1
					0	-	a			0	-	c
					1	b	-			1	d	-
					
					If the model is time reversible ("INDRev" or "DEPRev"), then a==b and c==d, and this is represented in the assigmentMatrixChar* matrices by having the same symbol for gain and loss rates (or multistate analog). 
					Thus, rate of 00->10 (a1) should equal 10->00 (b1): the gain/loss rates of character 1 when character 2 is in state 0 are equal. So, the first time we sample rate a or b, they should be the same variable label in assigmentMatrixChar*.
					If the model is not time-reversible ("INDNon" or "DEPNon"), a!=b and there will be different labels in assigmentMatrixChar*
					If the model is independent ("INDRev" or "INDNon"), the rate of 00->10 (a1) and 01->11 (a2) should be equal: the rate of the first character changing from 0 to 1 is the same no matter what state the second character is. Thus, the rates
					a1 and a2 should have the same label, that stored in assigmentMatrixChar1 (with label "a").
					If the model is dependent ("DEPRev" or "DEPNon"), a1 and a2 may be different, so we'll need a different label for a2 (this can be done by relabeling "a" in assigmentMatrixChar1). However, if the model is time-reversible, a2==b2, so in that case, when updating "a", we'd also have to update "b")

*/
				
				for (int i=0; i<char1numstates; i++) {
					for (int j=0; j<char2numstates; j++) {
						int rowstate=stateConversionMatrix[i][j];
						//The above loops gives rows: ij=00, 01, 02, .. 0N, 10, 11, 12, ...1N, ...MN, where M and N are the the maximum state number for each char (so, if both traits are binary, M=N=1: highest state you can have is 1 (the other state is 0).
						//Now, for each row, move across columns.
						for (int k=0; k<char1numstates; k++) {
							for (int l=0; l<char2numstates; l++) {
								int colstate=stateConversionMatrix[k][l];
								//kl is the label for the column. So, cell with row ij (say, 02) and column kl (say, 12) is the cell in the qmatrix for going from state pair 02 to state pair 12
								if (i==k && j==l) { //we must be on a diagonal. 
								}
								else if (i!=k && j!=l) {
									//we're at a cell representing a change of both characters simultaneously (01->23, for example). Set to zero (can't change both at same instant).
									assigmentMatrix[rowstate][colstate]="0";
								}
								else if (colstate>rowstate) { //upper diagonal part
									//there's a change in just one of the chars
									if (i!=k) { //change in char 1
										assigmentMatrix[rowstate][colstate]=assigmentMatrixChar1[i][k]; //is the rate of going from state i to state k
										if (pagelchosenmodel==1 || pagelchosenmodel==2) { // traits are independent, so a1==a2. If they are also non-reversible, a1==b1, but this is taken care of by the assigmentMatrixChar* matrix
											assigmentMatrix[colstate][rowstate]=assigmentMatrixChar1[k][i];											
										}
										else if (pagelchosenmodel==3 || pagelchosenmodel==4) { //traits are dependent, so a1 and a2 differ. So, we'll have to get a new label for a2
											freeparameterLabelIndex++;
											assigmentMatrixChar1[i][k]=freeparameterLabels[freeparameterLabelIndex];
											if (pagelchosenmodel==3) { //reversible, so backwards rate must have the new rate parameter
												assigmentMatrixChar1[k][i]=assigmentMatrixChar1[i][k];
											}
											else if (pagelchosenmodel==4) { //non-reversible, so backwards rate also needs a new parameter
												freeparameterLabelIndex++;
												assigmentMatrixChar1[k][i]=freeparameterLabels[freeparameterLabelIndex];
											}
											assigmentMatrix[colstate][rowstate]=assigmentMatrixChar1[k][i];
										}
									}
									else { //change in char2
										assigmentMatrix[rowstate][colstate]=assigmentMatrixChar2[j][l]; //is the rate of going from state i to state k
										if (pagelchosenmodel==1 || pagelchosenmodel==2) { // traits are independent, so a1==a2. If they are also non-reversible, a1==b1, but this is taken care of by the assigmentMatrixChar* matrix
											assigmentMatrix[colstate][rowstate]=assigmentMatrixChar2[l][j];											
										}
										else if (pagelchosenmodel==3 || pagelchosenmodel==4) { //traits are dependent, so a1 and a2 differ. So, we'll have to get a new label for a2
											freeparameterLabelIndex++;
											assigmentMatrixChar2[j][l]=freeparameterLabels[freeparameterLabelIndex];
											if (pagelchosenmodel==3) { //reversible, so backwards rate must have the new rate parameter
												assigmentMatrixChar2[l][j]=assigmentMatrixChar2[j][l];
											}
											else if (pagelchosenmodel==4) { //non-reversible, so backwards rate also needs a new parameter
												freeparameterLabelIndex++;
												assigmentMatrixChar2[l][j]=freeparameterLabels[freeparameterLabelIndex];
											}
											assigmentMatrix[colstate][rowstate]=assigmentMatrixChar2[l][j];
										}
										
									}
									
								}
								
							}
						}
						
					}
				}
				
				//Finally, we have our final matrix of rates. Now we just have to convert it to something ready for input to the DiscreteGeneralOptimization function
				/*message="Input matrix:\n";
				for (int i=0; i<char1numstates; i++) {
					for (int j=0; j<char2numstates; j++) {
						message+="\t(";
						message+=i;
						message+=",";
						message+=j;
						message+=")";
					}
				}*/
				for (int i=0; i<char1numstates; i++) {
					for (int j=0; j<char2numstates; j++) {
						int rowstate=stateConversionMatrix[i][j];
						message+="\n(";
						message+=i;
						message+=",";
						message+=j;
						message+=")";
						//The above loops gives rows: ij=00, 01, 02, .. 0N, 10, 11, 12, ...1N, ...MN, where M and N are the the maximum state number for each char (so, if both traits are binary, M=N=1: highest state you can have is 1 (the other state is 0).
						//Now, for each row, move across columns.
						for (int k=0; k<char1numstates; k++) {
							for (int l=0; l<char2numstates; l++) {
								int colstate=stateConversionMatrix[k][l];
								if (colstate!=rowstate) { //so, avoid diagonals
									string assigmentMatrixEntry=(assigmentMatrix[colstate][rowstate]);
									usermatrix+=assigmentMatrixEntry;
									usermatrix+=" ";
									if (isalpha(assigmentMatrixEntry[0])) { //Is a letter -- means that parameter is free to vary, but has same value as other rates with that value
										string::size_type loc = tempfreerateletterstring.find( assigmentMatrixEntry[0], 0 );
										if( loc != string::npos ) {
											temporaryratematassignvector.push_back(-1*(loc+1));
										} else {
											tempfreerateletterstring.append(1,assigmentMatrixEntry[0]);
											temporaryratematassignvector.push_back(-1*(tempfreerateletterstring.size()));         
										}
									}
									else { //is a number, which means it's a fixed value
										temporaryratematassignvector.push_back(temporaryratematfixedvector.size());
										temporaryratematfixedvector.push_back(atof( assigmentMatrixEntry.c_str()));
									}
								}
							}
						}
					}
				}
			
				for (int i=0; i<char1numstates; i++) {
					for (int j=0; j<char2numstates; j++) {
						int rowstate=stateConversionMatrix[i][j];
						cout<<endl<<"("<<i<<","<<j<<"):\t";
						for (int k=0; k<char1numstates; k++) {
							for (int l=0; l<char2numstates; l++) {
								int colstate=stateConversionMatrix[k][l];
								if (colstate!=rowstate) { //so, avoid diagonals
									cout<<assigmentMatrix[rowstate][colstate]<<"\t";
								}
								else {
									cout<<"-\t";
								}
							}
						}
					}
				}
				
				cout<<endl;
				cout<<"temporaryratematassignvector";
				for (int i=0; i<temporaryratematassignvector.size(); i++) {
					cout<<" "<<temporaryratematassignvector[i];
				}
				cout<<endl;
				
				cout<<endl;
				cout<<"temporaryratematfixedvector";
				for (int i=0; i<temporaryratematfixedvector.size(); i++) {
					cout<<" "<<temporaryratematfixedvector[i];
				}
				cout<<endl;
				
				cout<<"tempfreerateletterstring = "<<tempfreerateletterstring<<endl;
				
			freerateletterstring=tempfreerateletterstring;
			ratematfixedvector.swap(temporaryratematfixedvector);
			ratematassignvector.swap(temporaryratematassignvector);
			numberoffreefreqs=localnumbercharstates-1;
			numberoffreerates=freerateletterstring.size();
			numberoffreeparameters=numberoffreerates+numberoffreefreqs; //stored globally for later calculation of AIC/AICc
			
			tmessage="Tree\tTree weight\tTree name\tChar1\tChar2\tModel\tStateFreq\t\tneglnL\tK\tAIC\tAICc\t";
			for (int n=0; n<localnumbercharstates; n++) {
				tmessage+="P(";
				tmessage+=n;
				tmessage+=")\t";
			}
			for (int i=0; i<char1numstates; i++) {
				for (int j=0; j<char2numstates; j++) {
					int rowstate=stateConversionMatrix[i][j];
						//The above loops gives rows: ij=00, 01, 02, .. 0N, 10, 11, 12, ...1N, ...MN, where M and N are the the maximum state number for each char (so, if both traits are binary, M=N=1: highest state you can have is 1 (the other state is 0).
						//Now, for each row, move across columns.
					for (int k=0; k<char1numstates; k++) {
						for (int l=0; l<char2numstates; l++) {
							int colstate=stateConversionMatrix[k][l];
								//kl is the label for the column. So, cell with row ij (say, 02) and column kl (say, 12) is the cell in the qmatrix for going from state pair 02 to state pair 12
							if (colstate!=rowstate) {
								tmessage+="q_(";
								tmessage+=i;
								tmessage+=",";
								tmessage+=j;
								tmessage+=")_->_(";
								tmessage+=k;
								tmessage+=",";
								tmessage+=l;
								tmessage+=")";
								tmessage+="\t";
								
							}
						}
					}
				}
			}
			if (tablef_open && (!exists || !appending) ) {
				tablef<<tmessage;
			}
			
			
			
			
			
			
			gsl_vector* output=DiscreteGeneralOptimization();
			nxsstring treename=trees->GetTreeName(chosentree-1);
			message="Tree = ";
			message+=chosentree;
			message+=": ";
			message+=treename;
			message+="\n";
			assert(output->size>0);
			double likelihood=gsl_vector_get(output,-1+output->size);
			double K=1.0*numberoffreeparameters;
			double aicc=(2.0*likelihood)+2.0*K+2.0*K*(K+1.0)/(1.0*ntax-K-1.0); //AICc, n=1;
			double aic=(2.0*likelihood)+2.0*K;
			char outputstring[14];
			message+="\n  model = ";
			if (pagelchosenmodel==1) {
				message+="IndRev (no correlation between traits, equal forward-reverse rates of trait evolution within a character)";
			}
			else if (pagelchosenmodel==2) {
				message+="IndNon (no correlation between traits, forward-reverse rates of trait evolution within a character may vary)";
			}
			else if (pagelchosenmodel==3) {
				message+="DepRev (correlation between traits, equal forward-reverse rates of trait evolution within a character)";
			}
			else if (pagelchosenmodel==4) {
				message+="DepNon (correlation between traits, forward-reverse rates of trait evolution within a character may vary)";
			}
			message+="\n  state frequency = ";
			if (discretechosenstatefreqmodel==1) {
				message+="Uniform";
			}
			else if (discretechosenstatefreqmodel==2) {
				message+="Empirical";
			}
			else if (discretechosenstatefreqmodel==3) {
				message+="Equilibrium";
			}
			else if (discretechosenstatefreqmodel==4) {
				message+="Optimized";
			}
			else if (discretechosenstatefreqmodel==5) {
				message+="User";
			}
			
			message+="\n  -lnL = ";
			sprintf(outputstring,"%14.6f",likelihood);
			message+=outputstring;
			message+="\n  AIC  = ";
			sprintf(outputstring,"%14.6f",aic);
			message+=outputstring;
			message+="\n  AICc = ";
			sprintf(outputstring,"%14.6f",aicc);
			message+=outputstring;
			
			cout<<endl<<"output = (";
			for (int i=0; i<output->size; i++) {
					cout<<endl<<gsl_vector_get(output,i);
			}
			cout<<endl<<")"<<endl;
			
			 vectorposition=0;
			for (int i=0; i<char1numstates; i++) {
				for (int j=0; j<char2numstates; j++) {
					int rowstate=stateConversionMatrix[i][j];
						//The above loops gives rows: ij=00, 01, 02, .. 0N, 10, 11, 12, ...1N, ...MN, where M and N are the the maximum state number for each char (so, if both traits are binary, M=N=1: highest state you can have is 1 (the other state is 0).
						//Now, for each row, move across columns.
					for (int k=0; k<char1numstates; k++) {
						for (int l=0; l<char2numstates; l++) {
							int colstate=stateConversionMatrix[k][l];
								//kl is the label for the column. So, cell with row ij (say, 02) and column kl (say, 12) is the cell in the qmatrix for going from state pair 02 to state pair 12
							if (colstate!=rowstate) {
								message+="\n  q_(";
								message+=i;
								message+=",";
								message+=j;
								message+=")_->_(";
								message+=k;
								message+=",";
								message+=l;
								message+=")";
								message+=" = ";
								if (ratematassignvector[vectorposition]>=0) { //means there's an assigned rate
									message+=ratematfixedvector[(ratematassignvector[vectorposition])];
									message+=" FIXED";
								}
								else {
									position=-1*(1+ratematassignvector[vectorposition]);
									message+=gsl_vector_get(output,position);
									if (gsl_vector_get(output,position)<0.00000001 && !nonnegvariables) {
										message+="  Warning: an estimate near zero sometimes makes estimating other parameters, and therefore the lnL, very imprecise. Play with numopt or the model";
									}
									
								}
								vectorposition++;
							}
						}
					}
				}
			}
			
			 vectorposition=0;
			message+="\n\nQ Matrix (rounded)\n";
			for (int i=0; i<char1numstates; i++) {
				for (int j=0; j<char2numstates; j++) {
					message+="\t(";
					message+=i;
					message+=",";
					message+=j;
					message+=")";
				}
			}
			for (int i=0; i<char1numstates; i++) {
				for (int j=0; j<char2numstates; j++) {
					int rowstate=stateConversionMatrix[i][j];
					message+="\n(";
					message+=i;
					message+=",";
					message+=j;
					message+=")";
					for (int k=0; k<char1numstates; k++) {
						for (int l=0; l<char2numstates; l++) {
							int colstate=stateConversionMatrix[k][l];
							if (colstate!=rowstate) {
								message+="\t";
								if (ratematassignvector[vectorposition]>=0) { //means there's an assigned rate
									message+=0.001*floor(1000.0*ratematfixedvector[(ratematassignvector[vectorposition])]);
								}
								else {
									position=-1*(1+ratematassignvector[vectorposition]);
									message+=0.001*floor(1000.0*gsl_vector_get(output,position));
								}
								vectorposition++;
							}
							else {
								message+="\t-";
							}
						}
					}
				}
			}
			PrintMessage();
			}
			
            break;
        }
        else if (token.Abbreviation("Treeloop") ) {
            nxsstring yesnotreeloop=GetFileName(token);
            if (yesnotreeloop[0] == 'n') {
                treeloop=false;
            }
            else {
                treeloop=true;
            }
        }		
		else if( token.Abbreviation("CHAR1") ) {
            donenothing=false;
            nxsstring numbernexus = GetNumber(token);
            char1=-1+atoi( numbernexus.c_str() ); //convert to int
            message="You have chosen discrete character number ";
            message+=char1+1;
			message+=" for the first trait";
            PrintMessage();
            if (char1<0) {
                errormsg = "Error: must select a number greater than zero";
                char1=0;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            if (char1>-1+discretecharacters->GetNChar()) {
                errormsg = "Error: you chose char number ";
                errormsg += char1;
                errormsg += " but there are only ";
                errormsg += discretecharacters->GetNChar();
                errormsg += " discrete characters loaded.\n";
                char1=0;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
		else if( token.Abbreviation("CHAR2") ) {
            donenothing=false;
            nxsstring numbernexus = GetNumber(token);
            char2=-1+atoi( numbernexus.c_str() ); //convert to int
            message="You have chosen discrete character number ";
            message+=char2+1;
			message+=" for the second trait";
            PrintMessage();
            if (char2<0) {
                errormsg = "Error: must select a number greater than zero";
                char2=0;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            if (char2>-1+discretecharacters->GetNChar()) {
                errormsg = "Error: you chose char number ";
                errormsg += char2;
                errormsg += " but there are only ";
                errormsg += discretecharacters->GetNChar();
                errormsg += " discrete characters loaded.\n";
                char2=0;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
		
		
		else if (token.Abbreviation("Globalstates") ) {
            nxsstring yesnoglobalstates=GetFileName(token);
            if (yesnoglobalstates[0] == 'n') {
                globalstates=false;
            }
            else {
                globalstates=true;
            }
        }		
		else if( token.Abbreviation("Replace") ) {
            nxsstring yesnoreplace=GetFileName(token);
            if (yesnoreplace[0] == 'n') {
                replacing=false;
            }
            else {
                replacing=true;
				appending=false;
            }
        }
        else if( token.Abbreviation("APpend") ) {
            nxsstring yesnoappend=GetFileName(token);
            if (yesnoappend[0] == 'n') {
                appending=false;
            }
            else {
                appending=true;
            }
        }
		else if( token.Abbreviation("File") ) {
            tablefname = GetFileName(token);
            name_provided = true;
        }		
		else if( token.Abbreviation("Model") ) {
            nxsstring chosenmodelinput=GetFileName(token);
			//int numbercharstates=(discretecharacters->GetObsNumStates(discretechosenchar));
			int numberofrates=(localnumbercharstates*localnumbercharstates)-localnumbercharstates;
			int ntax=taxa->GetNumTaxonLabels();
            if (token.Abbreviation("INDRev")) {
				pagelchosenmodel=1;
				message="The two traits are independent and rates for each estimated under a time-reversible model (gain and loss rates for a given trait set to be equal)";
                PrintMessage();
            }
			else if (token.Abbreviation("INDNon")) {
				pagelchosenmodel=2;
				message="The two traits are independent and rates for each estimated under a non-time-reversible model (gain and loss rates for a given trait allowed to vary)";
                PrintMessage();
            }
			else if (token.Abbreviation("DEPRev")) {
				pagelchosenmodel=3;
				message="The two traits are dependent and rates are estimated under a time-reversible model";
                PrintMessage();
            }
			else if (token.Abbreviation("DEPNon")) {
				pagelchosenmodel=4;
				message="The two traits are dependent and rates are estimated under a non-time-reversible model";
                PrintMessage();
            }
			
			else if(token.Abbreviation("User")) {
				pagelchosenmodel=5;
                message="You have chosen a user-specified model. But this isn't implemented yet.";
				PrintMessage();
            }
			
            else {
                errormsg = "Unexpected option (";
                errormsg += chosenmodelinput;
                errormsg += ") encountered reading Model command";
                throw XNexus( errormsg);
            }
        }		
		else if( token.Abbreviation("Freq") ) {
            nxsstring chosenmodelinput=GetFileName(token);
            if (token.Abbreviation("Uniform")) {
				discretechosenstatefreqmodel=1;
				message="You have chosen equal frequencies for all states";
                PrintMessage();
            }
            else if(token.Abbreviation("EMpirical")) {
				discretechosenstatefreqmodel=2;
                message="You have chosen to use empirical state frequencies";
                PrintMessage();
            }
            else if(token.Abbreviation("EQuilibrium")) {
				discretechosenstatefreqmodel=3;
                message="You have chosen to use equilibrium state frequencies";
                PrintMessage();
            }
            else if(token.Abbreviation("Optimized")) {
				discretechosenstatefreqmodel=4;
                message="You have chosen to optimize state frequencies";
                PrintMessage();
            }
			else if(token.Abbreviation("User")) {
				discretechosenstatefreqmodel=5;
                message="You have chosen to use user-set state frequencies";
                PrintMessage();
            }
            else {
                errormsg = "Unexpected option (";
                errormsg += chosenmodelinput;
                errormsg += ") encountered reading Freq command";
                throw XNexus( errormsg);
            }
        }	
        else if( token.Abbreviation("?") ) {
			donesomething=true;
            message="Usage: PagelCorrelation char1= char2=[model=] [freq=] [ratemat=] [treeloop=] [breaknum=] [file=] [append=] [replace=] [globalstates=]\n\n";
			message+="This is a function to do the correlation tests of Pagel (1994), extended to allow different root frequencies and multistate traits\n";
            message+="Available options:\n\n";
            message+="Keyword ------- Option type ----------------------------- Current setting -----";
			message+="\n Char1          <int>                                     ?";
			message+="\n Char2          <int>                                     ?";
			
			message+="\n Model          <string>                                  ?";
			message+="\n RateMat        (<vector>)                                ";
				message+="Unspecified";			
            message+="\n Treeloop       No|Yes                                    *No";
 			message+="\n File           <file name>                               *None";
			message+="\n Append         No|Yes                                    *Yes";
			message+="\n Replace        No|Yes                                    *No";
			message+="\n GlobalStates   No|Yes                                    *No";
            message+="\n                                                        *Option is nonpersistent\n\n";
			message+="Char1 & 2: Allows you to specify the chosen characters. REQUIRED.\n";
			message+="Model: Allows you to specify the model. The model can be INDependent or DEPendent between the two characters. Within each character,\n";
			message+="       the forward and reverse rates can be the same (REVersible) or allowed to vary (NONreversible). The available options are therefore\n";
			message+="       IndRev, DepRev, IndNon, DepNon. IndNon and DepNon correspond to Pagel (1994) analyses (if state frequency is UNIFORM)\n";
			message+="Freq: The probability of each state at the root can be based on the EMPIRICAL distribution at the tips, can be SET by the user\n";
			message+="      (using the statevector command), can be OPTIMIZEd as part of the model, can be set to EQUILIBRIUM frequencies (the\n";
			message+="      frequencies expected with the optimized rate matrix given infinitely-long branches), or can be set to be UNIFORM (equal).\n";
			message+="Treeloop: Allows the analysis to be run across all the trees. A weighted average is returned (using tree weights such as posterior\n";
			message+="          probabilities or bootstrap frequencies if they are available) as well as values for the individual trees.\n";
			message+="File: Saves all output into a tab-delimited file\n";
			message+="Append: If the output file exists, appends to it rather than overwrites it. Turned on by default.\n";
			message+="Replace: If set to yes, if the output file already exists it will be quietly replaced.\n";
			message+="GlobalStates: If no, the number of character states assumed for each character is the maximum number of observed states for just that one character.\n";
			message+="              If yes, the number of states for each character is the maximum number of states observed for any character, even if the observed character\n";
			message+="              is lacking some of those states. This is useful if, for example, you have simulated a three state character on the tree for parametric\n";
			message+="              bootstrapping but some of the simulations result in characters with just states 0 and 1. Globalstates will automatically be set to yes\n";
			message+="              if allchar=y";
            PrintMessage();
        }
		else if( token.Abbreviation("Ratemat") ) {
			errormsg = "Custom user matrices are not yet available. But please feel free to submit a patch to fix this!";
			throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
		}
		
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading Discrete command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}

void BROWNIE::HandleDiscrete( NexusToken& token )
{
	// Retrieve all tokens for this command, stopping only in the event
	// of a semicolon or an unrecognized keyword
	//
	bool donesomething=false;
	bool treeloop=false;
	bool charloop=false;
	int breaknum=0;
	bool reconstruct=false;
	nxsstring tmessage;
    ofstream tablef;
    nxsstring tablefname;
    bool tablef_open=false;
    bool name_provided=false;	
	bool appending=true;
	bool replacing=false;
    for(;;)
    {
        token.GetNextToken();
		
        if( token.Equals(";") ) {
			if (donesomething==false) {
				if(allchar) { //If we are doing inferences from all char at once, use max number of states as nstates for each char
					globalstates=true;
				}
				if( appending && replacing ) {
					errormsg = "Cannot specify APPEND and REPLACE at the same time";
					throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
				}		
				bool exists = FileExists( tablefname.c_str() );
				bool userok = true;
				if (appending && name_provided) {
					tablef_open = true;
					tablef.open( tablefname.c_str(), ios::out | ios::app );
					message = "\nAppending to discrete model output file (creating it if need be) ";
					message += tablefname;
					PrintMessage();
				}
				else if (name_provided) {
					if( exists && !replacing && !UserSaysOk( "Ok to replace?", "Discrete model output file specified already exists" ) )
						userok = false;
					if( userok && !tablef_open) {
						tablef_open = true;
						tablef.open( tablefname.c_str() );
					}
					if( exists && userok ) {
						message = "\nReplacing discrete model output file ";
						message += tablefname;
					}
					else if( userok ) {
						message = "\nDiscrete model output file ";
						message += tablefname;
						message += " opened";
					}
					else {
						errormsg = "Aborting the discrete optimization so as not to overwrite the file.\n";
						throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
					}
					PrintMessage();
				}
				
				int origchosentree=chosentree;
				int origdiscretechosenchar=discretechosenchar;
				int charlooplimit=1;
				if(charloop==true) {
					charlooplimit=discretecharacters->GetNChar();
				}
				int looplimit=1;
				int ntax=taxa->GetNumTaxonLabels();
				if (treeloop==true) {
					looplimit=trees->GetNumTrees();
				}
				if (globalstates || charloop) {
					localnumbercharstates=numbercharstates;
				}
				tmessage="Tree\tTree weight\tTree name\tChar\tNo. states this char\tNo. states used\tModel\tStateFreq\tVariable Sites Only\tneglnL\tK\tAIC\tAICc\t";
				for (int n=0; n<localnumbercharstates; n++) {
					tmessage+="P(";
					tmessage+=n;
					tmessage+=")\t";
				}
				for (int i=0;i<localnumbercharstates;i++) {
					for (int j=0;j<localnumbercharstates;j++) {
						if (i!=j) {
							tmessage+="q_";
							tmessage+=i;
							tmessage+="_";
							tmessage+=j;
							tmessage+="\t";
						}
					}
				}
				if (tablef_open && (!exists || !appending) ) {
					tablef<<tmessage;
				}
				if (!globalstates && tablef_open) {
					message="WARNING: if some characters have fewer than ";
					message+=numbercharstates;
					message+=" character states, the output in the output table will be wrong for those chars";
					tmessage="\n";
					tmessage+=message;
					PrintMessage();
					tmessage+="\n";
					tablef<<tmessage;
				}
				for (int charnum=1;charnum<=charlooplimit;charnum++) {
					if (charloop==true) {
						discretechosenchar=charnum;
					}
					double weighttotal=0;
					localnumbercharstates=(discretecharacters->GetObsNumStates(discretechosenchar));
					if (globalstates) {
						localnumbercharstates=numbercharstates;
					}

					
					for (int treenum=1;treenum<=looplimit;treenum++) {
						if (treeloop==true) {
							chosentree=treenum;
						}
						double treeweight=trees->GetTreeWeight(chosentree-1);
						weighttotal+=treeweight;
						nxsstring treename=trees->GetTreeName(chosentree-1);
						gsl_vector* output=DiscreteGeneralOptimization();
						
						/*for (int i=0; i<output->size;i++) {
							cout<<gsl_vector_get(output,i)<<"\t";
						}
						cout<<endl;*/
						message="Tree = ";
						message+=chosentree;
						message+=": ";
						message+=treename;
						message+="\n";
						if (tablef_open) {
							tmessage="\n";
							tmessage+=chosentree;
							tmessage+="\t";
							tmessage+=treeweight;
							tmessage+="\t";
							tmessage+=treename;
							tmessage+="\t";
							if (allchar) {
								tmessage+="ALL";
							}
							else {
								tmessage+=discretechosenchar+1;
							}
							tmessage+="\t";
							tmessage+=discretecharacters->GetObsNumStates(discretechosenchar);
							tmessage+="\t";
							tmessage+=numbercharstates;
							tmessage+="\t";
							if (discretechosenmodel==1) {
								tmessage+="Equal";
							}
							else if (discretechosenmodel==2) {
								tmessage+="Rev";
							}
							else if (discretechosenmodel==3) {
								tmessage+="NonRev";
							}
							else if (discretechosenmodel==4) {
								tmessage+="User";
								tmessage+=": ( ";
								tmessage+=usermatrix;
								tmessage+=")";
							}
							tmessage+="\t";
							
							if (discretechosenstatefreqmodel==1) {
								tmessage+="Uniform";
							}
							else if (discretechosenstatefreqmodel==2) {
								tmessage+="Empirical";
							}
							else if (discretechosenstatefreqmodel==3) {
								tmessage+="Equilibrium";
							}
							else if (discretechosenstatefreqmodel==4) {
								tmessage+="Optimized";
							}
							else if (discretechosenstatefreqmodel==5) {
								tmessage+="User";
							}
							tmessage+="\t";
							if(variablecharonly) {
								tmessage+="Yes\t";
							}
							else {
								tmessage+="No\t";
							}
							tablef<<tmessage;
						}
						assert(output->size>0);
						double likelihood=gsl_vector_get(output,-1+output->size);
						double K=1.0*numberoffreeparameters;
						double aicc=(2.0*likelihood)+2.0*K+2.0*K*(K+1.0)/(1.0*ntax-K-1.0); //AICc, n=1;
						double aic=(2.0*likelihood)+2.0*K;
						if (tablef_open) {
							tmessage="";
							tmessage+=likelihood;
							tmessage+="\t";
							tmessage+=numberoffreeparameters;
							tmessage+="\t";
							tmessage+=aic;
							tmessage+="\t";
							tmessage+=aicc;
							tmessage+="\t";
							tablef<<tmessage;
						}
						char outputstring[14];
						message+="\n  -lnL = ";
						sprintf(outputstring,"%14.6f",likelihood);
						message+=outputstring;
						message+="\n  AIC  = ";
						sprintf(outputstring,"%14.6f",aic);
						message+=outputstring;
						message+="\n  AICc = ";
						sprintf(outputstring,"%14.6f",aicc);
						message+=outputstring;
						tmessage="";
						int vectorposition=0;
						int position=-1; //used only in user-set model
						message+="\n  Rates: ";
						if (discretechosenmodel==1) {
							message+="equal";
						}
						else if (discretechosenmodel==2) {
							message+="reversible";
						}
						else if (discretechosenmodel==3) {
							message+="nonreversible";
						}
						else if (discretechosenmodel==4) {
							message+="user, with entry ( ";
							message+=usermatrix;
							message+=")";
						}
						
						for (int i=0; i<localnumbercharstates;i++) {
							for (int j=0; j<localnumbercharstates;j++) {
								if (i!=j) {							
								//cout<<i<<" "<<j<<" "<<vectorposition<<" "<<ratematassignvector.size()<<" "<<ratematfixedvector.size()<<" "<<output->size<<" "<<numberoffreeparameters<<endl;
									if (discretechosenmodel==1) { //one rate
										message+="\n    q_";
										message+=i;
										message+="_";
										message+=j;
										message+=" = ";		
										message+=gsl_vector_get(output,0);
										message+=" +/- ";
										message+=gsl_vector_get(output,1);
										vectorposition=2;
										if (gsl_vector_get(output,0)<0.00000001 && !nonnegvariables) {
											message+="  Warning: an estimate near zero sometimes makes estimating other parameters, and therefore the lnL, very imprecise. Play with numopt or the model";
										}
									}
									if (discretechosenmodel==2) { //rev
										if (i<j) {
											message+="\n    q_";
											message+=i;
											message+="_";
											message+=j;
											message+=" = ";									
											message+=gsl_vector_get(output,vectorposition);
											message+=" +/- ";
											message+=gsl_vector_get(output,vectorposition+numberoffreeparameters);
										//since symmetric
											message+="\n    q_";
											message+=j;
											message+="_";
											message+=i;
											message+=" = ";									
											message+=gsl_vector_get(output,vectorposition);
											message+=" +/- ";
											message+=gsl_vector_get(output,vectorposition+numberoffreeparameters);
											if (gsl_vector_get(output,vectorposition)<0.00000001 && !nonnegvariables) {
												message+="  Warning: an estimate near zero sometimes makes estimating other parameters, and therefore the lnL, very imprecise. Play with numopt or the model";
											}										
											vectorposition++;
											
										}
									}
									if (discretechosenmodel==3) { //nonrev
										message+="\n    q_";
										message+=i;
										message+="_";
										message+=j;
										message+=" = ";									
										message+=gsl_vector_get(output,vectorposition);
										message+=" +/- ";
										message+=gsl_vector_get(output,vectorposition+numberoffreeparameters);
										if (gsl_vector_get(output,vectorposition)<0.00000001 && !nonnegvariables) {
											message+="  Warning: an estimate near zero sometimes makes estimating other parameters, and therefore the lnL, very imprecise. Play with numopt or the model";
										}									
										vectorposition++;
									}
									if (discretechosenmodel==4) { //user
										message+="\n    q_";
										message+=i;
										message+="_";
										message+=j;
										message+=" = ";	
										if (ratematassignvector[vectorposition]>=0) { //means there's an assigned rate
											message+=ratematfixedvector[(ratematassignvector[vectorposition])];
											message+=" FIXED";
										}
										else {
											position=-1*(1+ratematassignvector[vectorposition]);
											message+=gsl_vector_get(output,position);
											message+=" +/- ";
											assert((position+numberoffreeparameters)<output->size);
											message+=gsl_vector_get(output,position+numberoffreeparameters);	
											if (gsl_vector_get(output,position)<0.00000001 && !nonnegvariables) {
												message+="  Warning: an estimate near zero sometimes makes estimating other parameters, and therefore the lnL, very imprecise. Play with numopt or the model";
											}
											
										}
										vectorposition++;
									}
								}
							}
						}
						if (discretechosenmodel==4) {
							vectorposition=position+1; //since we care about position in the output vector, which just has variable parameters.
						}
						if (discretechosenstatefreqmodel==1) {
							message+="\n  Statefreqs: uniform";
						}
						else if (discretechosenstatefreqmodel==2) {
							message+="\n  Statefreqs: empirical";
						}
						else if (discretechosenstatefreqmodel==3) {
							message+="\n  Statefreqs: equilibrium";
						}
						else if (discretechosenstatefreqmodel==4) {
							message+="\n  Statefreqs: optimized";
						}					
						else if (discretechosenstatefreqmodel==5) {
							message+="\n  Statefreqs: user";
						}
						
						for (int i=0; i<localnumbercharstates; i++) { //do ancestralstatevector for freqs
							message+="\n    P(";
							message+=i;
							message+=") = ";							
						//cout<<i<<" "<<vectorposition<<" "<<optimaldiscretecharstatefreq->size<<" "<<output->size<<" "<<numberoffreeparameters<<endl;
							if (discretechosenstatefreqmodel==4) {
								if (i<(localnumbercharstates-1)) {
									message+=gsl_vector_get(output,vectorposition);
									message+=" +/- ";
									message+=gsl_vector_get(output,vectorposition+numberoffreeparameters);
									vectorposition++;
								}
								else { //last number must be 1-sum(other states)
									double frequencysum=0.0;
									for (int j=0; j<i;j++) {
										frequencysum+=gsl_vector_get(optimaldiscretecharstatefreq,j);
									}
									message+=1.0-frequencysum;
								}
							}
							else {
								message+=gsl_vector_get(optimaldiscretecharstatefreq,i);
							}
						}
						
						
						
						for (int n=0; n<localnumbercharstates; n++) {
							/*message+="\n\tP(";
							message+=n;
							message+=") = ";
							message+=gsl_vector_get(optimaldiscretecharstatefreq,n);*/
							char outputstring[14];
							sprintf(outputstring,"%E",gsl_vector_get(optimaldiscretecharstatefreq,n));
							tmessage+=outputstring;						
						//tmessage+=gsl_vector_get(optimaldiscretecharstatefreq,n);
							tmessage+="\t";
						}
						for (int i=0;i<localnumbercharstates;i++) {
							for (int j=0;j<localnumbercharstates;j++) {
								if (i!=j) {
									/*message+="\n\tq_";
									message+=i;
									message+="_";
									message+=j;
									message+=" = ";
									message+=gsl_matrix_get(optimaldiscretecharQmatrix,i,j);*/
									char outputstring[14];
									sprintf(outputstring,"%E",gsl_matrix_get(optimaldiscretecharQmatrix,i,j));
									tmessage+=outputstring;														
								//tmessage+=gsl_matrix_get(optimaldiscretecharQmatrix,i,j);
									tmessage+="\t";
								}
							}
						}
						if (tablef_open) {
							tablef<<tmessage;
						}
						PrintMessage();
						if (reconstruct) {
							NodePtr newroot=EstimateMLDiscreteCharJointAncestralStates(optimaldiscretecharQmatrix,optimaldiscretecharstatefreq,breaknum);
//(intrees.GetIthTree(chosentree-1)).SetRoot(newroot);
						}
						gsl_vector_free(output);
					}
				}
				chosentree=origchosentree;
				discretechosenchar=origdiscretechosenchar;
				if (tablef_open) {
					tablef.close();
				}
			}
			/*double rateA=0.0000000000000000001;
			while (rateA<100) {
				double neglnL=CalculateDiscreteLindy1(rateA);
				message="rateA=";
				message+=rateA;
				message+=" -lnL=";
				message+=neglnL;
				PrintMessage();
				rateA*=100.0;
			}
			double neglnL=CalculateDiscreteLindy1(0.022322);
			message="rateA=";
			message+=0.022322;
			message+=" -lnL=";
			message+=neglnL;
			PrintMessage();
			
			//BROWNIE my_fnA(maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput,trees, taxa, assumptions, characters);
			gsl_vector *optimalrateA=gsl_vector_calloc(3);
				//gsl_vector_memcpy(optimalrate,my_fn.OptimizeRateWithGivenTipVariance());
			gsl_vector_memcpy(optimalrateA,LindyGeneralOptimization(1));
			for(int i=0;i<optimalrateA->size;i++){
					message="i=";
				message+=i;
				message+=" val is ";
				message+=gsl_vector_get(optimalrateA,i);
				PrintMessage();
			}
			//BROWNIE my_fnB(maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput,trees, taxa, assumptions, characters);
			gsl_vector *optimalrateB=gsl_vector_calloc(5);
			gsl_vector_memcpy(optimalrateB,LindyGeneralOptimization(2));
			for(int i=0;i<optimalrateB->size;i++){
				message="i=";
				message+=i;
				message+=" val is ";
				message+=gsl_vector_get(optimalrateB,i);
				PrintMessage();
			}
			*/
			/*
			for (double rateA=0.000000001;rateA<1;rateA*=10) {
				gsl_matrix* ratematrixA=gsl_matrix_calloc(2,2);
				gsl_matrix_set(ratematrixA,0,0,0.0-rateA);
				gsl_matrix_set(ratematrixA,0,1,rateA);
				gsl_matrix_set(ratematrixA,1,0,rateA);
				gsl_matrix_set(ratematrixA,1,1,0.0-rateA);
				gsl_vector* basefreq=gsl_vector_calloc(2);
				gsl_vector_set(basefreq,0,0.5);
				gsl_vector_set(basefreq,1,0.5);
				double lnL=CalculateDiscreteCharLnL(ratematrixA,basefreq);
				cout<<"lnL="<<lnL<<" rate="<<rateA<<endl;
				NodePtr newroot=EstimateMLDiscreteCharJointAncestralStates(ratematrixA,basefreq,0);
				NodePtr newroot2=EstimateMLDiscreteCharJointAncestralStates(ratematrixA,basefreq,10);
			}*/
			
			
            break;
        }
		else if (token.Abbreviation("Breaknum")  ) {
            nxsstring breaknumchar;
            breaknumchar=GetFileName(token);
            breaknum=atoi(breaknumchar.c_str());
        }
        else if (token.Abbreviation("Treeloop") ) {
            nxsstring yesnotreeloop=GetFileName(token);
            if (yesnotreeloop[0] == 'n') {
                treeloop=false;
            }
            else {
                treeloop=true;
            }
        }		
        else if (token.Abbreviation("Charloop") ) {
            nxsstring yesnocharloop=GetFileName(token);
            if (yesnocharloop[0] == 'n') {
                charloop=false;
            }
            else {
                charloop=true;
            }
        }		
		else if (token.Abbreviation("Globalstates") ) {
            nxsstring yesnoglobalstates=GetFileName(token);
            if (yesnoglobalstates[0] == 'n') {
                globalstates=false;
            }
            else {
                globalstates=true;
            }
        }		
		
		else if (token.Abbreviation("ALlchar") ) {
            nxsstring yesnoallchar=GetFileName(token);
            if (yesnoallchar[0] == 'n' || yesnoallchar[0] == 'N') {
                allchar=false;
            }
            else {
                allchar=true;
            }
        }		
		else if (token.Abbreviation("Variable") ) {
            nxsstring yesnovarchar=GetFileName(token);
            if (yesnovarchar[0] == 'n' || yesnovarchar[0] == 'N') {
                variablecharonly=false;
            }
            else {
                variablecharonly=true;
            }
        }
		else if( token.Abbreviation("Replace") ) {
            nxsstring yesnoreplace=GetFileName(token);
            if (yesnoreplace[0] == 'n') {
                replacing=false;
            }
            else {
                replacing=true;
				appending=false;
            }
        }
        else if( token.Abbreviation("APpend") ) {
            nxsstring yesnoappend=GetFileName(token);
            if (yesnoappend[0] == 'n') {
                appending=false;
            }
            else {
                appending=true;
            }
        }
		
        else if (token.Abbreviation("Reconstruct") ) {
            nxsstring yesnoreconstruct=GetFileName(token);
            if (yesnoreconstruct[0] == 'n') {
                reconstruct=false;
            }
            else {
                reconstruct=true;
            }
        }	
		else if( token.Abbreviation("File") ) {
            tablefname = GetFileName(token);
            name_provided = true;
        }		
		else if( token.Abbreviation("Model") ) {
            nxsstring chosenmodelinput=GetFileName(token);
			//int numbercharstates=(discretecharacters->GetObsNumStates(discretechosenchar));
			int numberofrates=(localnumbercharstates*localnumbercharstates)-localnumbercharstates;
			int ntax=taxa->GetNumTaxonLabels();
            if (token.Abbreviation("Equal")) {
                    discretechosenmodel=1;
				message="You have chosen one rate for all discrete character transitions";
                PrintMessage();
            }
            else if(token.Abbreviation("Reversible")) {
				discretechosenmodel=2;
                message="You have chosen a time-reversible model: rates are free to vary, with the \nconstraint that forward and reverse rates for any two states are the same";
				if (ntax*10<(numberofrates/2) && !allchar) {
						message+="\n\nWARNING: You are trying to estimate ";
					message+=numberofrates/2;
					message+=" rates with only ";
					message+=ntax;
					message+=" taxa.";
				}
				else if (allchar && ntax*10*(discretecharacters->GetNChar())<(numberofrates/2)) {
					message+="\n\nWARNING: You are trying to estimate ";
					message+=numberofrates/2;
					message+=" rates with only ";
					message+=ntax;
					message+=" taxa and ";	
					message+=discretecharacters->GetNChar();
					message+=" characters simultaneously.";
				}
                PrintMessage();
            }
            else if((token.Abbreviation("Nonreversible")) || (token.Abbreviation("Irreversible"))) {
				discretechosenmodel=3;
                message="You have chosen a non-time-reversible model: rates are free to vary";
				if (ntax*10<(numberofrates) && !allchar) {
					message+="\n\nWARNING: You are trying to estimate ";
					message+=numberofrates;
					message+=" rates with only ";
					message+=ntax;
					message+=" taxa.";
				}
				else if (allchar && ntax*10*(discretecharacters->GetNChar())<(numberofrates)) {
					message+="\n\nWARNING: You are trying to estimate ";
					message+=numberofrates;
					message+=" rates with only ";
					message+=ntax;
					message+=" taxa and ";	
					message+=discretecharacters->GetNChar();
					message+=" characters simultaneously.";
				}
                PrintMessage();
            }
            else if(token.Abbreviation("User")) {
				discretechosenmodel=4;
                message="You have chosen a user-specified model";
				PrintMessage();
            }
			
            else {
                errormsg = "Unexpected option (";
                errormsg += chosenmodelinput;
                errormsg += ") encountered reading Model command";
                throw XNexus( errormsg);
            }
        }		
		else if( token.Abbreviation("Freq") ) {
            nxsstring chosenmodelinput=GetFileName(token);
            if (token.Abbreviation("Uniform")) {
				discretechosenstatefreqmodel=1;
				message="You have chosen equal frequencies for all states";
                PrintMessage();
            }
            else if(token.Abbreviation("EMpirical")) {
				discretechosenstatefreqmodel=2;
                message="You have chosen to use empirical state frequencies";
                PrintMessage();
            }
            else if(token.Abbreviation("EQuilibrium")) {
				discretechosenstatefreqmodel=3;
                message="You have chosen to use equilibrium state frequencies";
                PrintMessage();
            }
            else if(token.Abbreviation("Optimized")) {
				discretechosenstatefreqmodel=4;
                message="You have chosen to optimize state frequencies";
                PrintMessage();
            }
			else if(token.Abbreviation("User")) {
				discretechosenstatefreqmodel=5;
                message="You have chosen to use user-set state frequencies";
                PrintMessage();
            }
            else {
                errormsg = "Unexpected option (";
                errormsg += chosenmodelinput;
                errormsg += ") encountered reading Freq command";
                throw XNexus( errormsg);
            }
        }	
		else if( token.Abbreviation("Statevector") ) {
			if (debugmode) {
				cout<<"Now reading statevector"<<endl;
			}
            token.GetNextToken();
            token.GetNextToken(); //eat the equals sign
            vector<double> temporarystatevector;
            if (!token.Equals("(")) {
                errormsg="Expecting next token to be a \'(\' but instead got ";
                errormsg+=token.GetToken();
                throw XNexus( errormsg);
            }
            int inputcount=0;
            while (!token.Equals(")")) {
                nxsstring numbernexus;
                numbernexus=GetNumberOnly(token);
				if (debugmode) {
					cout<<"pushing back with "<<numbernexus<<endl;
				}
                if (numbernexus!=")") {
                    temporarystatevector.push_back(atof( numbernexus.c_str() ));
                    inputcount++;
                }
                else {
                    break;
                }
            }
			if (debugmode) {
				cout<<"finished with the pushback step"<<endl;
			}
			
			//int numbercharstates=(discretecharacters->GetObsNumStates(discretechosenchar));
            if (inputcount!=localnumbercharstates) {
                errormsg="You should have entered ";
				errormsg+=localnumbercharstates;
				errormsg+=" frequencies, you entered ";
                errormsg+=inputcount;
                throw XNexus( errormsg);
            }
            else {
				double sumoffreqs=0;
				for (int i=0; i<temporarystatevector.size(); i++) {
					sumoffreqs+=temporarystatevector[i];
				}
                userstatefreqvector.clear();
				message="Entering user frequencies of ( "; 
                for (int i=0; i<temporarystatevector.size(); i++) {
                    userstatefreqvector.push_back((temporarystatevector[i])/sumoffreqs);
					message+=userstatefreqvector[i];
					message+=" ";
                }
				message+=")";
				PrintMessage();
				
            }
        }		
        else if( token.Abbreviation("?") ) {
			donesomething=true;
            message="Usage: Discrete [model=] [freq=] [ratemat=] [statevector=] [treeloop=] [charloop=] [allchar=] [variable=] [reconstruct=] [breaknum=] [file=] [append=] [replace=] [globalstates=]\n\n";
			message+="This is a function to calculate the likelihood estimates of discrete character evolution parameters and scores for these models. This will allow\n";
			message+="you to do things like compare models with equal gain and loss rates with models which allow these to vary, evaluate models with a mixture of fixed\n";
			message+="and free rates, reconstruct the joint likelihood estimates of ancestral states at nodes and at various points within branches, and more.\n\n";
            message+="Available options:\n\n";
            message+="Keyword ------- Option type ----------------------------- Current setting -----";
			message+="\n Model          <string>                                  ";
			if (discretechosenmodel==1) {
					message+="Equal";
			}
			else if (discretechosenmodel==2) {
					message+="Rev";
			}
			else if (discretechosenmodel==3) {
					message+="NonRev";
			}
			else if (discretechosenmodel==4) {
					message+="User";
			}
			message+="\n Freq           <string>                                  ";
			if (discretechosenstatefreqmodel==1) {
				message+="Uniform";
			}
			else if (discretechosenstatefreqmodel==2) {
				message+="Empirical";
			}
			else if (discretechosenstatefreqmodel==3) {
				message+="Equilibrium";
			}
			else if (discretechosenstatefreqmodel==4) {
				message+="Optimized";
			}
			else if (discretechosenstatefreqmodel==5) {
				message+="User";
			}
			message+="\n RateMat        (<vector>)                                ";
			if (usermatrix.size()>0) {
				message+="( ";
				message+=usermatrix;
				message+=")";
			}
			else {
				message+="Unspecified";
			}			
			message+="\n StateVector    (<vector>)                                ";
			if ((userstatefreqvector.size())>1) {
				message+="( ";
				for (int i=0; i++; i<userstatefreqvector.size()) {
					message+=userstatefreqvector[i];
					message+=" ";
				}
				message+=")";
			}
			else {
				message+="Unspecified";
			}
            message+="\n Treeloop       No|Yes                                    *No";
            message+="\n Charloop       No|Yes                                    *No";
			message+="\n AllChar        No|Yes                                     ";
			if (allchar) {
				message+="Yes";
			}
			else {
				message+="No";
			}
			message+="\n Variable       No|Yes                                     ";
			if (variablecharonly) {
				message+="Yes";
			}
			else {
				message+="No";
			}
            message+="\n Reconstruct    No|Yes                                    *No";
            message+="\n Breaknum       <integer>                                 ";
			message+=breaknum;
			message+="\n File           <file name>                               *None";
			message+="\n Append         No|Yes                                    *Yes";
			message+="\n Replace        No|Yes                                    *No";
			message+="\n GlobalStates   No|Yes                                    *No";
            message+="\n                                                        *Option is nonpersistent\n";
			message+="\nModel: Allows you to specify whether to use a USER-specified model, a model where all rates are EQUAL, a REVersible model where q_ij=q_ji\n";
			message+="       for all states i and j but are otherwise free to vary, or an NONREVersible model where all rates can vary independently.";
			message+="\nRateMat: A vector containing information about rate parameters. Note that only off-diagonal entries should be included.\n";
			message+="         The model specification, except for the built-in types (equal, rev, nonrev), is grossly similar to PAUP's method for specifying\n";
			message+="         which rates are constrained to be equal, plus allows fixing of certain values. For example, the following rate matrix:\n";
			message+="                to->  0     1     2\n";
			message+="               from ------------------\n";
			message+="                 0 |  -     a    0.5\n";
			message+="                 1 |  b     -     c\n";
			message+="                 2 | 0.0    a     -\n";
			message+="         means that rate q01 (instantaneous rate going from state 0 to state 1), with value a, must also equal rate q21 but is otherwise unconstrained\n";
			message+="         (so they are both optimized, but forced to take the same value), while q20 is forced to a rate value of 0.0. We could have specified any\n";
			message+="         non-negative fixed value. Basically, all rates sharing a letter take the same optimized rate value, while those assigned a number have a fixed\n";
			message+="         rate value. To specify the above rate matrix, the command would be\n";
			message+="		     \"discrete model=user ratemat=(a 0.5 b c 0.0 a)\"\n";
			message+="         Letters are case-sensitive, so there are 52 (26*2) possible free rate parameters you can use in a user model.\n";
			message+="Freq: The probability of each state at the root can be based on the EMPIRICAL distribution at the tips, can be SET by the user\n";
			message+="      (using the statevector command), can be OPTIMIZEd as part of the model, can be set to EQUILIBRIUM frequencies (the\n";
			message+="      frequencies expected with the optimized rate matrix given infinitely-long branches), or can be set to be UNIFORM (equal).\n";
			message+="StateVector: Contains the user-specified probabilities of the ancestral states.\n";
			message+="             Example: \"discrete freq=user statevector=(0.4 0.6)\" for a binary trait\n";
			message+="Treeloop: Allows the analysis to be run across all the trees. A weighted average is returned (using tree weights such as posterior\n";
			message+="          probabilities or bootstrap frequencies if they are available) as well as values for the individual trees.\n";
			message+="Charloop: Allows the analysis to be run across all the characters individually.\n";
			message+="Allchar: Estimates the model parameters using all the characters simultaneously, not based on a single character.\n";
			message+="Variable: If true, performs a correction to correct rates for only examining variable characters.\n";
			message+="Reconstruct: Using the likelihood rate matrix and state frequencies, reconstructs the joint estimates of the ancestral states\n";
			message+="             at internal nodes and, optionally, along the branches. This uses the Pupko et al 2000 algorithm, which is fast but\n";
			message+="             only returns the estimated states, not the confidence in these states.\n";
			message+="Breaknum: Setting this value >0 allows the program to estimate the likeliest state at breaknum points along each branch. This can\n";
			message+="          be useful in estimating when on a branch a character changed. Note that this may underestimate the number of changes on a\n";
			message+="          branch (for example, on a very long branch which starts and ends in state 0, with high enough transition rates it may be\n";
			message+="          probable that the character has changed  multiple times on that branch, but at any given instant on that branch, the likeliest\n";
			message+="          state it will be in is state 0).\n";
			message+="File: Saves all output into a tab-delimited file\n";
			message+="Append: If the output file exists, appends to it rather than overwrites it. Turned on by default.\n";
			message+="Replace: If set to yes, if the output file already exists it will be quietly replaced.\n";
			message+="GlobalStates: If no, the number of character states assumed for each character is the maximum number of observed states for just that one character.\n";
			message+="              If yes, the number of states for each character is the maximum number of states observed for any character, even if the observed character\n";
			message+="              is lacking some of those states. This is useful if, for example, you have simulated a three state character on the tree for parametric\n";
			message+="              bootstrapping but some of the simulations result in characters with just states 0 and 1. Globalstates will automatically be set to yes\n";
			message+="              if allchar=y";
            PrintMessage();
        }
		else if( token.Abbreviation("Ratemat") ) {
            token.GetNextToken();
            token.GetNextToken(); //eat the equals sign
            vector<double> temporaryratematfixedvector; //Will be a vector containing JUST the fixed values
			vector<int> temporaryratematassignvector; //Will be a vector containing ints corresponding to "pointers" to either  fixed or variable values. If entries are non-negative,
													//they point to entries in temporaryratematfixedvector (i.e., value of 2 means the rate is whatever is stored at temporaryratematfixedvector[2])
													//negative values point to entries in a yet-to-be created vector of numbers to vary.
			string tempfreerateletterstring;		//Allows mapping of letters on input to negative values in temporaryratematassignvector vector. New letters are appended, old ones are looked up
			usermatrix=""; //just a string to store the description
            if (!token.Equals("(")) {
                errormsg="Expecting next token to be a \'(\' but instead got ";
                errormsg+=token.GetToken();
                throw XNexus( errormsg);
            }
            int inputcount=0;
            while (!token.Equals(")")) {
                nxsstring nextitem;
                nextitem=GetNumberOnly(token);
                if (nextitem!=")") {
					usermatrix+=nextitem;
					usermatrix+=" ";
					if (isalpha(nextitem[0])) { //Is a letter -- means that parameter is free to vary, but has same value as other rates with that value
						string::size_type loc = tempfreerateletterstring.find( nextitem[0], 0 );
						if( loc != string::npos ) {
							temporaryratematassignvector.push_back(-1*(loc+1));
						} else {
							tempfreerateletterstring.append(1,nextitem[0]);
							temporaryratematassignvector.push_back(-1*(tempfreerateletterstring.size()));         
						}
					}
					else { //is a number, which means it's a fixed value
						temporaryratematassignvector.push_back(temporaryratematfixedvector.size());
						temporaryratematfixedvector.push_back(atof( nextitem.c_str()));
					}
                    inputcount++;
                }
                else {
                    break;
                }
            }
			//int numbercharstates=(discretecharacters->GetObsNumStates(discretechosenchar));
			int expectednumberofentries=(localnumbercharstates*localnumbercharstates)-localnumbercharstates;
            if (inputcount!=expectednumberofentries) {
                errormsg="You should have entered ";
				errormsg+=expectednumberofentries;
				errormsg+=" values (the current character has ";
				errormsg+=localnumbercharstates;
				errormsg+="and so ";
				errormsg+=expectednumberofentries;
				errormsg+=" off-diagonal rates, you entered ";
                errormsg+=inputcount;
                throw XNexus( errormsg);
            }
            else {
				freerateletterstring=tempfreerateletterstring;
				ratematfixedvector.swap(temporaryratematfixedvector);
				ratematassignvector.swap(temporaryratematassignvector);
				if (debugmode) {
					cout<<"usermatrix is "<<usermatrix<<endl;
					cout<<"freerateletterstring is "<<freerateletterstring<<endl;
					cout<<"ratematfixedvector = ( ";
					for (int i=0;i<ratematfixedvector.size();i++) {
						cout<<ratematfixedvector[i]<<" ";
					}
					cout<<")"<<endl;
					cout<<"ratematassignvector = ( ";
					for (int i=0;i<ratematassignvector.size();i++) {
						cout<<ratematassignvector[i]<<" ";
					}
					cout<<")"<<endl;
					
				}
            }
        }
		
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading Discrete command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
    if (debugmode) {
    	cout<<endl<<"Now leaving Handle Discrete"<<endl;
    }
}

//Gets user commands for simulating characters on the tree and outputting them to a file
void BROWNIE::HandleSimulateCharacters( NexusToken& token ) {
	bool donenothing=true;
	int n=100;
	bool simtreeloop=false;
	int chartype=0; //0=discrete, 1=continuous
	nxsstring outputfilename="SimulatedChars.nex";
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (donenothing) {
				if (chartype==0 && optimaldiscretecharQmatrix->size1<2) {
					errormsg="Error: You must first input or optimize a model using the Discrete or Opt command";
					throw XNexus (errormsg);
				}
                SimulateCharacters(n,chartype,outputfilename,simtreeloop);
				//And  get output
            }
            break;
        }
        else if( token.Abbreviation("N") ) {
            nxsstring numbernexus = GetNumber(token);
            n=atoi( numbernexus.c_str() ); //convert to int
            message="You have chosen to simulate ";
            message+=n;
			message+=" characters";
            PrintMessage();
            if (n<1) {
                errormsg = "Error: must select a number greater than zero";
                n=100;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
		else if( token.Abbreviation("File") ) {
            outputfilename = GetFileName(token);
			message="Setting output file to ";
			message+=outputfilename;
			PrintMessage();
        }	
		else if (token.Abbreviation("Treeloop") ) {
            nxsstring yesnotreeloop=GetFileName(token);
            if (yesnotreeloop[0] == 'y' || yesnotreeloop[0] == 'Y') {
                simtreeloop=true;
				message="Setting to loop over trees";
				PrintMessage();
            }
            else {
                simtreeloop=false;
				message="Not looping over trees";
				PrintMessage();
            }
        }		
		else if (token.Abbreviation("Chartype") ) {
            nxsstring inputchartype=GetFileName(token);
            if (inputchartype[0] == 'd' || inputchartype[0] == 'D') {
                chartype=0; //discrete char
				message="Simulating discrete characters";
				PrintMessage();
            }
            else {
                chartype=1; //continuous char
				message="WARNING: Simulating continuous characters does not work yet";
				PrintMessage();
            }
        }				
        else if( token.Abbreviation("?") ) {
            donenothing=false;
            message="Usage: Simulate n=<integer> chartype=<discrete/continuous> treeloop=<yes/no> file=<output file>\n\n";
            message+="Simulates n discrete or continuous characters using the last optimized (or user-set) character model\n";
			message+="and saves these into a nexus file, along with the tree used to generate them. If treeloop=yes, trees\n";
			message+="are sampled for use in simulation based on their proportion of the total tree weight.\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --\n";
            message+="n            <integer>                            100*\n";
			message+="chartype     <discrete | continuous>              Discrete*\n";
			message+="file         <output file name>                   SimulatedChars.nex*\n";
			message+="treeloop     <yes | no>                           No*\n\n";
			message+="                                                *means option is nonpersistent";                
            message+="\n\n";
            PrintMessage();
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading SIMULATE command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}


//Simulates discrete or continuous characters and puts them into a nexus file, suitable for later reading by Brownie
void BROWNIE::SimulateCharacters(int n, int chartype, nxsstring outputfilename, bool treeloop) {
	ofstream simulationf;
	bool simulationf_open=false;
	bool exists = FileExists( outputfilename.c_str() );
	bool userok = true;
	if( exists && !UserSaysOk( "Ok to replace?", "Simulated character output file specified already exists" ) )
		userok = false;
	if( userok && !simulationf_open) {
		simulationf_open = true;
		simulationf.open( outputfilename.c_str() );
	}
	if( exists && userok ) {
		message = "\nReplacing simulation output file ";
		message += outputfilename;
		PrintMessage();

	}
	else if( userok ) {
		message = "\nSimulation output file ";
		message += outputfilename;
		message += " opened";
		PrintMessage();

	}
	else {
		errormsg = "Aborting the simulation so as not to overwrite the file.\n";
		throw XNexus( errormsg);
	}
	ProgressBar(n);
	int ntax=taxa->GetNumTaxonLabels();
	vector<nxsstring> charactermatrixvector;
	vector<double> startingfreqcumulativevector;
	int ntrees=trees->GetNumTrees();
	simulationf<<"#nexus\n";
	if(chartype==0) {
		simulationf<<"[current best rate matrix = \n";
	//int numbercharstates=optimaldiscretecharQmatrix->size1;
		localnumbercharstates=optimaldiscretecharQmatrix->size1;
		for (int i=0;i<localnumbercharstates;i++) {
			for (int j=0;j<localnumbercharstates;j++) {
				if (i!=j) {
					simulationf<<gsl_matrix_get(optimaldiscretecharQmatrix,i,j);
					simulationf<<"\t";
				}
				else {
					simulationf<<"-\t";
				}
			}
			simulationf<<"\n";
		}
		double frequencysum=0.0;
		for (int i=0; i<localnumbercharstates; i++) { //do ancestralstatevector for freqs
			simulationf<<"\n    P(";
			simulationf<<i;
			simulationf<<") = ";							
			if (i<(localnumbercharstates-1)) {
				frequencysum+=gsl_vector_get(optimaldiscretecharstatefreq,i);
				simulationf<<gsl_vector_get(optimaldiscretecharstatefreq,i);
				startingfreqcumulativevector.push_back(frequencysum);
				
			}
			else {
				startingfreqcumulativevector.push_back(1);
				simulationf<<1.0-frequencysum;
			}
		}
		simulationf<<"\n]\n\n";
	}
	simulationf<<"begin taxa;\ndimensions ntax="<<ntax<<";\ntaxlabels\n";
	for (int i=0;i<ntax;i++) {
		simulationf<<GetTaxonLabel(i)<<endl;
		nxsstring newstringforvector="";
		newstringforvector+=GetTaxonLabel(i);
		newstringforvector+=" ";
		charactermatrixvector.push_back(newstringforvector);
	}
	simulationf<<"\n;\nend;\n\n";
	simulationf<<"begin trees;\n";
	vector<double> treeweightvector;
	double totaltreeweight=0;
	for (int i=0;i<ntrees;i++) {
		simulationf<<"tree "<<trees->GetTreeName(i);
		simulationf<<" = [&W "<<trees->GetTreeWeight(i)<<" ] ";
		Tree t=intrees.GetIthTree(i);
		totaltreeweight+=trees->GetTreeWeight(i);
		treeweightvector.push_back(totaltreeweight);
		t.Write(simulationf);
		simulationf<<"\n";
	}
	for (int i=0;i<ntrees;i++) {
		treeweightvector[i]/=totaltreeweight; //standardize tree weights
	}
	simulationf<<"end;\n\n";
	simulationf<<"begin characters;\ndimensions nchar="<<n<<" ntax="<<ntax<<";\nformat  datatype=";
	if (chartype==0) {
		simulationf<<"standard";
	}
	else if (chartype==1) {
		simulationf<<"continuous";
	}
	simulationf<<";\nmatrix\n";
	for (int i=0;i<n;i++) {
		int oldchosentree=chosentree;
		//gsl_vector *newtips=gsl_vector_calloc(ntax);
		if (treeloop) {
			double desiredweight=gsl_ran_flat (r,0.0,1.0);
			for (int i=0;i<ntrees;i++) {
				if (desiredweight<treeweightvector[i]) {
					chosentree=i+1;
					break;
				}
			}			
		}
		Tree T=intrees.GetIthTree(chosentree-1);
				//Simulate up tree
		double desiredstartingfreq=gsl_ran_flat (r,0.0,1.0);
		int startingchar=0;
		if (chartype==0) {
			for (startingchar=0;startingchar<localnumbercharstates;startingchar++) {
				if (desiredstartingfreq<startingfreqcumulativevector[startingchar]) {
					break;
				}
			}	
		}
		//we've chosen a starting char
		if(chartype==0) {
			PreorderIterator <Node> m (T.GetRoot()); //Goes from root up
			NodePtr currentnode = m.begin();
			while (currentnode)
			{
				if (currentnode==T.GetRoot() ) {
					currentnode->SetLabelNumber(startingchar); //just use label numbers to store ancestral states
				}
				if (currentnode!=T.GetRoot() ) {
					int startstate=(currentnode->GetAnc())->GetLabelNumber();
					int nextstate=0;
					gsl_matrix * Pmatrix=ComputeTransitionProb(optimaldiscretecharQmatrix,currentnode->GetEdgeLength());
					vector<double> Pancstatetopossiblenext;
					double cumulativeP=0;
					for (int endstate=0;endstate<localnumbercharstates;endstate++) {
						cumulativeP+=gsl_matrix_get(Pmatrix,startstate,endstate);
						Pancstatetopossiblenext.push_back(cumulativeP);
					}
					double randomprob=gsl_ran_flat(r,0.0,1.0);
					for (nextstate=0;nextstate<localnumbercharstates;nextstate++) {
						if (randomprob<Pancstatetopossiblenext[nextstate]) {
							break;
						}
					}
					if (currentnode->IsLeaf()) {
						nxsstring newstate="";
						newstate+=nextstate;
						charactermatrixvector[taxa->FindTaxon(currentnode->GetLabel())]+=newstate;
					//gsl_vector_set(newtips,taxa->FindTaxon(currentnode->GetLabel()),nextstate);
					}
					else {
						currentnode->SetLabelNumber(nextstate);
					}
					gsl_matrix_free(Pmatrix);
				}
				currentnode=m.next();
			}
		}
		if (chartype==1) {
			if ((oldchosentree!=chosentree) || (i==0)) { //Either a first run or a new tree, so have to recalculate VCV and expected values
				GetOptimalVCVAndTraitsContinuous();
			}
			gsl_vector *tipsfromthissim=gsl_vector_calloc(ntax);
			if (debugmode) {
				cout<<"optimal VCV = "<<endl;
				PrintMatrix(optimalVCV);
			}
			tipsfromthissim=SimulateTips(optimalVCV, 1.0, optimalTraitMeans);
			for (int taxonpos=0;taxonpos<ntax;taxonpos++) {
				charactermatrixvector[taxonpos]+=gsl_vector_get(tipsfromthissim,taxonpos);
				charactermatrixvector[taxonpos]+="\t";
			}
		}
		
		
		
		//for (int taxon=0;taxon<ntax;taxon++) {
		//	charactermatrixvector[taxon]+=gsl_vector_get(newtips,taxon);
		//}
		//gsl_vector_free(newtips);
		chosentree=oldchosentree;
		ProgressBar(0);
	}
	for (int taxon=0;taxon<ntax;taxon++) {
		simulationf<<charactermatrixvector[taxon]<<"\n";
	}
	simulationf<<";\nend;\n\n";
	simulationf.close();
}


void BROWNIE::FindFixedDiscreteModel() {
	
	
}




/**Make combined VCV
*/
void BROWNIE::MakeCombinedVCV(gsl_matrix *VCVcombined, gsl_matrix *VCVtoadd, int ntaxprocessed)
{
    int ntaxtoadd=VCVtoadd->size1;
    gsl_matrix_view newview= gsl_matrix_submatrix(VCVcombined,ntaxprocessed,ntaxprocessed,ntaxtoadd,ntaxtoadd);
    gsl_matrix_memcpy(&newview.matrix,VCVtoadd);
}

//Simulate continuous tip values
gsl_vector* BROWNIE::SimulateTips(gsl_matrix * VCV, double rate, gsl_vector *MeanValues)
{
    //Code inspired by John Burkardt, also based on code from Handbook of Simulation: Principles, Methodology, Advances, Applications, and Practice, Jerry Banks, ed. 1998.
	//Code later changed to use ideas from http://www.mail-archive.com/help-gsl@gnu.org/msg00631.html by Ralph dos Santos Silva
    int ntax=VCV->size1;
	//cout<<"Now simulating tips\n";
	//PrintMatrix(VCV);
    gsl_vector *newtips=gsl_vector_calloc(ntax);
    gsl_matrix *A=gsl_matrix_calloc(ntax,ntax);
    int CopyResult= gsl_matrix_memcpy(A, VCV);
	gsl_matrix_scale (A,rate);
    int CholResult= gsl_linalg_cholesky_decomp( A); //A has LU in upper and lower diagonals
	gsl_vector *randomvect=gsl_vector_calloc(ntax);
	for (int i=0;i<ntax;i++) {
		gsl_vector_set(randomvect,i,gsl_ran_ugaussian(r));
	}
	gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, A, randomvect);
	
	//gsl_blas_dgemv (CblasNoTrans,1, A, randomvect,0, newtips); //Old method
	for (int i=0;i<ntax;i++) {
		//gsl_vector_set(newtips,i,gsl_vector_get(newtips,i)+gsl_vector_get(MeanValues,i));
		gsl_vector_set(newtips,i,gsl_vector_get(randomvect,i)+gsl_vector_get(MeanValues,i));
	}
	gsl_matrix_free(A);
	gsl_vector_free(randomvect);
	/*if(1==1) {
		gsl_vector *resid=gsl_vector_calloc(ntax);
		gsl_vector_memcpy(resid,newtips);
		double newancstate=GetAncestralState(DeleteStem(GetVCV(globalchosentaxset)),resid);
		gsl_vector_add_constant(resid, -1.0*newancstate);
		double newrate=EstimateRate(DeleteStem(GetVCV(globalchosentaxset)),resid);
		cout<<"VCV rate = "<<(gsl_matrix_get(VCV,0,0))/(gsl_matrix_get(GetVCV(globalchosentaxset),0,0))<<endl;
		cout<<"newrate="<<newrate<<endl<<"ancstate = "<<newancstate<<endl;
		cout<<"newtips=( ";
		for (int i=0;i<ntax;i++) {
			cout<<gsl_vector_get(newtips,i)<<" ";
		}
		cout<<")\n";
		
		cout<<"mean values vect=( ";
		for (int i=0;i<ntax;i++) {
			cout<<gsl_vector_get(MeanValues,i)<<" ";
		}
		cout<<")\n\n";
	}*/
    return newtips;

}

void BROWNIE::GetOptimalVCVAndTraitsContinuous()
{
	nxsstring chosentaxset=globalchosentaxset;
	if (chosenmodel==1 && tipvariancetype!=1) {
		int ntax=GetVCV(chosentaxset)->size1;
		double rate=gsl_vector_get(optimalvaluescontinuouschar,0);
	/*	if(1==1) {
			cout<<"original rate = "<<rate<<endl;
		}*/
		gsl_matrix *VCV=gsl_matrix_calloc(ntax,ntax);
		gsl_vector * tipvariance=gsl_vector_calloc(ntax);
		gsl_vector * observedtips=gsl_vector_calloc(ntax);
		observedtips=GetTipValues(chosentaxset,chosenchar);
		if (tipvariancetype==2) {
			tipvariance=GetTipValues(chosentaxset,chosenchar+1);
		}				
		VCV=DeleteStem(GetVCV(chosentaxset));
		gsl_matrix *RateTimesVCV=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix_memcpy(RateTimesVCV, VCV);
		gsl_matrix_scale(RateTimesVCV,rate);
		optimalVCV=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix_memcpy(optimalVCV,AddTipVarianceVectorToRateTimesVCV(DeleteStem(RateTimesVCV),tipvariance));
		//gsl_vector_free(optimalTraitMeans);
		optimalTraitMeans=gsl_vector_calloc(ntax);
		double ancestralstate=GetAncestralState(optimalVCV,observedtips);
		for (int taxon=0;taxon<ntax;taxon++) {
			gsl_vector_set(optimalTraitMeans,taxon,ancestralstate);
		}
		gsl_vector_free(observedtips);
		gsl_matrix_free(VCV);
		gsl_vector_free(tipvariance);
		gsl_matrix_free(RateTimesVCV);
	}
	else if (chosenmodel==5 && tipvariancetype!=1) { //BMS
		int ntax=GetVCV(chosentaxset)->size1;
		gsl_matrix * Matrix0=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix1=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix2=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix3=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix4=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix5=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix6=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix7=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix8=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix9=gsl_matrix_calloc(ntax,ntax);				
		Matrix0=GetVCVforOneModel(chosentaxset,0);
		Matrix1=GetVCVforOneModel(chosentaxset,1);
		Matrix2=GetVCVforOneModel(chosentaxset,2);
		Matrix3=GetVCVforOneModel(chosentaxset,3);
		Matrix4=GetVCVforOneModel(chosentaxset,4);
		Matrix5=GetVCVforOneModel(chosentaxset,5);
		Matrix6=GetVCVforOneModel(chosentaxset,6);
		Matrix7=GetVCVforOneModel(chosentaxset,7);
		Matrix8=GetVCVforOneModel(chosentaxset,8);
		Matrix9=GetVCVforOneModel(chosentaxset,9);
		gsl_matrix *VCVtotal=gsl_matrix_calloc(ntax,ntax);
		gsl_vector *tipvariance=gsl_vector_calloc(ntax);
		if (tipvariancetype==2) {
			tipvariance=GetTipValues(chosentaxset,chosenchar+1);
		}				
		int numberofmodels=-2+(optimalvaluescontinuouschar->size); //Do minus 2 because we have the lnL as an entry
		
		gsl_matrix *RateTimesVCV=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix_memcpy(RateTimesVCV, Matrix0);
		gsl_matrix_scale(RateTimesVCV,gsl_vector_get(optimalvaluescontinuouschar,1));
		gsl_matrix_add(VCVtotal,RateTimesVCV);
		if (numberofmodels>1) {
			gsl_matrix_memcpy(RateTimesVCV, Matrix1);
			gsl_matrix_scale(RateTimesVCV,gsl_vector_get(optimalvaluescontinuouschar,2));
			gsl_matrix_add(VCVtotal,RateTimesVCV);
		}
		if (numberofmodels>2) {
			gsl_matrix_memcpy(RateTimesVCV, Matrix2);
			gsl_matrix_scale(RateTimesVCV,gsl_vector_get(optimalvaluescontinuouschar,3));
			gsl_matrix_add(VCVtotal,RateTimesVCV);
		}	
		if (numberofmodels>3) {
			gsl_matrix_memcpy(RateTimesVCV, Matrix3);
			gsl_matrix_scale(RateTimesVCV,gsl_vector_get(optimalvaluescontinuouschar,4));
			gsl_matrix_add(VCVtotal,RateTimesVCV);
		}		
		if (numberofmodels>4) {
			gsl_matrix_memcpy(RateTimesVCV, Matrix4);
			gsl_matrix_scale(RateTimesVCV,gsl_vector_get(optimalvaluescontinuouschar,5));
			gsl_matrix_add(VCVtotal,RateTimesVCV);
		}	
		if (numberofmodels>5) {
			gsl_matrix_memcpy(RateTimesVCV, Matrix5);
			gsl_matrix_scale(RateTimesVCV,gsl_vector_get(optimalvaluescontinuouschar,6));
			gsl_matrix_add(VCVtotal,RateTimesVCV);
		}	
		if (numberofmodels>6) {
			gsl_matrix_memcpy(RateTimesVCV, Matrix6);
			gsl_matrix_scale(RateTimesVCV,gsl_vector_get(optimalvaluescontinuouschar,7));
			gsl_matrix_add(VCVtotal,RateTimesVCV);
		}	
		if (numberofmodels>7) {
			gsl_matrix_memcpy(RateTimesVCV, Matrix7);
			gsl_matrix_scale(RateTimesVCV,gsl_vector_get(optimalvaluescontinuouschar,8));
			gsl_matrix_add(VCVtotal,RateTimesVCV);
		}	
		if (numberofmodels>8) {
			gsl_matrix_memcpy(RateTimesVCV, Matrix8);
			gsl_matrix_scale(RateTimesVCV,gsl_vector_get(optimalvaluescontinuouschar,9));
			gsl_matrix_add(VCVtotal,RateTimesVCV);
		}	
		if (numberofmodels>9) {
			gsl_matrix_memcpy(RateTimesVCV, Matrix9);
			gsl_matrix_scale(RateTimesVCV,gsl_vector_get(optimalvaluescontinuouschar,10));
			gsl_matrix_add(VCVtotal,RateTimesVCV);
		}	
		optimalVCV=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix_memcpy(optimalVCV,AddTipVarianceVectorToRateTimesVCV(DeleteStem(VCVtotal),tipvariance));
		optimalTraitMeans=gsl_vector_calloc(ntax);
		for (int taxon=0;taxon<ntax;taxon++) {
			gsl_vector_set(optimalTraitMeans,taxon,gsl_vector_get(optimalvaluescontinuouschar,0));
		}
		gsl_matrix_free(VCVtotal);
		gsl_matrix_free(RateTimesVCV);
		gsl_vector_free(tipvariance);
		gsl_matrix_free(Matrix0);
		gsl_matrix_free(Matrix1);
		gsl_matrix_free(Matrix2);
		gsl_matrix_free(Matrix3);
		gsl_matrix_free(Matrix4);
		gsl_matrix_free(Matrix5);
		gsl_matrix_free(Matrix6);
		gsl_matrix_free(Matrix7);
		gsl_matrix_free(Matrix8);
		gsl_matrix_free(Matrix9);		
	}
	else if (chosenmodel==6 && tipvariancetype!=1) { //BMC
		int ntax=GetVCV(chosentaxset)->size1;
		gsl_matrix * Matrix0=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix1=gsl_matrix_calloc(ntax,ntax);
		Matrix0=GetVCVforChangeNoChange(chosentaxset,false); //branches with no changes
		Matrix1=GetVCVforChangeNoChange(chosentaxset,true); //branches with changes
		gsl_matrix *VCVtotal=gsl_matrix_calloc(ntax,ntax);
		gsl_vector *tipvariance=gsl_vector_calloc(ntax);
		if (tipvariancetype==2) {
			tipvariance=GetTipValues(chosentaxset,chosenchar+1);
		}						
		gsl_matrix *RateTimesVCV=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix_memcpy(RateTimesVCV, Matrix0);
		gsl_matrix_scale(RateTimesVCV,gsl_vector_get(optimalvaluescontinuouschar,1));
		gsl_matrix_add(VCVtotal,RateTimesVCV);
		gsl_matrix_memcpy(RateTimesVCV, Matrix1);
		gsl_matrix_scale(RateTimesVCV,gsl_vector_get(optimalvaluescontinuouschar,2));
		gsl_matrix_add(VCVtotal,RateTimesVCV);
		optimalVCV=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix_memcpy(optimalVCV,AddTipVarianceVectorToRateTimesVCV(DeleteStem(VCVtotal),tipvariance));
		optimalTraitMeans=gsl_vector_calloc(ntax);
		for (int taxon=0;taxon<ntax;taxon++) {
			gsl_vector_set(optimalTraitMeans,taxon,gsl_vector_get(optimalvaluescontinuouschar,0));
		}
		gsl_matrix_free(VCVtotal);
		gsl_matrix_free(RateTimesVCV);
		gsl_vector_free(tipvariance);
		gsl_matrix_free(Matrix0);
		gsl_matrix_free(Matrix1);
		
	}
	else if(chosenmodel==12 && tipvariancetype!=1) {
		int ntax=GetVCV(chosentaxset)->size1;
		gsl_vector * Vector1=gsl_vector_calloc(ntax);
		gsl_vector * Vector2=gsl_vector_calloc(ntax);
		gsl_vector * tips=gsl_vector_calloc(ntax);
		gsl_vector * variance=gsl_vector_calloc(ntax);
		Vector1=GetTipValues(chosentaxset,chosenchar);
		tips=GetTipValues(chosentaxset,chosenchar);
		if (tipvariancetype==2) {
			variance=GetTipValues(chosentaxset,chosenchar+1);
		}		
		gsl_matrix * Matrix0=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix1=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix2=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix3=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix4=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix5=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix6=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix7=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix8=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix * Matrix9=gsl_matrix_calloc(ntax,ntax);		
		Matrix0=DeleteStem(GetVCV(chosentaxset));
		gsl_matrix_free (Matrix1);
		Matrix1=gsl_matrix_calloc(ntax,maxstartstops*ntax); //assumes that states change fewer than maxstartstops/2 times from root to any tip
		Matrix1=GetStartStopTimesforOneState(chosentaxset,0);
		//cout<<"Startstoptimes matrix 1"<<endl;
		//PrintMatrix(Matrix1);
		gsl_matrix_free (Matrix2);
		Matrix2=gsl_matrix_calloc(ntax,maxstartstops*ntax);
		Matrix2=GetStartStopTimesforOneState(chosentaxset,1);
		gsl_matrix_free (Matrix3);
		Matrix3=gsl_matrix_calloc(ntax,maxstartstops*ntax);
		Matrix3=GetStartStopTimesforOneState(chosentaxset,2);
		gsl_matrix_free (Matrix4);
		Matrix4=gsl_matrix_calloc(ntax,maxstartstops*ntax);
		Matrix4=GetStartStopTimesforOneState(chosentaxset,3);
		gsl_matrix_free (Matrix5);
		Matrix5=gsl_matrix_calloc(ntax,maxstartstops*ntax);
		Matrix5=GetStartStopTimesforOneState(chosentaxset,4);
		gsl_matrix_free (Matrix6);
		Matrix6=gsl_matrix_calloc(ntax,maxstartstops*ntax);
		Matrix6=GetStartStopTimesforOneState(chosentaxset,5);
		gsl_matrix_free (Matrix7);
		Matrix7=gsl_matrix_calloc(ntax,maxstartstops*ntax);
		Matrix7=GetStartStopTimesforOneState(chosentaxset,6);
		gsl_matrix_free (Matrix8);
		Matrix8=gsl_matrix_calloc(ntax,maxstartstops*ntax);
		Matrix8=GetStartStopTimesforOneState(chosentaxset,7);
		gsl_matrix_free (Matrix9);
		Matrix9=gsl_matrix_calloc(ntax,maxstartstops*ntax);
		Matrix9=GetStartStopTimesforOneState(chosentaxset,8);
		
		double rate=gsl_vector_get(optimalvaluescontinuouschar,0);
		double attraction=gsl_vector_get(optimalvaluescontinuouschar,1);
		int numberofmeans=-4+(optimalvaluescontinuouschar->size); //DO minus 4 here because we include the lnL
		double rootmean=gsl_vector_get(optimalvaluescontinuouschar,2);
		gsl_matrix *VCVtotal=gsl_matrix_calloc(ntax,ntax);
		gsl_vector *tipvariance=gsl_vector_calloc(ntax);
		gsl_vector *observedtips=gsl_vector_calloc(ntax);
		gsl_vector *expectedtips=gsl_vector_calloc(ntax);
		gsl_matrix *W_BK_A7=gsl_matrix_calloc(ntax,numberofmeans+1); //W matrix based on equation A7 of Butler and King
		gsl_vector_memcpy(observedtips,Vector1);
		gsl_vector_memcpy(tipvariance,Vector2);
			//	if (detailedoutput) {
			//		brownie.message="rate = ";
			//		brownie.message+=rate;
			//		brownie.message+=" attraction = ";
			//		brownie.message+=attraction;
			//		brownie.PrintMessage();
			//	}
		
			//roottotiptime calculation (and probably this OU model in general) assumes the taxa are coeval.
		double roottotiptime=gsl_matrix_get(Matrix0,0,0);	
		gsl_matrix * BranchingTimes=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix_memcpy(BranchingTimes, Matrix0);
		gsl_matrix * ScaledVCV=gsl_matrix_calloc(ntax,ntax);
		double exptonegalphaT=gsl_sf_exp(-1.0*attraction*roottotiptime);
		for (int rowtaxon=0;rowtaxon<ntax;rowtaxon++) {
			for (int coltaxon=0;coltaxon<ntax;coltaxon++) {
				gsl_matrix_set(ScaledVCV,rowtaxon,coltaxon,(0.5*rate/attraction)*(gsl_sf_exp(-2.0*attraction*(roottotiptime-gsl_matrix_get(BranchingTimes,rowtaxon,coltaxon))))*(1.0-(gsl_sf_exp(-2.0*attraction*(gsl_matrix_get(BranchingTimes,rowtaxon,coltaxon))))));
			}
			gsl_matrix_set(W_BK_A7,rowtaxon,0,exptonegalphaT);
			if (numberofmeans>0) { //Here's where we do Butler and King A7
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<maxstartstops*ntax) {
					runningtotal+=gsl_sf_exp(attraction*gsl_matrix_get(Matrix1,rowtaxon,chosencolumn)); //if we don't have entries, we'll be taking e^0-e^0=0
					chosencolumn++;
					runningtotal+=-1.0*gsl_sf_exp(attraction*gsl_matrix_get(Matrix1,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,1,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>1) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<maxstartstops*ntax) {
					runningtotal+=gsl_sf_exp(attraction*gsl_matrix_get(Matrix2,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*gsl_sf_exp(attraction*gsl_matrix_get(Matrix2,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,2,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>2) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<maxstartstops*ntax) {
					runningtotal+=gsl_sf_exp(attraction*gsl_matrix_get(Matrix3,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*gsl_sf_exp(attraction*gsl_matrix_get(Matrix3,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,3,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>3) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<maxstartstops*ntax) {
					runningtotal+=gsl_sf_exp(attraction*gsl_matrix_get(Matrix4,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*gsl_sf_exp(attraction*gsl_matrix_get(Matrix4,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,4,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>4) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<maxstartstops*ntax) {
					runningtotal+=gsl_sf_exp(attraction*gsl_matrix_get(Matrix5,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*gsl_sf_exp(attraction*gsl_matrix_get(Matrix5,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,5,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>5) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<maxstartstops*ntax) {
					runningtotal+=gsl_sf_exp(attraction*gsl_matrix_get(Matrix6,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*gsl_sf_exp(attraction*gsl_matrix_get(Matrix6,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,6,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>6) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<maxstartstops*ntax) {
					runningtotal+=gsl_sf_exp(attraction*gsl_matrix_get(Matrix7,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*gsl_sf_exp(attraction*gsl_matrix_get(Matrix7,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,7,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>7) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<maxstartstops*ntax) {
					runningtotal+=gsl_sf_exp(attraction*gsl_matrix_get(Matrix8,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*gsl_sf_exp(attraction*gsl_matrix_get(Matrix8,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,8,exptonegalphaT*runningtotal);
			}
		}
		
		gsl_matrix *VCVfinal=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix_memcpy(VCVfinal,AddTipVarianceVectorToRateTimesVCV(DeleteStem(ScaledVCV),tipvariance));
		gsl_vector *tipresiduals=gsl_vector_calloc(ntax);
		gsl_vector *tipexpectations=gsl_vector_calloc(ntax);
		gsl_vector * OUmeans=gsl_vector_calloc(numberofmeans+1);
		gsl_vector_set(OUmeans,0,rootmean);
		for (int position=1;position<=numberofmeans;position++) {
			gsl_vector_set(OUmeans,position,gsl_vector_get(optimalvaluescontinuouschar,position+2));
		}
		gsl_blas_dgemv (CblasNoTrans,1, W_BK_A7, OUmeans,0, tipexpectations); 
		optimalVCV=gsl_matrix_calloc(ntax,ntax); 
		gsl_matrix_memcpy(optimalVCV,VCVfinal);
		optimalTraitMeans=gsl_vector_calloc(ntax); 
		gsl_vector_memcpy(optimalTraitMeans,tipexpectations);
		cout<<"W matrix"<<endl;
		PrintMatrix(W_BK_A7);
		cout<<"\nVCVfinal"<<endl;
		PrintMatrix(VCVfinal);
		cout<<"\nOUmeans = (";
		for (int i=0;i<OUmeans->size;i++) {
			cout<<" "<<gsl_vector_get(OUmeans,i);
		}
		cout<<" )"<<endl;
		//calculate real optimized means
		gsl_matrix *Wtranspose=gsl_matrix_calloc(W_BK_A7->size2,W_BK_A7->size1);
		gsl_matrix_transpose_memcpy (Wtranspose,W_BK_A7);
		gsl_matrix *scaledVCVtilde=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix_memcpy(scaledVCVtilde,ScaledVCV);
		gsl_matrix_scale(scaledVCVtilde,1.0/rate);
		gsl_matrix *InversescaledVCVtildestart=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix *InversescaledVCVtilde=gsl_matrix_calloc(ntax,ntax);
		gsl_vector *A8OUmeans=gsl_vector_calloc(numberofmeans+1);
		gsl_permutation * p = gsl_permutation_alloc (ntax);
		int signum;
		gsl_matrix_memcpy (InversescaledVCVtildestart, scaledVCVtilde);
		gsl_linalg_LU_decomp (InversescaledVCVtildestart,p, &signum);
		gsl_linalg_LU_invert (InversescaledVCVtildestart,p, InversescaledVCVtilde);
		//OUmeans=(W'(V^-1)W)W'(V^-1)tips
		//VecA=(V^-1)tips
		//VecB=W'VecA
		//MatC=(V^-1)W
		//MatD=W'MatC
		//OUmeans=MatD*VecB
		gsl_vector *VecA=gsl_vector_calloc(ntax);
		gsl_vector *VecB=gsl_vector_calloc(numberofmeans+1);
		gsl_matrix *MatC=gsl_matrix_calloc(W_BK_A7->size1,W_BK_A7->size2);
		gsl_matrix *MatD=gsl_matrix_calloc(numberofmeans+1,numberofmeans+1);
		gsl_blas_dgemv (CblasNoTrans, 1.0, InversescaledVCVtilde, observedtips, 0.0, VecA);
		cout<<"Wtranspose row = "<<Wtranspose->size1<<" col = "<<Wtranspose->size2<<endl;
		cout<<"W_BK_A7 row = "<<W_BK_A7->size1<<" col = "<<W_BK_A7->size2<<endl;
		cout<<"VecA size = "<<VecA->size<<endl<<"VecB size = "<<VecB->size<<endl;
		gsl_blas_dgemv (CblasNoTrans, 1.0, Wtranspose, VecA, 0.0, VecB);
		gsl_blas_dgemm (CblasNoTrans,CblasNoTrans, 1.0, InversescaledVCVtilde, W_BK_A7, 0.0, MatC);
		gsl_blas_dgemm (CblasNoTrans,CblasNoTrans, 1.0, Wtranspose, MatC, 0.0, MatD);
		gsl_blas_dgemv (CblasNoTrans, 1.0, InversescaledVCVtilde, observedtips, 0.0, VecA);
		gsl_blas_dgemv (CblasNoTrans, 1.0, MatD, VecB, 0.0, A8OUmeans);
		cout<<"VecA"<<endl;
		PrintVector(VecA);
		cout<<"VecB"<<endl;
		PrintVector(VecB);
		cout<<"MatC"<<endl;
		PrintMatrix(MatC);
		cout<<"MatD"<<endl;
		PrintMatrix(MatD);
		
		cout<<"Recontructed OU means using Butler King A8\n   (";
		for (int i=0;i<A8OUmeans->size;i++) {
			cout<<" "<<gsl_vector_get(A8OUmeans,i);
		}
		cout<<" )"<<endl;
		gsl_vector_free(A8OUmeans);
		gsl_vector_free(VecA);
		gsl_vector_free(VecB);
		gsl_matrix_free(MatC);
		gsl_matrix_free(MatD);
		gsl_permutation_free (p);
		gsl_matrix_free(Wtranspose);
		gsl_matrix_free(scaledVCVtilde);
		gsl_matrix_free(InversescaledVCVtildestart);
		gsl_matrix_free(InversescaledVCVtilde);
		gsl_matrix_free (VCVfinal);
		gsl_vector_free (tipresiduals);
		gsl_vector_free (tipexpectations);
		gsl_vector_free (OUmeans);
		gsl_matrix_free(BranchingTimes);
		gsl_matrix_free(ScaledVCV);
		gsl_matrix_free(VCVtotal);
		gsl_vector_free(tipvariance);
		gsl_vector_free(observedtips);
		gsl_vector_free(expectedtips);
		gsl_matrix_free(W_BK_A7);		
		gsl_matrix_free(Matrix0);
		gsl_matrix_free(Matrix1);
		gsl_matrix_free(Matrix2);
		gsl_matrix_free(Matrix3);
		gsl_matrix_free(Matrix4);
		gsl_matrix_free(Matrix5);
		gsl_matrix_free(Matrix6);
		gsl_matrix_free(Matrix7);
		gsl_matrix_free(Matrix8);
		gsl_matrix_free(Matrix9);
		gsl_vector_free(Vector1);
		gsl_vector_free(Vector2);
		gsl_vector_free(tips);
		gsl_vector_free(variance);
	}
}


//Convert VCV matrix (n x n) to VCV vector (length n^2) and append it to a pre-existing vector
// if you don't have a pre-existing vector, create one: gsl_vector OrigVector(0,0);
//gsl_vector* BROWNIE::ConvertVCVMatrixToVector(gsl_matrix *VCV, gsl_vector *OrigVector)
//{
//    int OrigVectorLength=OrigVector->size;
//   int VCVMatrixArea=(VCV->size1)*(VCV->size2);
//    gsl_matrix *VCVcomb=gsl_matrix_alloc(0,0);
//    gsl_matrix *OrigMatrix=gsl_matrix_alloc(0,0);
//    gsl_vector *OutVector;
//   OutVector=gsl_vector_calloc(0);
//    OrigMatrix=ConvertVCVVectorToMatrix(OrigVector);
//    VCVcomb=MakeCombinedVCV(OrigMatrix,VCV);
//    int ntax=VCVcomb->size1;
//    for (int r=0;r<ntax;r++) {
//        for (int c=0;c<ntax;c++) {
//            OutVector.push_back(VCVcomb[r][c]);
//        }
//    }
//    return OutVector;
//}

//Convert VCV vector (length n^2) to VCV matrix (n x n)
//gsl_matrix* BROWNIE::ConvertVCVVectorToMatrix(gsl_vector *VCVvector)
//{
//    double ntaxsquared=VCVvector->size;
//    int ntax=int(sqrt(ntaxsquared));
//    gsl_matrix *VCV=gsl_matrix_calloc (ntax,ntax);
//    for (int r=0;r<ntax;r++) {
//        for (int c=0;c<ntax;c++) {
//            VCV[r][c]=VCVvector[(r*ntax)+c];
//        }
//    }
//    return VCV;
//}

    //Extract one VCV matrix from a VCV combined matrix that is represented as a vector
    //gsl_matrix * BROWNIE::ExtractMatrixFromVector(gsl_vector *VCVvector, int ntaxprocessed, int currentntax) {
    //    gsl_matrix *VCVentiremat=gsl_matrix_calloc(0,0);
    //    gsl_matrix *VCVextracted=gsl_matrix_calloc(currentntax,currentntax);
    //    VCVentiremat=ConvertVCVVectorToMatrix(VCVvector);
    //    for (int r=0;r<currentntax;r++) {
    //        for (int c=0;c<currentntax;c++) {
    //            VCVextracted[r][c]=VCVentiremat[ntaxprocessed+r][ntaxprocessed+c];
    //        }
    //    }
    //    return VCVextracted;
    //}


    /**Make combined tips
    */
    //void BROWNIE::MakeCombinedTips(gsl_vector *tipscombined, gsl_vector *tipstoadd, int ntaxprocessed)
    void BROWNIE::MakeCombinedTips(gsl_vector *tipscombined, gsl_vector *tipstoadd, int ntaxprocessed)
{
        int ntaxtoadd=tipstoadd->size;
        gsl_vector_view newview=gsl_vector_subvector(tipscombined,ntaxprocessed,ntaxtoadd);
        gsl_vector_memcpy(&newview.vector,tipstoadd);
}

//Returns log likelihood. If the VCV matrix includes other components (like tip variance), deal with these AND THE RATE first, and just pass a rate of 1 to this function.
double BROWNIE::GetLScore(gsl_matrix *VCV,gsl_vector *tipresid,double rate){
    int ntax=VCV->size1;
    double lscore;
    gsl_matrix *RateTimesVCV=gsl_matrix_calloc(ntax,ntax);
    //gsl_matrix RateTimesVCV(ntax,ntax,0);
    gsl_matrix *InverseRateTimesVCVstart=gsl_matrix_calloc(ntax,ntax);
    gsl_matrix *InverseRateTimesVCV=gsl_matrix_calloc(ntax,ntax);
    //gsl_matrix InverseRateTimesVCV(ntax,ntax,0);
    //long double DetRateTimesVCV;
    gsl_matrix_memcpy (RateTimesVCV, VCV);

/*    if(debugmode) {
        message="debugging line 957\n";
        for (int currentrow=0;currentrow<VCV->size1;currentrow++) {
            for (int currentcol=0;currentcol<VCV->size2;currentcol++) {
                message+=gsl_matrix_get(VCV,currentrow,currentcol);
                message+="\t";
            }
            message+="\n";
        }
        message+="\n\n RateTimesVCV\n";
        for (int currentrow=0;currentrow<RateTimesVCV->size1;currentrow++) {
            for (int currentcol=0;currentcol<RateTimesVCV->size2;currentcol++) {
                message+=gsl_matrix_get(RateTimesVCV,currentrow,currentcol);
                message+="\t";
            }
            message+="\n";
        }

        PrintMessage();
    } */


    gsl_matrix_scale (RateTimesVCV,rate);
/*    if(debugmode) {
        message="debugging line 980\n";
        message+="\nrate is ";
        message+=rate;
        for (int currentrow=0;currentrow<RateTimesVCV->size1;currentrow++) {
            for (int currentcol=0;currentcol<RateTimesVCV->size2;currentcol++) {
                message+=gsl_matrix_get(RateTimesVCV,currentrow,currentcol);
                message+="\t";
            }
            message+="\n";

        }

        PrintMessage();
    } */


    //RateTimesVCV=rate*VCV;
    //matrixsingular=TestSingularity(RateTimesVCV);
    //if (matrixsingular) {
    //    errormsg="Singular matrix (RateTimesVCV) during GetLScore";
    //    throw XNexus(errormsg);
    //}
/*    if (debugmode) {
        message="\nNow calculating inverse VCV for estimating lnL.\n";
        PrintMessage();
    } */
//USE THESE

gsl_permutation * p = gsl_permutation_alloc (ntax);
int signum;
gsl_matrix_memcpy (InverseRateTimesVCVstart, RateTimesVCV);
/*if(debugmode) {
    message="debugging\n";
    for (int currentrow=0;currentrow<RateTimesVCV->size1;currentrow++) {
        for (int currentcol=0;currentcol<RateTimesVCV->size2;currentcol++) {
            message+=gsl_matrix_get(RateTimesVCV,currentrow,currentcol);
            message+="\t";
        }
        message+="\n";
    }
    PrintMessage();
}
if(debugmode) {
    message="ntax is ";
    message+=ntax;
    message+=" and size of InverseRateTimesVCVstart is ";
    int irtvssize=InverseRateTimesVCVstart->size1;
    message+=irtvssize;
    int permsize=p->size;
    message+=" size of perm is ";
    message+=permsize;
    PrintMessage();
}*/
gsl_linalg_LU_decomp (InverseRateTimesVCVstart,p, &signum);
gsl_linalg_LU_invert (InverseRateTimesVCVstart,p, InverseRateTimesVCV);
//gsl_permutation_free (p);
//InverseRateTimesVCV=Inverse(RateTimesVCV);
/*if (debugmode) {
    message="\nFinished calculating inverse VCV for estimating rate.\n";
    PrintMessage();
    cout<<"RateTimesVCV\n";
    cout<<RateTimesVCV;
    cout<<"\n";
}

if(debugmode) {
    message="ntax is ";
    message+=ntax;
    message+=" and size of InverseRateTimesVCVstart is ";
    int irtvssize=InverseRateTimesVCVstart->size1;
    message+=irtvssize;
    int permsize=p->size;
    message+=" size of perm is ";
    message+=permsize;
    PrintMessage();
}
*/
gsl_matrix *tipsasmatrixprime=gsl_matrix_calloc(1,ntax);
//gsl_matrix tipsasmatrixprime(1,ntax,0);
for (int i=0; i<ntax; i++) {
    gsl_matrix_set(tipsasmatrixprime, 0, i, gsl_vector_get(tipresid,i));
    // tipsasmatrixprime[0][i]=tipresid[i];
}
gsl_vector *step1vect;
step1vect=gsl_vector_calloc(ntax);
//gsl_vector step1vect(ntax,0);
gsl_blas_dgemv (CblasNoTrans,1, InverseRateTimesVCV, tipresid,0, step1vect); ///TEST THIS
                                                                             //step1vect=MatrixTimesVector(InverseRateTimesVCV,tipresid);
gsl_vector *step2vect;
step2vect=gsl_vector_calloc(1);
//gsl_vector step2vect(1,0);
gsl_blas_dgemv (CblasNoTrans,1, tipsasmatrixprime, step1vect,0, step2vect); ///TEST THIS
                                                                            //step2vect=MatrixTimesVector(tipsasmatrixprime,step1vect);
                                                                            //DetRateTimesVCV=Determinant(RateTimesVCV);
                                                                            //lscore=-1*log((exp(-0.5*step2vect[0]))/(sqrt(DetRateTimesVCV*(pow((2*PI),ntax)))));
gsl_matrix *RateTimesVCVLU=gsl_matrix_calloc(ntax,ntax);
gsl_matrix_memcpy (RateTimesVCVLU, RateTimesVCV);
/*if(debugmode) {
    message="debugging\n";
    for (int currentrow=0;currentrow<RateTimesVCV->size1;currentrow++) {
        for (int currentcol=0;currentcol<RateTimesVCV->size2;currentcol++) {
            message+=gsl_matrix_get(RateTimesVCV,currentrow,currentcol);
            message+="\t";
        }
        message+="\n";
    }
    PrintMessage();
}*/
gsl_linalg_LU_decomp (RateTimesVCVLU,p, &signum);
lscore=0.5*gsl_vector_get(step2vect,0)+0.5*(gsl_linalg_LU_lndet(RateTimesVCVLU))+0.5*ntax*log(2*PI);
/*if (debugmode) {
    message="Likelihood score is ";
    message+=lscore;
    PrintMessage();
}*/
gsl_matrix_free(RateTimesVCV);
gsl_matrix_free(InverseRateTimesVCVstart);
gsl_matrix_free(InverseRateTimesVCV);
gsl_matrix_free(tipsasmatrixprime);
gsl_vector_free(step1vect);
gsl_vector_free(step2vect);
gsl_matrix_free(RateTimesVCVLU);
gsl_permutation_free (p);
return lscore; //-lnL actually
}




/** @method HandleRateTest
*
*
*/
void BROWNIE::HandleRateTest( NexusToken& token )
{
	citationarray[0]=true;
    bool noerror=true;
    //int seed=time(NULL);
    //int *seedptr=&seed;
    //srand(time(0));
    typedef std::set<nxsstring, std::less < nxsstring> > nxsstring_set;
    nxsstring tmessage;
    ofstream tablef;
    nxsstring tablefname;
    bool tablef_open=false;
    bool name_provided=false;
    nxsstring_set chosentaxSETS;
    //gsl_vector *chosentaxsetNTAXvector=gsl_vector_calloc(0);
    //gsl_vector chosentaxsetNTAXvector(0,0);
    set<int> chosentaxNTAX;
    map<nxsstring,int> chosentaxsetntaxmap;
    //set<gsl_vector> chosentaxCholVCVvectors;
    // gsl_vector vectorofCholVCVvectors;
    // map<int,gsl_vector> mapofCholVCVvectors;
    //gsl_vector *CholVCVvectorvector;
    // CholVCVvectorvector=gsl_vector_calloc(0);
    //gsl_vector CholVCVvectorvector(0,0);
    //gsl_vector *VCVvectorvector;
    //VCVvectorvector=gsl_vector_calloc(0);
    //gsl_vector VCVvectorvector(0,0);

    nxsstring chosentaxset;
    int repsnumber=0;
    int listedtaxsets=0;
    int ntaxcomb=0;
    bool treeloop=false;
    bool charloop=false;
    bool adequateinput=false;
    bool notquietmode=true;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (adequateinput==false) {
                message="Insufficient input: type \"ratetest ?\" for help";
                PrintMessage();
            }
            break;
        }
        else if( token.Equals("?") ) {
            adequateinput=true;
            message="Usage: RateTest taxset=<taxset 1> [taxset=<taxset 2>...] [options...]\n\n";
            message+="Performs a censored rate test for two or more subtrees.\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --";
            message+="\nTaxset       <taxset name>                        *None";
            message+="\nReps         <integer>                            *0";
            message+="\nTreeloop     No|Yes                               *No";
            message+="\nCharloop     No|Yes                               *No";
            message+="\nQuiet        No|Yes                               *No";
            message+="\nFile         <file name>                          *None";
            message+="\n                                                 *Option is nonpersistent\n\n";
            message+="This will compare the likelihood under a single rate Brownian motion\nmodel with the likelihood under a multiple parameter Brownian motion\nmodel using the 'censored' test using the current tree (use 'choose' to\nchange the current tree). If just one taxset is given, it will assign\none rate to the pruned tree containing just those taxa and another rate\nto the pruned tree containing all the other taxa. If multiple taxsets\nare given, it will assign one rate to each pruned tree and test this\nagainst the simple model where there is one rate for all these pruned\ntrees. \n\nNote that this method currently does not do much error checking: it is\nup to the user to check that the pruned trees represented by the taxsets\nare clades or paraphyletic with respect to the other pruned trees and\nthat no taxon is in multiple taxsets. If you are just entering one\ntaxset, and its pruned tree is a clade, you should have nothing to worry\nabout, even if its complement is paraphyletic.\n\nThe 'reps' statement is also optional. This tells the program how many\ntimes it should simulate data under the null model when doing parametric\nbootstrapping to test for significance using a likelihood ratio test.\nThe default is 0. When set to zero, it avoids parametric\nbootstrapping altogether.\n\nTreeloop and charloop tell the program whether to loop across all\ntrees and/or characters or just use the currently-selected ones.\nThe default for both is no.\nIf quiet=yes, only the summary is printed.\n\n\nYou may want to have logging activated to record the output from ratetest.";
            PrintMessage();
        }
        else if( token.Abbreviation("File") ) {
            tablefname = GetFileName(token);
            name_provided = true;
            bool exists = FileExists( tablefname.c_str() );
            bool userok = true;
            if( exists && !UserSaysOk( "Ok to replace?", "Ratetest output file specified already exists" ) )
                userok = false;
            if( userok ) {
                tablef_open = true;
                tablef.open( tablefname.c_str() );
            }

            if( exists && userok ) {
                message = "\nReplacing ratetest output file ";
                message += tablefname;
            }
            else if( userok ) {
                message = "\nRatetest output file ";
                message += tablefname;
                message += " opened";
            }
            else {
                errormsg = "Aborting the ratetest so as not to overwrite the file.\n";
                throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );

            }
            PrintMessage();

        }

        else if (token.Abbreviation("Taxset") ) {
            adequateinput=true;
            int ntax=0;
            if (trees->GetNumTrees()<1) {
                errormsg = "Error: No valid trees are loaded.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            chosentaxset=GetFileName(token);
            chosentaxSETS.insert(chosentaxset);
			
            listedtaxsets++;
            IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
            if (taxonlist.empty()) {
				errormsg="Error: Taxset ";
                errormsg+=chosentaxset.c_str();
                errormsg+=" does not exist.\nYou can define it using the taxset command.";
				assumptions->Report(cerr);
				if( logf_open )
					assumptions->Report(logf);		
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
IntSet::const_iterator xi;
            for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
                ntax++;
            }
            chosentaxsetntaxmap[chosentaxset]=ntax;
            ntaxcomb+=ntax;
        }
        else if (token.Abbreviation("Reps")  ) {
            nxsstring repschar;
            repschar=GetFileName(token);
            repsnumber=atoi(repschar.c_str());
        }
        else if (token.Abbreviation("TReeloop") ) {
            nxsstring yesnotreeloop=GetFileName(token);
            if (yesnotreeloop[0] == 'n') {
                treeloop=false;
            }
            else {
                treeloop=true;
            }
        }
        else if (token.Abbreviation("Quiet") ) {
            nxsstring yesnoquiet=GetFileName(token);
            if (yesnoquiet[0] == 'n') {
                notquietmode=true;
            }
            else {
                notquietmode=false;
            }
        }
        else if (token.Abbreviation("CHarloop") ) {
            nxsstring yesnocharloop=GetFileName(token);
            if (yesnocharloop[0] == 'n') {
                charloop=false;
            }
            else {
                charloop=true;
            }
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading RateTest command. Type \"RateTest ?\" for help.";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }


    }
    if (listedtaxsets>0) {
        if (listedtaxsets==1) {
            listedtaxsets++;
            nxsstring chosentaxset2="NOT";
            chosentaxset2+=chosentaxset; //so if one taxset is listed, its complement is automatically used.
            IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset2 );
            if (taxonlist.empty()) {
                errormsg= "Error: Taxset ";
                errormsg+=chosentaxset2.c_str();
                errormsg+=" does not exist.\nYou can define it using the taxset command.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            chosentaxSETS.insert(chosentaxset2);
            int ntax=0;
IntSet::const_iterator xi;
            for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
                ntax++;
            }
            chosentaxsetntaxmap[chosentaxset2]=ntax;
            ntaxcomb+=ntax;
        }



        if (tablef_open) {
            tmessage="Output in tab-delimited table:\nModel A is constrained to have one rate for all groups, Model B has one rate for each group.\n";
            tablef<<tmessage;
        }
        nxsstring_set chosentaxSETS;
        map<nxsstring, int>::const_iterator iter;
        int n=1;
        for (iter=chosentaxsetntaxmap.begin();iter!=chosentaxsetntaxmap.end();++iter) {
            //nxsstring_set::const_iterator cti;
            //for( cti = chosentaxSETS.begin(); cti != chosentaxSETS.end(); cti++ ) {
            tmessage="Taxset ";
            tmessage+=n;
            n++;
            tmessage+="=";
            tmessage+=iter->first;
            //message+=chosentaxSETS(*cti);
            tmessage+="\n";
            if (tablef_open) {
                tablef<<tmessage;
            }
            }
        tmessage="Tree\tTree weight\tTree name\tChar\t";
        for (int n=0; n<listedtaxsets; n++) {
            tmessage+="anc_";
            tmessage+=n+1;
            tmessage+="\trate_";
            tmessage+=n+1;
            tmessage+="\t -lnL_";
            tmessage+=n+1;
            tmessage+="\t";
        }
        tmessage+="rate_A\tparam_A\tparam_B\tAIC_A\tAIC_B\tAICc_A\tAICc_B\t -lnL_A\t -lnL_B\tAIC dif\tAICc diff\tchi p\tparam p\tchosen model under AIC, AICc, chi, param\n";
        if (tablef_open) {
            tablef<<tmessage;
        }
        int starttree, stoptree, startchar, stopchar;
        int originalchosentree=chosentree;
        int originalchosenchar=chosenchar;
        if (treeloop) {
            starttree=1;
            stoptree=trees->GetNumTrees();
        }
        else {
            starttree=chosentree;
            stoptree=chosentree;
        }
        if (charloop) {
            startchar=1;
            stopchar=continuouscharacters->GetNChar();
        }
        else {
            startchar=chosenchar;
            stopchar=chosenchar;
        }
        //Loop across trees and chars
        nxsstring summaryofresults="Summary of results\n(B/b = strong/weak support for multiple rate model)\nTree\tChar\tAIC\tAICc\tChi\tParam\n";
        if (repsnumber>0) {
            message="\nExpect a delay while doing parametric simulation.\nIf this is taking too long, use the option reps=0 when using the\nratetest command in the future. Note that parametric simulation\nis very important if you are using a p-value approach and do not\nhave many taxa (when the chi-square test is non-conservative).\n";
            PrintMessage();
        }
        double weighttotal=0;
        double weightedAIC1=0;
        double weightedAIC2=0;
        double weightedAICc1=0;
        double weightedAICc2=0;
        gsl_vector *weightedratevector;
        weightedratevector=gsl_vector_calloc(listedtaxsets+1);
        //gsl_vector weightedratevector(listedtaxsets+1,0); //So that the first rate is for the single rate model
        gsl_vector *weightedancstatevector;
        weightedancstatevector=gsl_vector_calloc(listedtaxsets+1);
        //gsl_vector weightedancstatevector(listedtaxsets+1,0); // Just for consistency with above
        double weightedchip=0;
        double weightedparamp=0;


        for (chosentree=starttree;chosentree<=stoptree;chosentree++) {
            int ntaxprocessedtreeloop=0;
            bool goodtree=true;
            //gsl_vector *CholVCVvectorvector;
            //CholVCVvectorvector=gsl_vector_calloc(0);
            //            gsl_vector CholVCVvectorvector(0,0);
            //gsl_vector *VCVvectorvector;
            //VCVvectorvector=gsl_vector_calloc(0);
            //            gsl_vector VCVvectorvector(0,0);
            gsl_matrix *VCVcomb=gsl_matrix_calloc(ntaxcomb,ntaxcomb);
            //gsl_matrix VCVcomb(0,0,0);
            double treeweight=trees->GetTreeWeight(chosentree-1);
            weighttotal+=treeweight;
            nxsstring treename=trees->GetTreeName(chosentree-1);
            map<nxsstring, int>::const_iterator iter;
            for (iter=chosentaxsetntaxmap.begin();iter!=chosentaxsetntaxmap.end();++iter) {
                nxsstring currenttaxset=iter->first;
                int currentntax=iter->second;
                gsl_matrix *VCV=gsl_matrix_calloc(currentntax,currentntax);
                //gsl_matrix VCV(currentntax,currentntax,0);
                VCV=DeleteStem(GetVCV(currenttaxset));
                //cout<<"taxset "<<currenttaxset<<"  VCV"<<endl<<VCV<<endl;
                // matrixsingular=TestSingularity(VCV);
                //if (matrixsingular) {
                //   goodtree=false;
                //    noerror=false;
                //}
                //else {
                //CholVCVvectorvector=ConvertVCVMatrixToVector(Chol(VCV),CholVCVvectorvector);
                //VCVvectorvector=ConvertVCVMatrixToVector(VCV,VCVvectorvector);
                MakeCombinedVCV(VCVcomb,VCV,ntaxprocessedtreeloop);
                ntaxprocessedtreeloop+=currentntax;


                //if (iter==chosentaxsetntaxmap.begin()) { // It's our first loop
                //    VCVcomb=VCV;
                //}
                //else {
                //    VCVcomb=MakeCombinedVCV(VCVcomb,VCV);
                //}
                //}
				gsl_matrix_free(VCV);
            }
if (goodtree) {
	
    //cout<<"VCV Comb="<<endl<<VCVcomb<<endl<<endl<<"VCVvectorvector (as matrix)"<<endl<<ConvertVCVVectorToMatrix(VCVvectorvector)<<endl<<endl;
    //cout<<"VCVvectorvector"<<endl;
    //cout<<"Length of VCVvectorvector "<<VCVvectorvector->size<<endl<<endl;
    //for (int m=0;m<VCVvectorvector->size;m++) {
    //    cout<<VCVvectorvector[m]<<endl;
    //}
	
for (chosenchar=startchar;chosenchar<=stopchar;chosenchar++) {
    int ntaxprocessedcharloop=0;
    gsl_vector *tipscombresid=gsl_vector_calloc(ntaxcomb);
    //tipscombresid=gsl_vector_calloc(0);
    //gsl_vector tipscombresid(0,0);
    message="Tree = ";
    message+=chosentree;
    message+=": ";
    message+=treename;
    message+=", character = ";
    message+=chosenchar;
    message+="\n";
    if (tablef_open) {
        tmessage="";
        tmessage+=chosentree;
        tmessage+="\t";
        tmessage+=treeweight;
        tmessage+="\t";
        tmessage+=treename;
        tmessage+="\t";
        tmessage+=chosenchar;
        tmessage+="\t";
        tablef<<tmessage;
    }
    int ntaxprocessed=0;
    //Do the test, spit out the values
    map<nxsstring, int>::const_iterator iter2;
    double likelihoodmultiparametermodel=0;
    int taxsetnumber=0;
    for (iter2=chosentaxsetntaxmap.begin();iter2!=chosentaxsetntaxmap.end();++iter2) {
        nxsstring currenttaxset=iter2->first;
        int currentntax=iter2->second;
        taxsetnumber++;
        message+="  Taxset = ";
        message+=currenttaxset;
        message+="\n";
        // gsl_vector *currentVCVvect;
        //currentVCVvect=gsl_vector_calloc(currentntax*currentntax);
        //                        gsl_vector currentVCVvect(currentntax*currentntax,0);
        gsl_matrix *currentVCVmat=gsl_matrix_calloc(currentntax,currentntax);
        currentVCVmat=DeleteStem(GetVCV(currenttaxset));
        //gsl_matrix currentVCVmat(currentntax,currentntax,0);
        gsl_vector *tips=gsl_vector_calloc(currentntax);
        //                        gsl_vector tips(currentntax,0);
        gsl_vector *tipsresid=gsl_vector_calloc(currentntax);
        //gsl_vector tipsresid(currentntax,0);
        double ancstate;
        double rate;
        double likelihood;
        //cout<<"Ntaxprocessed "<<ntaxprocessed<<endl;
        //for (int i=0;i<(currentntax*currentntax);i++) {
        //    currentVCVvect[i]=VCVvectorvector[(ntaxprocessed*ntaxprocessed)+i];
        //    cout<<"currentVCVvect[i] "<<currentVCVvect[i]<<endl;
        // }
        //cout<<"Length of currentVCVvect = "<<currentVCVvect->size<<endl;
        //currentVCVmat=ConvertVCVVectorToMatrix(currentVCVvect);
        //currentVCVmat=ExtractMatrixFromVector(VCVvectorvector,ntaxprocessed,currentntax);
        //matrixsingular=TestSingularity(currentVCVmat);
        //cout<<endl<<endl<<"Current taxset "<<currenttaxset<<endl<<"Current ntax "<<currentntax<<endl<<"CurrentVCVmat:\n"<<currentVCVmat<<endl<<endl;

        //if (matrixsingular) {
        //    cerr<<"Matrix currentVCVmat is singular, ratetest 1"<<endl;
        //    break;
        //}
        ntaxprocessed+=currentntax;
        tips=GetTipValues(currenttaxset,chosenchar);
        //cout<<"Tips of length "<<tips->size<<endl<<"Line 1297"<<endl;
        ancstate=GetAncestralState(currentVCVmat,tips);
        //cout<<endl<<endl<<endl<<"CurrentVCVmat:\n"<<currentVCVmat<<endl;
        //cout<<"Anc state = "<<ancstate<<endl;
        tipsresid=GetTipResiduals(tips,ancstate);
        //for (int debugtaxon=0;debugtaxon<currentntax;debugtaxon++) {
        //   cout<<"Taxon "<<debugtaxon+1<<": "<<tips[debugtaxon]<<"\t"<<tipsresid[debugtaxon]<<endl;
        //}

        rate=EstimateRate(currentVCVmat,tipsresid);
        //cout<<"Rate: "<<rate<<endl;
        if (rate<0) {
            message="Tree ";
            message+=chosentree;
            message+=" excluded: one rate estimate (";
            message+=rate;
            message+=") for taxset ";
			message+=currenttaxset;
			message+=" was negative\n";
            PrintMessage();
            summaryofresults+=chosentree;
            summaryofresults+="\t--Has a negative rate estimate (";
            summaryofresults+=rate;
            summaryofresults+="). No output.--\n";
            if (tablef_open) {
                tmessage="";
                tmessage+=chosentree;
                tmessage+="\t--Has a negative rate estimate. No output.--\n";
                tablef<<tmessage;
            }
            goodtree=false;
            noerror=false;
        }
		else if (rate==0) {
			message="Tree ";
            message+=chosentree;
            message+=" excluded: one rate estimate (";
            message+=rate;
            message+=") for taxset ";
			message+=currenttaxset;
			message+=" was zero (this often happens if there is a single taxon in a clade OR if all taxa have the same value)\n";
            PrintMessage();
            summaryofresults+=chosentree;
            summaryofresults+="\t--Has a zero rate estimate (";
            summaryofresults+=rate;
            summaryofresults+="). No output.--\n";
            if (tablef_open) {
                tmessage="";
                tmessage+=chosentree;
                tmessage+="\t--Has a zero rate estimate. No output.--\n";
                tablef<<tmessage;
            }
            goodtree=false;
            noerror=false;
			
		}
        else {
            gsl_matrix *RateTimesVCVfortest=gsl_matrix_calloc(currentVCVmat->size1,currentVCVmat->size2);
            gsl_matrix_memcpy(RateTimesVCVfortest, currentVCVmat);
            gsl_matrix_scale(RateTimesVCVfortest,rate);
            //RateTimesVCVfortest=rate*currentVCVmat;
            //matrixsingular=TestSingularity(RateTimesVCVfortest);
            //if (matrixsingular==false) {
            likelihood=GetLScore(currentVCVmat,tipsresid,rate);
            likelihoodmultiparametermodel+=likelihood;
            //tipscombresid=MakeCombinedTips(tipscombresid,tipsresid);
            MakeCombinedTips(tipscombresid,tipsresid,ntaxprocessedcharloop);
            ntaxprocessedcharloop+=currentntax;
            message+="    Anc state = ";
            message+=ancstate;
            message+="\n    Rate = ";
            message+=rate;
            message+="\n    -lnL = ";
            message+=likelihood;
            message+="\n";
            if (tablef_open) {
                tmessage="";
                tmessage+=ancstate;
                tmessage+="\t";
                tmessage+=rate;
                tmessage+="\t";
                tmessage+=likelihood;
                tmessage+="\t";
                tablef<<tmessage;
            }
            gsl_vector_set(weightedratevector,taxsetnumber,gsl_vector_get(weightedratevector,taxsetnumber)+(treeweight*rate));
            gsl_vector_set(weightedancstatevector,taxsetnumber,gsl_vector_get(weightedancstatevector,taxsetnumber)+(treeweight*ancstate));
            // }
            // else {
            //     goodtree=false;
            //     noerror=false;
            // }
			gsl_matrix_free(RateTimesVCVfortest);
        }
		gsl_matrix_free(currentVCVmat);
		gsl_vector_free(tips);
		gsl_vector_free(tipsresid);		
    }
    if (goodtree) {
        double ratecomb;
        double likelihoodsingleparametermodel;
        ratecomb=EstimateRate(VCVcomb,tipscombresid);
        likelihoodsingleparametermodel=GetLScore(VCVcomb,tipscombresid,ratecomb); //-lnL actually
        double K1=listedtaxsets+1; //One ancestral state parameter for each taxset, plus one rate parameter
        double K2=2*listedtaxsets; //One ancestral state parameter and one rate parameter for each taxset.
        double model1aicc=(2*likelihoodsingleparametermodel)+2.0*K1+2.0*K1*(K1+1)/(ntaxcomb-K1-1); //AICc, n=1;
        double model2aicc=(2*likelihoodmultiparametermodel)+2.0*K2+2.0*K2*(K2+1)/(ntaxcomb-K2-1);
        double model1aic=(2*likelihoodsingleparametermodel)+2.0*K1;
        double model2aic=(2*likelihoodmultiparametermodel)+2.0*K2;
        if ((ntaxcomb-K1-1)==0) {
            message+="\nWARNING: In the single rate parameter model,\n   there are ";
            message+=K1;
            message+=" parameters to estimate and only\n   ";
            message+=ntaxcomb;
            message+=" datapoints. You need more taxa.\nAICc=Inf [division by zero].\n\n";
        }
        if ((ntaxcomb-K2-1)==0) {
            message+="\n\nWARNING: In the multiple rate parameter model,\n   there are ";
            message+=K2;
            message+=" parameters to estimate but only\n   ";
            message+=ntaxcomb;
            message+=" datapoints. You need more taxa.\nAICc=Inf [division by zero].\n\n";
        }
        tmessage="";
        tmessage+=ratecomb;
        tmessage+="\t";
        tmessage+=K1;
        tmessage+="\t";
        tmessage+=K2;
        tmessage+="\t";
        double aicdif=model1aic-model2aic;
        double absaicdif=fabs(aicdif);
        double aiccdif=model1aicc-model2aicc;
        double absaiccdif=fabs(aiccdif);
        double plrtchi=1-gsl_cdf_chisq_P(2*(likelihoodsingleparametermodel-likelihoodmultiparametermodel), ((2*listedtaxsets)-(listedtaxsets+1))); //TEST THIS
        char outputstring[14];
        // if (debugmode) {
        //     cout<<"One param aic="<<model1aic<<" 2 param aic="<<model2aic<<endl;
        //    cout<<"One param aicc="<<model1aicc<<" 2 param aicc="<<model2aicc<<endl;
        //     cout<<"listedtaxsets="<<listedtaxsets<<" K1="<<K1<<" K2="<<K2<<endl;
        //    cout<<"ntaxcomb="<<ntaxcomb<<endl;
        //cout<<"likelihoodsingleparametermodel="<<likelihoodsingleparametermodel<<endl;
        //cout<<"likelihoodmultiparametermodel="<<likelihoodmultiparametermodel<<endl;
        //    cout<<"(ntaxcomb-K1-1)="<<(ntaxcomb-K1-1)<<endl;
        //    cout<<"(ntaxcomb-K2-1)="<<(ntaxcomb-K2-1)<<endl;

        // }
        message+="\n  Single rate parameter model (model A):\n     -lnL = ";
        sprintf(outputstring,"%14.6f",likelihoodsingleparametermodel);
        message+=outputstring;
        message+="\n    AIC =  ";
        sprintf(outputstring,"%14.6f",model1aic);
        message+=outputstring;
        message+="\n    AICc = ";
        sprintf(outputstring,"%14.6f",model1aicc);
        message+=outputstring;
        message+="\n    rate = ";
        message+=ratecomb;
        gsl_vector_set(weightedratevector,0,gsl_vector_get(weightedratevector,0)+(treeweight*ratecomb));
        //weightedratevector[0]+=ratecomb*treeweight;
        message+="\n\n  Multiple rate parameter model (model B):\n     -lnL = ";
        sprintf(outputstring,"%14.6f",likelihoodmultiparametermodel);
        message+=outputstring;
        message+="\n    AIC =  ";
        sprintf(outputstring,"%14.6f",model2aic);
        message+=outputstring;
        message+="\n    AICc = ";
        sprintf(outputstring,"%14.6f",model2aicc);
        message+=outputstring;
        message+="\n\n";
        weightedAIC1+=treeweight*model1aic;
        weightedAIC2+=treeweight*model2aic;
        weightedAICc1+=treeweight*model1aicc;
        weightedAICc2+=treeweight*model2aicc;


        if (tablef_open) {
            sprintf(outputstring,"%14.6f",model1aic);
            tmessage+=outputstring;
            tmessage+="\t";
            sprintf(outputstring,"%14.6f",model2aic);
            tmessage+=outputstring;
            tmessage+="\t";
            sprintf(outputstring,"%14.6f",model1aicc);
            tmessage+=outputstring;
            tmessage+="\t";
            sprintf(outputstring,"%14.6f",model2aicc);
            tmessage+=outputstring;
            tmessage+="\t";
            sprintf(outputstring,"%14.6f",likelihoodsingleparametermodel);
            tmessage+=outputstring;
            tmessage+="\t";
            sprintf(outputstring,"%14.6f",likelihoodmultiparametermodel);
            tmessage+=outputstring;
            tmessage+="\t";
            tmessage+=aicdif;
            tmessage+="\t";
            tmessage+=aiccdif;
            tmessage+="\t";
            sprintf(outputstring,"%3g",plrtchi);
            tmessage+=outputstring;
            tmessage+="\t";
            tablef<<tmessage;
        }
        double parametricbootstrappingpvalue=0;
        if (repsnumber==0) {
            tmessage="NA\t";
        }
        else { //Do parametric simulation
               //int *seedptr=&seed;
			ProgressBar(repsnumber);
            for(int i=0; i<repsnumber; i++) {
				ProgressBar(0);
                int ntaxprocessed=0;
                gsl_vector *simtipscombined=gsl_vector_calloc(ntaxcomb);
                //simtipscombined=gsl_vector_calloc(0);
                //gsl_vector simtipscombined(0,0);
                double simlikelihoodmultiparametermodel=0;
                gsl_vector *simtipresidcomb=gsl_vector_calloc(ntaxcomb);
                //simtipresidcomb=gsl_vector_calloc(0);
                //gsl_vector simtipresidcomb(0,0);
                map<nxsstring, int>::const_iterator iter3;
                for (iter3=chosentaxsetntaxmap.begin();iter3!=chosentaxsetntaxmap.end();++iter3) {
                    nxsstring taxset=iter3->first;
                    int ntax=iter3->second;
                    double *simtipoutput;
                    //double cholvect[ntax*ntax];
                    //gsl_vector *cholvect;
                    // cholvect=gsl_vector_calloc(0);
                    // gsl_vector cholvect(0,0);
                    //gsl_vector currentVCVvect(ntax*ntax,0);
                    gsl_matrix *currentVCVmat=gsl_matrix_calloc(ntax,ntax);
                    currentVCVmat=DeleteStem(GetVCV(taxset)); //NOTE that this will take a lot of time (generate the VCV for each taxset for each rep). Better to store the VCVs once.
                                                              //currentVCVmat=ExtractMatrixFromVector(VCVvectorvector,ntaxprocessed,ntax);
                                                              //cholvect=ConvertVCVMatrixToVector(ExtractMatrixFromVector(CholVCVvectorvector,ntaxprocessed,ntax),cholvect);
                                                              //double cholvectdouble[ntax*ntax];
                                                              //for (int i=0;i<(ntax*ntax);i++) {
                                                              //  cholvectdouble[i]=cholvect[i];
                                                              //      currentVCVvect[i]=VCVvectorvector[(ntaxprocessed*ntaxprocessed)+i];
                                                              // }
                                                              //  currentVCVmat=ConvertVCVVectorToMatrix(currentVCVvect);

                    gsl_vector *nullmean=gsl_vector_calloc(ntax);
                    gsl_vector *simtipvector=gsl_vector_calloc(ntax);
                    simtipvector=SimulateTips(currentVCVmat, ratecomb, nullmean);
                    gsl_vector *simtipsresid=gsl_vector_calloc(ntax);
                    //gsl_vector simtipsresid(ntax,0);
                    double simancstate;
                    double simrate;
                    double simlikelihood;
                    simancstate=GetAncestralState(currentVCVmat,simtipvector);
                    simtipsresid=GetTipResiduals(simtipvector,simancstate);
                    simrate=EstimateRate(currentVCVmat,simtipsresid);
					//cout<<simrate<<"\t";
                    if (simrate==0) {
                        cerr<<"Warning: Rate in one parametric simulation was zero.\nChanging to a very tiny number";
                        simrate=1.0e-10;
                    }
                    simlikelihood=GetLScore(currentVCVmat,simtipsresid,simrate);
                    simlikelihoodmultiparametermodel+=simlikelihood;
                    MakeCombinedTips(simtipresidcomb,simtipsresid,ntaxprocessed);
                    ntaxprocessed+=ntax;
					gsl_matrix_free(currentVCVmat);
					gsl_vector_free(nullmean);
					gsl_vector_free(simtipvector);
					gsl_vector_free(simtipsresid);					
                }
                double simratecomb;
                double simlikelihoodsingleparametermodel;
                simratecomb=EstimateRate(VCVcomb,simtipresidcomb);
                simlikelihoodsingleparametermodel=GetLScore(VCVcomb,simtipresidcomb,simratecomb);
				//cout<<"-\tsim like (sim, then mult) "<<simlikelihoodsingleparametermodel<<"\t"<<simlikelihoodmultiparametermodel<<"\t-\treal like (single, multi) "<<likelihoodsingleparametermodel<<"\t"<<likelihoodmultiparametermodel<<"\t-\tsim rate: "<<simratecomb<<"\test rate: "<<ratecomb<<endl;
				//cout<<endl;
				//for (int i=0;i<simtipresidcomb->size;i++) {
				//	cout<<gsl_vector_get(simtipresidcomb,i)<<" ";
				//}
				//cout<<endl;
				//cout<<"\trates\t"<<simratecomb<<"\t"<<ratecomb<<endl;
                if ((simlikelihoodsingleparametermodel-simlikelihoodmultiparametermodel)>=(likelihoodsingleparametermodel-likelihoodmultiparametermodel)) {
                    parametricbootstrappingpvalue+=(1.0/repsnumber); //Need 1.0 rather than 1 or C++ assumes you want an int output
                }
				gsl_vector_free(simtipscombined);
				gsl_vector_free(simtipresidcomb);
            }
            tmessage="";
			sprintf(outputstring,"%3g",parametricbootstrappingpvalue);
            tmessage+=outputstring;
            tmessage+="\t";
            //PrintMessage();
            if (tablef_open) {
                tablef<<tmessage;
            }
        }

        /////////
        //print selected model under various methods
        summaryofresults+=chosentree;
        summaryofresults+="\t";
        summaryofresults+=chosenchar;
        summaryofresults+="\t";

        if (model1aic<model2aic) {
            if (absaicdif<10) {
                message+="a: AIC dif of ";
                message+=absaicdif;
                message+=" somewhat favors the single rate parameter model.\n";
                summaryofresults+="a\t";
            }
            else {
                message+="A: AIC dif of ";
                message+=absaicdif;
                message+=" strongly favors the single rate parameter model.\n";
                summaryofresults+="A\t";

            }
        }
        else {
            if (absaicdif<10) {
                message+="b: AIC dif of ";
                message+=absaicdif;
                message+=" somewhat favors the multiple rate parameter model.\n";
                summaryofresults+="b\t";
            }
            else {
                message+="B: AIC dif of ";
                message+=absaicdif;
                message+=" strongly favors the multiple rate parameter model.\n";
                summaryofresults+="B\t";
            }
        }

        if (model1aicc<model2aicc) {
            if (absaiccdif<10) {
                message+="a: AICc dif of ";
                message+=absaiccdif;
                message+=" somewhat favors the single rate parameter model.\n";
                summaryofresults+="a\t";
            }
            else {
                message+="A: AICc dif of ";
                message+=absaiccdif;
                message+=" strongly favors the single rate parameter model.\n";
                summaryofresults+="A\t";
            }
        }
        else {
            if (absaiccdif<10) {
                message+="b: AICc dif of ";
                message+=absaiccdif;
                message+=" somewhat favors the multiple rate parameter model.\n";
                summaryofresults+="b\t";
            }
            else {
                message+="B: AICc dif of ";
                message+=absaiccdif;
                message+=" strongly favors the multiple rate parameter model.\n";
                summaryofresults+="B\t";

            }
        }

        if (plrtchi>0.05) {
            message+="a: Chi-square p of ";
            message+=plrtchi;
            message+=" does not reject the single rate parameter model.\n";
            summaryofresults+="a\t";
        }
        else {
            if (plrtchi<0.01) {
                message+="B: Chi-square p of ";
                message+=plrtchi;
                message+=" rejects the single rate parameter model.\n";
                summaryofresults+="B\t";
            }
            else {
                message+="b: Chi-square p of ";
                message+=plrtchi;
                message+=" weakly rejects the single rate parameter model.\n";
                summaryofresults+="b\t";

            }
        }
        weightedchip+=plrtchi*treeweight;
        if (repsnumber==0) {
            message+="?: Parametric bootstrapping not done, despite bias in chi-square test.\n";
            summaryofresults+="?\t";
        }
        else {
            weightedparamp+=parametricbootstrappingpvalue*treeweight;
            if (parametricbootstrappingpvalue>0.05) {
                message+="a: Parametric bootstrap p of ";
                message+=parametricbootstrappingpvalue;
                message+=" does not reject the single rate parameter model.\n";
                summaryofresults+="a\t";
            }
            else {
                if (parametricbootstrappingpvalue==0) {
                    message+="B: Parametric bootstrap p of <";
                    message+=1.0/repsnumber;
                    message+=" rejects the single rate parameter model.\n";
                    summaryofresults+="B\t";
                }
                else if (parametricbootstrappingpvalue<0.01) {
                    message+="B: Parametric bootstrap p of ";
                    message+=parametricbootstrappingpvalue;
                    message+=" rejects the single rate parameter model.\n";
                    summaryofresults+="B\t";

                }
                else {
                    message+="b: Parametric bootstrap p of ";
                    message+=parametricbootstrappingpvalue;
                    message+=" weakly rejects the single rate parameter model.\n";
                    summaryofresults+="b\t";
                }
            }
        }


        //print selected model under various methods
        if (model1aic<model2aic) {
            if (absaicdif<10) {
                tmessage="a";
            }
            else {
                tmessage="A";
            }
        }
        else {
            if (absaicdif<10) {
                tmessage="b";
            }
            else {
                tmessage="B";
            }
        }

        if (model1aicc<model2aicc) {
            if (absaiccdif<10) {
                tmessage+="a";
            }
            else {
                tmessage+="A";
            }
        }
        else {
            if (absaiccdif<10) {
                tmessage+="b";
            }
            else {
                tmessage+="B";
            }
        }

        if (plrtchi>0.05) {
            tmessage+="a";
        }
        else {
            if (plrtchi<0.01) {
                tmessage+="B";
            }
            else {
                tmessage+="b";
            }
        }

        if (repsnumber==0) {
            tmessage+="?";
        }
        else {
            if (parametricbootstrappingpvalue>0.05) {
                tmessage+="a";
            }
            else {
                if (parametricbootstrappingpvalue<0.01) {
                    tmessage+="B";
                }
                else {
                    tmessage+="b";
                }
            }
        }
        if (notquietmode) {
            PrintMessage();
        }
        summaryofresults+="\n";
        if (tablef_open) {
            tmessage+="\n";
            tablef<<tmessage;
        }

    }
	gsl_vector_free(tipscombresid);	
}
}
else {
    weighttotal-=treeweight;
    message="Tree ";
    message+=chosentree;
    message+=" excluded: has at least one subtree matrix that is singular\n";
    PrintMessage();
    summaryofresults+=chosentree;
    summaryofresults+="\t--Has a singular subtree VCV. No output.--\n";
    if (tablef_open) {
        tmessage="";
        tmessage+=chosentree;
        tmessage+="\t--Has a singular subtree VCV. No output.--\n";
        tablef<<tmessage;
    }
}
gsl_matrix_free(VCVcomb);
        }
chosentree=originalchosentree; //restore initial values.
chosenchar=originalchosenchar;
message=summaryofresults;
PrintMessage();
if (charloop==false && treeloop==true) {
    message="\n\nAverage result (weighted by tree weights)\n";
    if (noerror==false) {
        message+="\tError: Some trees had errors\n";
        PrintMessage();
    }
    else {
        map<nxsstring, int>::const_iterator iter5;
        int taxsetnumber=0;
        for (iter5=chosentaxsetntaxmap.begin();iter5!=chosentaxsetntaxmap.end();++iter5) {
            nxsstring taxset=iter5->first;
            taxsetnumber++;
            message+="\tTaxset ";
            message+=taxset;
            message+=": rate = ";
            message+=gsl_vector_get(weightedratevector,taxsetnumber)/weighttotal;
            message+=", anc state = ";
            message+=gsl_vector_get(weightedancstatevector,taxsetnumber)/weighttotal;
            message+="\n";
        }
        message+="\tSingle rate parameter model rate = ";
        message+=gsl_vector_get(weightedratevector,0)/weighttotal;
        message+="\n\n        Single rate parameter           Multiple rate parameters\nAIC =     ";
        char outputstring[14];
        sprintf(outputstring,"%14.6f",weightedAIC1/weighttotal);
        message+=outputstring;
        message+="                    ";
        sprintf(outputstring,"%14.6f",weightedAIC2/weighttotal);
        message+=outputstring;
        message+="\nAICc =    ";
        sprintf(outputstring,"%14.6f",weightedAICc1/weighttotal);
        message+=outputstring;
        message+="                    ";
        sprintf(outputstring,"%14.6f",weightedAICc2/weighttotal);
        message+=outputstring;
        message+="\n\nChi-square p value = ";
        sprintf(outputstring,"%14.6f",weightedchip/weighttotal);
        message+=outputstring;
        message+="\nParametric bootstrap p value = ";
        if (repsnumber>0) {
            sprintf(outputstring,"%14.6f",weightedparamp/weighttotal);
            message+=outputstring;
        }
        else {
            message+="NA: Not done";
        }
        PrintMessage();
    }
}
gsl_vector_free(weightedratevector);
gsl_vector_free(weightedancstatevector);
        }
if (tablef_open) {
    tablef.close();
}
    }



/**@method NumOpt
*Allows setting of parameters for numerical optimization
*/
void BROWNIE::NumOpt( NexusToken& token)
{
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else if( token.Abbreviation("?") ) {
            message="Usage: NumOpt [options...]\n\n";
            message+="This sets options for numerical optimization\n";
            message+="Available options:\n\n";
            message+="Keyword ------ Option type ------------------------ Current setting --";
            message+="\n Iter          <integer>                            ";
            message+=maxiterations;
            message+="\n Toler         <double>                             ";
            char outputstring[10];
            sprintf(outputstring,"%1.9f",stoppingprecision);
            message+=outputstring;
            message+="\n RandStart     <integer>                            ";
            message+=randomstarts;
            message+="\n Seed          <integer>                            ";
            message+=gslseedtoprint;
            message+="\n StepSize      <double>                             ";
            sprintf(outputstring,"%1.9f",stepsize);
            message+=outputstring;
            message+="\n Detail        Yes|No                               ";
            if (detailedoutput) {
                message+="Yes";
            }
            else {
                message+="No";
            }
			message+="\n Redo          Yes|No                               ";
            if (redobad) {
                message+="Yes";
            }
            else {
                message+="No";
            }
			message+="\n GiveUpFactor  <integer>                            ";
            message+=giveupfactor;

            message+="\n\nIter sets the maximum number of iterations of the Nelder-Mead simplex algorithm.\n\nToler sets the precision of the stopping criterion: what amount\nof change in the likelihood is considered small enough to count as zero change.\n\nRandStart sets the number of random starts to use.\n\nStepSize sets the NM step size.\n\nDetail specifies whether or not to have detailed output from numerical optimization\n\nRedo specifies whether to redo reps which stop due to iteration limits\n\nGiveUpFactor, when redo=yes, is used to tell the software when to stop restarting: when the ratio of unsuccessful to successful starts is > giveupfactor";
            PrintMessage();
        }
        else if( token.Abbreviation("Iter") ) {
            nxsstring numbernexus = GetNumber(token);
            maxiterations=atoi( numbernexus.c_str() ); //convert to int
        }
        else if( token.Abbreviation("TOler") ) {
            nxsstring numbernexus = GetNumber(token);
            stoppingprecision=atof( numbernexus.c_str() );
        }
        else if( token.Abbreviation("RAndstart") ) {
            nxsstring numbernexus = GetNumber(token);
            randomstarts=atoi( numbernexus.c_str() ); //convert to int
        }
        else if( token.Abbreviation("SEed") ) {
            nxsstring numbernexus = GetNumber(token);
            gslseedtoprint=atoi( numbernexus.c_str() );
            gslseed=gslseedtoprint;
            gsl_rng_set(r,gslseed); //convert to int
            message="Sorry, setting the seed doesn't quite work yet\n";
            PrintMessage();
        }
        else if( token.Abbreviation("Detail") ) {
            nxsstring yesnodetail=GetFileName(token);
            if (yesnodetail[0] == 'n') {
                detailedoutput=false;
            }
            else {
                detailedoutput=true;
            }
        }
		else if( token.Abbreviation("REdo") ) {
            nxsstring yesnodetail=GetFileName(token);
            if (yesnodetail[0] == 'n') {
                redobad=false;
            }
            else {
                redobad=true;
            }
        }
		else if( token.Abbreviation("Giveupfactor") ) {
            nxsstring numbernexus = GetNumber(token);
            giveupfactor=atoi( numbernexus.c_str() ); //convert to int
        }		
        else if( token.Abbreviation("STepsize") ) {
            nxsstring numbernexus = GetNumber(token);
            stepsize=atof( numbernexus.c_str() );
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading NUMOPT command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }

}


void BROWNIE::HandleOrderByTree( NexusToken& token)
{
	bool donesomething=false;
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else if( token.Abbreviation("?") ) {
            message="This will reorder taxa and data so that taxa adjacent on the tree\nare adjacent in the data matrix\n\n";
            PrintMessage();
			donesomething=true;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading GARLAND command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
	if (!donesomething) { //Means we want to output the matrix
		if (intrees.GetNumTrees()==0) {
			errormsg="First load a tree (include in the input file)";
			throw XNexus( errormsg );
		}
		else {
			nxsstring exportfname="ReorderedTaxonNames.nex";
			exportf.open( exportfname.c_str(), ios::out | ios::app );
			Tree t = intrees.GetIthTree(chosentree-1);
			NodeIterator <Node> n (t.GetRoot());
			NodePtr cur = n.begin();
			while (cur) {
				if (cur->IsLeaf()) {
					exportf<<cur->GetLabel()<<" ";
					nxsstring currenttaxonlabel=(cur->GetLabel()).c_str();
					int taxonnumber=TaxonLabelToNumber( currenttaxonlabel );
					for (int charnumber=0; charnumber<(discretecharacters->GetNChar()); charnumber++) {
						discretecharacters->ShowStateLabels(exportf, taxonnumber-1, charnumber);
					}
					exportf<<endl;
				}
				cur = n.next();
			}
			exportf.close();
		}
	}
}

/**@method Garland
*Spits out values to use to compare Brownie's performance with other programs
*/
void BROWNIE::HandleGarland( NexusToken& token)
{
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else if( token.Abbreviation("?") ) {
            message="This will perform various tests. It assumes you have set a taxset named ALL\n\n";
            PrintMessage();
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading GARLAND command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
    //Do something
    int ntax=0;
    nxsstring chosentaxset="ALL";
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }
IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntax++;
    }
    gsl_matrix* VCV=gsl_matrix_calloc(ntax,ntax);
    VCV=DeleteStem(GetVCV(chosentaxset));
    gsl_permutation * p = gsl_permutation_alloc (ntax);
    int signum;
    gsl_matrix * VCVLU=gsl_matrix_calloc(ntax,ntax);
    gsl_matrix_memcpy (VCVLU, VCV);
    gsl_linalg_LU_decomp (VCVLU,p, &signum);
    int sig=signum;
    double determinant =  gsl_linalg_LU_det (VCVLU,sig);
    message="Determinant is ";
    message+=determinant;
    PrintMessage();
    gsl_matrix * VCVrescaled=gsl_matrix_calloc(ntax,ntax);
    gsl_matrix_memcpy (VCVrescaled, VCV);
    gsl_matrix_scale (VCVrescaled,1/determinant);
    ofstream vcvf;
    vcvf.open( "vcv.txt" );
    for (int currentrow=0;currentrow<ntax;currentrow++) {
        for (int currentcol=0;currentcol<ntax;currentcol++) {
            vcvf<<gsl_matrix_get(VCVrescaled,currentrow,currentcol);
            vcvf<<"\t";
        }
        vcvf<<"\n";
    }
    vcvf.close();
    message+="The rescaled VCV matrix has been saved to vcv.txt";
    PrintMessage();
    gsl_vector * tips=gsl_vector_calloc(ntax);
    tips=GetTipValues(chosentaxset,chosenchar);
    double meantips=0;
    for (int taxonnum=0;taxonnum<ntax;taxonnum++) {
        meantips+=gsl_vector_get(tips,taxonnum);
    }
    meantips=meantips/ntax;
    double ancestralstate=GetAncestralState(VCVrescaled,tips);
    message="Mean tip value is ";
    message+=meantips;
    message+="\nReconstructed ancestral state is ";
    message+=ancestralstate;
    PrintMessage();
    gsl_vector *variance =gsl_vector_calloc(ntax);
    int originalmodel=chosenmodel;
    chosenmodel=1;




    OptimizationFn my_fnA(VCVrescaled,tips,variance,maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput);


    gsl_vector *optimalrateA=gsl_vector_calloc(2);
				//gsl_vector_memcpy(optimalrate,my_fn.OptimizeRateWithGivenTipVariance());
    gsl_vector_memcpy(optimalrateA,my_fnA.GeneralOptimization(1));
    message="\nOptimal rate = ";
    message+=gsl_vector_get(optimalrateA,0);
    message+=" +/- ";
    message+=gsl_vector_get(optimalrateA,1);
    PrintMessage();
    OptimizationFn my_fnB(VCVrescaled,tips,variance,maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput);
    chosenmodel=3;
    message="OU model chosen";
    PrintMessage();
    gsl_vector *optimalrateB=gsl_vector_calloc(6);
    gsl_vector_memcpy(optimalrateB,my_fnB.GeneralOptimization(3));
    message="\nOptimal rate = ";
    message+=gsl_vector_get(optimalrateB,0);
    message+=" +/- ";
    message+=gsl_vector_get(optimalrateB,3);
    message+="\nAncestral state = ";
    message+=gsl_vector_get(optimalrateB,1);
    message+=" +/- ";
    message+=gsl_vector_get(optimalrateB,4);
    message+="\nd = ";
    message+=gsl_vector_get(optimalrateB,2);
    message+=" +/- ";
    message+=gsl_vector_get(optimalrateB,5);
    chosenmodel=4;
    PrintMessage();
    message="ACDC model chosen";
    PrintMessage();
    OptimizationFn my_fnC(VCVrescaled,tips,variance,maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput);

    gsl_vector *optimalrateC=gsl_vector_calloc(6);
    gsl_vector_memcpy(optimalrateC,my_fnC.GeneralOptimization(4));
    message="\nOptimal rate = ";
    message+=gsl_vector_get(optimalrateC,0);
    message+=" +/- ";
    message+=gsl_vector_get(optimalrateC,3);
    message+="\nAncestral state = ";
    message+=gsl_vector_get(optimalrateC,1);
    message+=" +/- ";
    message+=gsl_vector_get(optimalrateC,4);
    message+="\ng = ";
    message+=gsl_vector_get(optimalrateC,2);
    message+=" +/- ";
    message+=gsl_vector_get(optimalrateC,5);
    PrintMessage();



    chosenmodel=originalmodel;
	gsl_matrix_free(VCV);
	gsl_matrix_free(VCVLU);
	gsl_matrix_free(VCVrescaled);
	gsl_vector_free(tips);
	gsl_vector_free(variance);
	gsl_vector_free(optimalrateA);
	gsl_vector_free(optimalrateB);
	gsl_vector_free(optimalrateC);
	gsl_permutation_free (p);
}


/**
* @method HandleTipVariance [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the TipVariance command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleTipVariance( NexusToken& token )
{
    bool starting = false;
    bool stopping = false;
    bool appending = false;
    bool replacing = false;
    bool name_provided = false;
    bool help_provided = false;

    // Retrieve all tokens for this command, stopping only in the event
    // of a semicolon or an unrecognized keyword
    //
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else if( token.Abbreviation("?") ) {
            help_provided = true;
            message="Usage: TipVariance [options...]\n\n";
            message+="Observations of taxon means have uncertainty. By default, the program assumes that these values are known exactly.\nIf you have measured these variances for each taxon (so that, for each measured character,\nthe data matrix has the mean followed by the variance), choose Given.\nIf you want estimate one tip variance across all the taxa, choose Same.\nIf you want to assume no tip variance (one fewer parameter), choose None\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --";
            message+="\n Type        Given|Same|None                      ";
            if (tipvariancetype==0) {
                message+="None";
            }
            else if(tipvariancetype==1) {
                message+="Same";
            }
            else if(tipvariancetype==2) {
                message+="Given";
            }
            else {
                message+=tipvariancetype;
            }
            message+="\n\nNote that if the tip variance type is set to given, values in the matrix are\nassumed to consist of a taxon average followed by the\ncorresponding variance, so the matrix might be\n";
            message+=" taxon1   char1_mean  char1_variance  char2_mean  char2_variance  char3_mean  char3_variance...\n";
            message+=" taxon2   char1_mean  char1_variance  char2_mean  char2_variance  char3_mean  char3_variance...\n";
            message+=" taxon3   char1_mean  char1_variance  char2_mean  char2_variance  char3_mean  char3_variance...\n";
            message+="\nThis differs somewhat from the NEXUS specification (Maddison et al, 1997), which would use an\nITEMS=(AVERAGE VARIANCE) command and a somewhat different matrix.\n";
            PrintMessage();
        }
        else if( token.Abbreviation("Type") ) {
            nxsstring tipvarianceinput=GetFileName(token);
            if (tipvarianceinput[0] == 'N' || tipvarianceinput[0] == 'n') {
                tipvariancetype=0;
            }
            else if(tipvarianceinput[0] == 'S' || tipvarianceinput[0] == 's') {
                tipvariancetype=1;
            }
            else if(tipvarianceinput[0] == 'G' || tipvarianceinput[0] == 'g') {
                tipvariancetype=2;
            }
            else {
                errormsg = "Unexpected option (";
                errormsg += tipvarianceinput;
                errormsg += ") encountered reading TipVariance command";
                throw XNexus( errormsg);
            }
            message="\n You have chosen ";
            if (tipvariancetype==0) {
                message+="to assume that the tip variance is zero.";
            }
            else if(tipvariancetype==1) {
                message+="to estimate one tip variance parameter.";
            }
            else if(tipvariancetype==2) {
                message+="to use given variances from the data matrix.\nThis will assume that odd characters (char 1, 3, 5...) are taxon averages\nand even characters are the corresponding variances.\n";
            }
            else {
                message+=tipvariancetype;
            }
            PrintMessage();
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading TipVariance command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}





/**
* @method HandleEcho [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the ECHO command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleEcho( NexusToken& token )
{
    bool starting = false;
    bool stopping = false;
    bool appending = false;
    bool replacing = false;
    bool name_provided = false;
    bool help_provided = false;

    // Retrieve all tokens for this command, stopping only in the event
    // of a semicolon or an unrecognized keyword
    //
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else if( token.Abbreviation("?") ) {
            help_provided = true;
            message="Usage: Echo [options...]\n\n";
            message+="Records everything you type into a Brownie-executable batch file,\nsuitable for pasting at the end of your data file\n(delete the #nexus if you do this)\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --";
            message+="\nFile         <echo-file-name>                     ";
            if (echof_open) {
                message+=echofname;
            }
            else {
                message+="[no file open]";
            }
            message+="\nReplace      No|Yes                              *No";
            message+="\nAppend       No|Yes                              *No";
            message+="\n          Start|Stop";
            message+="\n                                                 *Option is nonpersistent\n\nNote that you should specify =yes or =no for Replace \nor Append, but not for Start or Stop.\nFor example,\n\echo start file=test.nex replace=yes;\n\n";
            message+="\n                                                 *Option is nonpersistent\n\n";
            PrintMessage();
        }
        else if( token.Abbreviation("STOp") ) {
            stopping=true;
        }
        else if( token.Abbreviation("STArt") ) {
            starting=true;
        }
        else if( token.Abbreviation("Replace") ) {
            nxsstring yesnoreplace=GetFileName(token);
            if (yesnoreplace[0] == 'n') {
                replacing=false;
            }
            else {
                replacing=true;
            }
        }
        else if( token.Abbreviation("Append") ) {
            nxsstring yesnoappend=GetFileName(token);
            if (yesnoappend[0] == 'n') {
                appending=false;
            }
            else {
                appending=true;
            }
        }
        else if( token.Abbreviation("File") ) {
            echofname = GetFileName(token);
            name_provided = true;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading ECHO command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }

    // Check for incompatible combinations of keywords
    //
    if( stopping && ( starting || appending || replacing || name_provided ) ) {
        errormsg = "Cannot specify STOP with any of the following START, APPEND, REPLACE, FILE";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    if( appending && replacing ) {
        errormsg = "Cannot specify APPEND and REPLACE at the same time";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    if( echof_open && ( starting || name_provided || appending || replacing ) ) {
        errormsg = "Cannot start echo file since echo file is already open";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    // Is user closing an open echo file?
    //
    if( stopping ) {
        echof<<"end;"<<endl;
        echof.close();
        echof_open = false;

        message = "\nEcho file closed";
        PrintMessage();

        return;
    }

    // If this far, must be attempting to open a echo file
    //
    if( !name_provided && !help_provided) {
        errormsg = "Must provide a file name when opening a echo file\n";
        errormsg += "e.g., echo file=batch.nex start replace;";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    if( appending ) {
        echof_open = true;
        echof.open( echofname.c_str(), ios::out | ios::app );

        message = "\nAppending to echo file ";
        message += echofname;
        PrintMessage();
        echof<<"Begin Brownie;"<<endl;
    }
    else if( replacing ) {
        echof_open = true;
        echof.open( echofname.c_str() );

        message = "\nReplacing echo file ";
        message += echofname;
        PrintMessage();
        echof<<"#nexus"<<endl<<"Begin Brownie;"<<endl;
    }
    else  {
        bool exists = FileExists( echofname.c_str() );
        bool userok = true;
        if( exists && !UserSaysOk( "Ok to replace?", "Echo file specified already exists" ) )
            userok = false;
        if( userok ) {
            echof_open = true;
            echof.open( echofname.c_str() );
            echof<<"#nexus"<<endl<<"Begin Brownie;"<<endl;
        }

        if( exists && userok ) {
            message = "\nReplacing echo file ";
            message += echofname;
        }
        else if( userok ) {
            message = "\nEcho file ";
            message += echofname;
            message += " opened";
        }
        else {
            message = "\nEcho command aborted";
        }
        PrintMessage();
    }
}




/**
* @method HandleLog [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * Called when the LOG command needs to be parsed
 * from within the BROWNIE block.
 */
void BROWNIE::HandleLog( NexusToken& token )
{
    bool starting = false;
    bool stopping = false;
    bool appending = false;
    bool replacing = false;
    bool name_provided = false;
    bool help_provided = false;


    // Retrieve all tokens for this command, stopping only in the event
    // of a semicolon or an unrecognized keyword
    //
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            break;
        }
        else if( token.Abbreviation("?") ) {
            help_provided=true;
            message="Usage: Log [options...]\n\n";
            message+="Records output into a log file.\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --";
            message+="\nFile         <log-file-name>                      ";
            if (logf_open) {
                message+=logfname;
            }
            else {
                message+="[no file open]";
            }
            message+="\nReplace      No|Yes                              *No";
            message+="\nAppend       No|Yes                              *No";
            message+="\n          Start|Stop";
            message+="\n                                                 *Option is nonpersistent\n\nNote that you should specify =yes or =no for Replace \nor Append, but not for Start or Stop.\nFor example,\n\tlog start file=test.log replace=yes;\n\n";
            PrintMessage();
        }
        else if( token.Abbreviation("STOp") ) {
            stopping=true;
        }
        else if( token.Abbreviation("STArt") ) {
            starting=true;
        }
        else if( token.Abbreviation("Replace") ) {
            nxsstring yesnoreplace=GetFileName(token);
            if (yesnoreplace[0] == 'n') {
                replacing=false;
            }
            else {
                replacing=true;
            }
        }
        else if( token.Abbreviation("Append") ) {
            nxsstring yesnoappend=GetFileName(token);
            if (yesnoappend[0] == 'n') {
                appending=false;
            }
            else {
                appending=true;
            }
        }
        else if( token.Abbreviation("File") ) {
            logfname = GetFileName(token);
            name_provided = true;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading LOG command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }

    // Check for incompatible combinations of keywords
    //
    if( stopping && ( starting || appending || replacing || name_provided ) ) {
        errormsg = "Cannot specify STOP with any of the following START, APPEND, REPLACE, FILE";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    if( appending && replacing ) {
        errormsg = "Cannot specify APPEND and REPLACE at the same time";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    if( logf_open && ( starting || name_provided || appending || replacing ) ) {
        errormsg = "Cannot start log file since log file is already open";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    // Is user closing an open log file?
    //
    if( stopping ) {
        logf.close();
        logf_open = false;

        message = "\nLog file closed";
        PrintMessage();

        return;
    }

    // If this far, must be attempting to open a log file
    //
    if( !name_provided && !help_provided) {
        errormsg = "Must provide a file name when opening a log file\n";
        errormsg += "e.g., log file=doofus.txt start replace;";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    if( appending ) {
        logf_open = true;
        logf.open( logfname.c_str(), ios::out | ios::app );

        message = "\nAppending to log file ";
        message += logfname;
        PrintMessage();
    }
    else if( replacing ) {
        logf_open = true;
        logf.open( logfname.c_str() );

        message = "\nReplacing log file ";
        message += logfname;
        PrintMessage();
    }
    else  {
        bool exists = FileExists( logfname.c_str() );
        bool userok = true;
        if( exists && !UserSaysOk( "Ok to replace?", "Log file specified already exists" ) )
            userok = false;
        if( userok ) {
            logf_open = true;
            logf.open( logfname.c_str() );
        }

        if( exists && userok ) {
            message = "\nReplacing log file ";
            message += logfname;
        }
        else if( userok ) {
            message = "\nLog file ";
            message += logfname;
            message += " opened";
        }
        else {
            message = "\nLog command aborted";
        }
        PrintMessage();
    }
}

/**
* @method HandleNextCommand [void:public]
 *
 * Accepts a string in the form of a BROWNIE block containing one command
 * and processes it just like a real BROWNIE block in a NEXUS data file.
 */
void BROWNIE::HandleNextCommand()
{
std::istrstream cmdin( next_command );
    //std::istringstream cmdin( next_command );
    //cerr << "next_command is: " << next_command <<endl;
    //istream cmdin( next_command );

    NexusToken token(cmdin);
    try {
        Read( token );
    }
    catch( XNexus x )
    {
        NexusError( errormsg, x.pos, x.line, x.col );
        Reset();
    }
}

/**
* @method NexusError [virtual void:public]
 * @param msg [nxsstring&] the error message
 * @param pos [streampos] the point in the NEXUS file where the error occurred
 * @param line [long] the line in the NEXUS file where the error occurred
 * @param col [long] the column in the NEXUS file where the error occurred
 *
 * Called when an error is encountered in a NEXUS file. Allows program to
 * give user details of the error as well as the precise location of the
 * error. Virtual function that overrides the pure virtual function in the
 * base class Nexus.
 */
void BROWNIE::NexusError( nxsstring& msg, streampos /* pos */, long line, long col )
{
    message = "\n";
    message += msg;
    PrintMessage();

    if( inf_open )
    {
        message = "Line:   ";
        message += line;
        PrintMessage();

        message = "Column: ";
        message += col;
        PrintMessage();
    }
	if (quit_onerr) {
		message = "Quitting on error. To change this, type \"noquitonerr\"";
		PrintMessage();
		quit_now=true;
		exit(EXIT_FAILURE);
	}	
}

/**
* @method PreprocessNextCommand [void:public]
 *
 * Begins with the command just entered by the user, which is stored in
 * the data member next_command, adds a semicolon (if the user failed
                                                   * to supply one), and then adds "end;" so the whole bundle looks
 * like a very short BROWNIE block.  This is then passed to HandleNextCommand,
 * which processes it just like a real BROWNIE block in a NEXUS data file.
 */
void BROWNIE::PreprocessNextCommand()
{
    // If user failed to add the terminating semicolon,
    // we'll do it now. We will also remove the line feed
    // at the end and add the command "end;" to the end
    // of the line (see explanation below).
    //
    int len = strlen(next_command);
    //assert( len > 0 );

    // Remove any whitespace characters from end of string entered by user
    //
    int i = len;
    while( i > 0 && next_command[i-1] == ' ' || next_command[i-1] == '\t' || next_command[i-1] == '\n' )
        i--;


    // If character at position i-1 is a semicolon, put '\0' terminator at position i;
    // otherwise, put a semicolon at position i and terminator at i+1
    //
    if( next_command[i-1] != ';' ) {
        next_command[i] = ';';
        i++;
    }
	if (i> COMMAND_MAXLEN) {
		cout<<"The maximum command length is "<<COMMAND_MAXLEN<<" but the following line is "<<i<<" characters long"<<endl;
		nxsstring toolong="";
		toolong+=next_command;
		cout<<toolong<<endl;
	}
    assert( i <= COMMAND_MAXLEN );
    next_command[i] = '\0';

    if( echof_open ) {
        echof<<next_command<<endl;
    }

    // Now add a semicolon at the beginning and terminate with an "END;" command
    // so that we can pretend this is simply a very short private NEXUS block
    // containing only one command.  This allows us to simply use the Read
    // function we inherited from the base class BstBase to process the command.
    //
    len = strlen(next_command);
    assert( len < COMMAND_MAXLEN-2 );
    nxsstring tmp = ";";
    tmp += next_command;
    tmp += "end;";
    strcpy( next_command, tmp.c_str() );
}

/**
* @method PrintMessage [void:public]
 * @param linefeed [bool] if true, places newline character after message
 *
 * All output handled here.  Writes string currently stored in message
 * (a nxsstring data member) to the output file stream, if open, and also
 * to the console via cerr. Places newline after string if linefeed is true.
 */
void BROWNIE::PrintMessage( bool linefeed /* = true */ )
{
    cerr << message;
    if( linefeed )
        cerr << endl;

    if( logf_open ) {
        logf << message;
        if( linefeed )
            logf << endl;
    }
}

/**
* @method PurgeBlocks [void:protected]
 *
 * Detaches all blocks, deletes them, creates new blocks, and
 * finally adds the new blocks.
 */
void BROWNIE::PurgeBlocks()
{
    Detach( taxa );
    Detach( trees );
    Detach( assumptions );
    Detach( characters );

    delete characters;
    delete assumptions;
    delete trees;
    delete taxa;

    taxa = new TaxaBlock();
    trees = new TreesBlock(*taxa);
    assumptions = new AssumptionsBlock( *taxa );
    characters = new CharactersBlock( *taxa, *assumptions );
	characters2 = new CharactersBlock2( *taxa, *assumptions );

    Add( taxa );
    Add( trees );
    Add( assumptions );
    Add( characters );
	Add( characters2 );
	
	
}

void BROWNIE::Assign(NexusToken& token)
{
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            message="Current assignment vector is:\n";
            for (int i=0; i<convertsamplestospecies.size(); i++) {
                message+=convertsamplestospecies[i];
                message+=" ";
            }
            PrintMessage();
            break;
        }
		if( token.Abbreviation("?") ) {
            message="Usage: Assign [manual | onetoone | automatic]";
			PrintMessage();
			
        }
        else if (token.Abbreviation("Manual")) {
            int ntax=taxa->GetNumTaxonLabels();
            message="Number of taxa is ";
            message+=ntax;
            PrintMessage();
            convertsamplestospecies.clear();
            for (int taxon=0; taxon<ntax; taxon++) {
                int assignment;
                while ((cout<<"Taxon "<<taxa->GetTaxonLabel(taxon)<<" Assign to species# = ") && (!(cin >> assignment) || assignment < 1 || assignment > ntax)) {
                    cout << "Incorrect input: enter a number between 1 and "<<ntax<<": ";
                    cin.clear();
                    //cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }
                convertsamplestospecies.push_back(assignment);
				
            }
            bool changedmaxnum=false;
            for (int i=0;i<convertsamplestospecies.size();i++) {
                if (convertsamplestospecies[i]>maxnumspecies) {
                    maxnumspecies=convertsamplestospecies[i];
                    changedmaxnum=true;
                }
            }
            if (changedmaxnum) {
                maxnumspecies++;
                message="Increased the maximum number of species allowed to ";
                message+=maxnumspecies;
                PrintMessage();
            }
        }
        else if (token.Abbreviation("Onetoone")) {
            message="Now assigning each sample to its own species";
            PrintMessage();
            convertsamplestospecies.clear();
            int ntax=taxa->GetNumTaxonLabels();
            for (int i=0; i<ntax; i++) {
                convertsamplestospecies.push_back(i+1);
            }
            bool changedmaxnum=false;
            for (int i=0;i<convertsamplestospecies.size();i++) {
                if (convertsamplestospecies[i]>maxnumspecies) {
                    maxnumspecies=convertsamplestospecies[i];
                    changedmaxnum=true;
                }
            }
            if (changedmaxnum) {
                maxnumspecies++;
                message="Increased the maximum number of species allowed to ";
                message+=maxnumspecies;
                PrintMessage();
            }
        }
        else if (token.Abbreviation("Automatic")) {
            message="Now using the first two symbols in the taxon names to assign species";
            PrintMessage();
            map<nxsstring, int> assignmentmap;
            convertsamplestospecies.clear();
            int ntax=taxa->GetNumTaxonLabels();
            int nspecies=0;
            for (int taxon=0; taxon<ntax; taxon++) {
                nxsstring maplabel=taxa->GetTaxonLabel(taxon);
                message="Taxon ";
                message+=maplabel;
                message+=" abbreviated to ";
                nxsstring shortenedlabel;
                char* s =new char[3];
                strncpy( s, maplabel.c_str(), 2 );
                s[2]='\0';
                shortenedlabel = s;
                delete(s);
                message+=shortenedlabel;
                //message+=" or ";
                //message+=maplabel.ShortenTo(5);
                message+=" assigned to species: taxon";

                if (assignmentmap.count(shortenedlabel)==0) { //No match
                    nspecies++;
                    assignmentmap[shortenedlabel]=nspecies;
                    convertsamplestospecies.push_back(nspecies);
                    message+=nspecies;
                    PrintMessage();
                }
                else { //Already in there
                    typedef map<nxsstring, int>::const_iterator ITER;
                    ITER found=assignmentmap.find(shortenedlabel);
                    int assignedtospeciesnum=found->second;
                    convertsamplestospecies.push_back(assignedtospeciesnum);
                    message+=assignedtospeciesnum;
                    PrintMessage();
                }
            }
            bool changedmaxnum=false;
            for (int i=0;i<convertsamplestospecies.size();i++) {
                if (convertsamplestospecies[i]>maxnumspecies) {
                    maxnumspecies=convertsamplestospecies[i];
                    changedmaxnum=true;
                }
            }
            if (changedmaxnum) {
                maxnumspecies++;
                message="Increased the maximum number of species allowed to ";
                message+=maxnumspecies;
                PrintMessage();
            }

        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading Assign command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }

}

char* BROWNIE::OutputForGTP(ContainingTree *SpeciesTreePtr)
{
    nxsstring TreeFile="#nexus\nbegin trees;";
    //first the species tree
    //SpeciesTree.RandomTree(taxa->GetNumTaxonLabels());
    TreeFile+="\ntree speciestree = ";
    TreeFile+=PipeSpeciesTree(SpeciesTreePtr);
    //now the gene trees
    int originalchosentree=chosentree;
    for (chosentree = 0; chosentree < trees->GetNumTrees(); chosentree++) {
		      Tree t = intrees.GetIthTree(chosentree);
		      TreeFile+="\ntree ";
                      //  if (t.GetName() != "")
                      //		  TreeFile+=NEXUSString (t.GetName());
                      //else
                      //{
                      TreeFile+="genetree_";
                          TreeFile+=chosentree+1;
                          //}
                          TreeFile+=" = ";
                          // Tree
                          TreeFile+=PipeGTP(t);
    }
    chosentree=originalchosentree;
    TreeFile+="\nend;";
    char* OutputForGTPChar=(char*)TreeFile.c_str();
    // if (debugmode) {
    // cout<<"OutputForGTP output =\n"<<OutputForGTPChar<<endl<<endl;
    // }
    return OutputForGTPChar;
}

//Take the current best trees vector, format each tree and copy it into a new vector
vector<ContainingTree> BROWNIE::MakePrettyForOutput()
{
    vector<ContainingTree> PrettyTrees;
    for (int i=0; i<RawBestTrees.size(); i++) {
        ContainingTree PrettyTree=RawBestTrees[i];
        PrettyTree.FindAndSetRoot();
        PrettyTree.Update();
        PrettyTree.GetNodeDepths();
        NodeIterator <Node> n (PrettyTree.GetRoot());
        cur = n.begin();
        while (cur) {
            if (cur->IsLeaf()) {
                nxsstring newlabel=PipeLeafFinalSpeciesTree();
                cur->SetLabel(newlabel);
            }
            else {
                cur->SetLabel("");
            }
            cur = n.next();
        }
        PrettyTrees.push_back(PrettyTree);
    }
    return PrettyTrees;
}

void BROWNIE::FormatAndStoreBestTree(ContainingTree *NewBestTree,vector<double> scorevector)
{
    ContainingTree FormattedNewBestTree=*NewBestTree;
    FormattedNewBestTree.FindAndSetRoot();
    FormattedNewBestTree.Update();
    FormattedNewBestTree.GetNodeDepths();
	if (useCOAL || useMS) {
		FormattedNewBestTree.SetEdgeLengths(true);	
	}
    NodeIterator <Node> n (FormattedNewBestTree.GetRoot());
    cur = n.begin();
    while (cur) {
        if (cur->IsLeaf()) {
            nxsstring newlabel=PipeLeafFinalSpeciesTree();
            cur->SetLabel(newlabel);
        }
        else {
            cur->SetLabel("");
        }
        cur = n.next();
    }
    bool newtree=true;
    for (int j=0; j<FormattedBestTrees.size(); j++) {
        ContainingTree t1=FormattedNewBestTree;
        ContainingTree t2=FormattedBestTrees[j];
        int ntaxt1=t1.GetNumLeaves();
        int ntaxt2=t2.GetNumLeaves();
        int ntaxincommon=PrepareTreesForTriplet(&t1,&t2);
        if (ntaxincommon==ntaxt1 && ntaxincommon==ntaxt2) {
            if (ntaxincommon>=3) { //so, all leaves are in common, but is the topology the same? If there are only two or one species, yes, so newtree=false. If there are at least 3 taxa, there are different possible topologies, so maybe, thus check below
                vector<int> tripletoverlapoutput=GetTripletOverlap(&t1,&t2,ntaxincommon);
                int numberdisagree=tripletoverlapoutput[1];
                if (numberdisagree==0) {
                    newtree=false;
                }
            }
            else {
                newtree=false; //all leaves in common, must have same topology, so same tree
            }
        }
    }
	if ((useCOAL || useMS) && exportalltrees) {
		newtree=true; //since brlen might be different
	}
    if (newtree) {
        FormattedBestTrees.push_back(FormattedNewBestTree);
		BestConversions.push_back(convertsamplestospecies);
		TotalScores.push_back(scorevector[0]);
		GTPScores.push_back(scorevector[1]);
		StructScores.push_back(scorevector[2]);
		if (useCOAL && ContourSearchVector.size()>0 && contourBrlenToExport>0) {
			vector<double> particularcombo=ContourSearchVector[0];
			int numberofbranches=-1+particularcombo.size();
			nxsstring ContourSearchStringBest="";
			nxsstring ContourSearchStringOther="";
			for (int branchcombination=0; branchcombination<ContourSearchVector.size(); branchcombination++) {
				if (ContourSearchVector[branchcombination][numberofbranches]<=0) { //if the score is the best
					ContourSearchStringBest+="\n";
					ContourSearchStringBest+=ContourSearchVector[branchcombination][numberofbranches];
					for (int branch=0; branch<numberofbranches; branch++) {
						ContourSearchStringBest+="\t";
						ContourSearchStringBest+=ContourSearchVector[branchcombination][branch];
					}
				}
				else {
					ContourSearchStringOther+="\n";
					ContourSearchStringOther+=ContourSearchVector[branchcombination][numberofbranches];
					for (int branch=0; branch<numberofbranches; branch++) {
						ContourSearchStringOther+="\t";
						ContourSearchStringOther+=ContourSearchVector[branchcombination][branch];
					}
				}
			}
			if (contourBrlenToExport==1) {
				ContourSearchDescription.push_back(ContourSearchStringBest);
			}
			else if (contourBrlenToExport==2) {
				ContourSearchStringBest+=ContourSearchStringOther;
				ContourSearchDescription.push_back(ContourSearchStringBest);
			}
		}
    }
    for (int k=0; k<FormattedBestTrees.size(); k++) {
        (FormattedBestTrees[k]).Update();
        (FormattedBestTrees[k]).GetNodeDepths();
        (FormattedBestTrees[k]).ClearInternalLabels();
    }
    if (newtree) {
		if (jackknifesearch) {
				jackknifetreestooutput="";
		}
        ofstream outtreef;
        outtreef.open(treefilename.c_str());
        outtreef<<"#nexus\nbegin trees;\n";
        //outtreef<<"[heuristic search, best results after replicate "<<replicate<<"\nSearch options: \n]\n";
		
        for (int i=0; i<FormattedBestTrees.size(); i++) {
			outtreef<<"[ assignment vector: (";
			int numspecies=0;
			for (int j=0;j<convertsamplestospecies.size();j++) {
				outtreef<<" "<<BestConversions[i][j];
				numspecies=GSL_MAX(numspecies,BestConversions[i][j]);
			}
			outtreef<<" ) ]\n";
			if (useCOAL && contourBrlenToExport && ContourSearchDescription.size()>i) {
				outtreef<<"["<<endl<<ContourSearchDescription[i]<<endl<<"]"<<endl;
			}
			if (jackknifesearch) {
				jackknifetreestooutput+="tree jackrep";
				jackknifetreestooutput+=jackrep;
				jackknifetreestooutput+=".";
				jackknifetreestooutput+=i+1;
				jackknifetreestooutput+=" = [&W 1/";
				int ntrees=FormattedBestTrees.size();
				jackknifetreestooutput+=ntrees;
				jackknifetreestooutput+="] [&R] ";
				jackknifetreestooutput+=ReturnFinalSpeciesTree(FormattedBestTrees[i]);
				jackknifetreestooutput+="\n";
			}
            if (!useCOAL && !useMS) {
				outtreef<<"tree sptre"<<i+1<<" = [ "<<"Number of species: "<<numspecies<<"; Score: "<<TotalScores[i]<<" = "<<(1-structwt)<<" x GTP ("<<GTPScores[i]<<") + "<<structwt<<" x Struct ("<<StructScores[i]<<") ] [&R] ";
			}
			else {
				outtreef<<"tree sptre"<<i+1<<" = [ "<<"Number of species: "<<numspecies;
				if (COALaicmode==0) {
					outtreef<<"; NegLnL: ";
				}
				else if (COALaicmode==1) {
					outtreef<<"; AIC: ";
				}
				else {
					outtreef<<"; AICc: ";
				}
				
				outtreef<<TotalScores[i]<<" ] [&R] ";
			}
            outtreef<<ReturnFinalSpeciesTree(FormattedBestTrees[i]);
            outtreef<<endl;
        }
        outtreef<<"end;";
        outtreef.close();
		
		//output for mesquite
		nxsstring mesquitefilename=treefilename;
		mesquitefilename+=".mesquite.nex";
		ofstream outmesqf;
		outmesqf.open(mesquitefilename.c_str());
		outmesqf<<"#nexus\n";
		
		outmesqf<<"\nbegin taxa;\ntitle SamplesBlock;\ndimensions ntax="<<taxa->GetNumTaxonLabels()<<";\ntaxlabels\n\t"; //only need one samples block and one samples tree block
		for (int k=0;k<taxa->GetNumTaxonLabels();k++) {
			outmesqf<<taxa->GetTaxonLabel(k)<<" ";
		}
		outmesqf<<"\n;\nend;\n";

		outmesqf<<"\nbegin trees;\nTitle 'Trees of individual loci';\nLINK Taxa=SamplesBlock;\n";
		//output gene trees
		for (int selectedtree = 0; selectedtree < intrees.GetNumTrees(); selectedtree++) {
			outmesqf<<"tree genetree"<<selectedtree+1<<" = ";
			Tree t=intrees.GetIthTree(selectedtree);
			t.Write(outmesqf);
			outmesqf<<endl;
		}
		//outmesqf<<"OUTPUT GENE TREES\n";
		outmesqf<<"end;\n";
		
		for (int i=0; i<FormattedBestTrees.size(); i++) {
			int numspecies=0;
			for (int j=0;j<convertsamplestospecies.size();j++) {
				numspecies=GSL_MAX(numspecies,BestConversions[i][j]);
			}
			outmesqf<<"\nbegin taxa;\ntitle SpeciesBlock"<<i<<";\ndimensions ntax="<<numspecies<<";\ntaxlabels\n\t";
			for (int k=1;k<=numspecies;k++) {
					outmesqf<<"taxon"<<k<<" ";
			}
			outmesqf<<"\n;\nend;\n";
			
						
			outmesqf<<"\nbegin trees;\ntitle 'Species tree block "<<i<<"';\nLINK Taxa=SpeciesBlock"<<i<<";\n";
			//output tree
			outmesqf<<"tree sptree"<<i+1<<" = ";
			RawBestTrees[i].Write(outmesqf);
			outmesqf<<"\n";
			outmesqf<<"end;\n";
			
			outmesqf<<"\nbegin taxaassociation;\ntitle Assignment"<<i<<";\nTaxa SpeciesBlock"<<i<<" , SamplesBlock;\nAssociates\n\n";
			for (int k=1;k<=numspecies;k++) {
				outmesqf<<"taxon"<<k<<" / ";
				for(int l=0;l<convertsamplestospecies.size();l++) {
					if (convertsamplestospecies[l]==k) {
						outmesqf<<GetTaxonLabel(l)<<" ";
					}
				}
				if (k<numspecies) {
					outmesqf<<",";
				}
				outmesqf<<"\n";
			}
			outmesqf<<";\nend;\n";
			
        }
		
    }
}



double BROWNIE::GetGTPScoreNew(ContainingTree *SpeciesTreePtr)
{
    //cout<<"SpeciesTree"<<endl;
    //(*SpeciesTreePtr).Write(cout);
    //cout<<endl;
    //uses algorithm from Zmasek & Eddy 2001, with modification by Sanderson for inferring only "strong" duplications (those that haven't happened after the last speciation event)
    SpeciesTreePtr->SetLeafNumbers();
    //cout<<"has set leafnumbers\n";
    PreorderIterator <Node> n (SpeciesTreePtr->GetRoot());
    //cout<<"set preorderiterator\n";
    //cout<<"Species Tree = "<<endl;
    // SpeciesTreePtr->ReportTreeHealth();
    //FormatAndStoreBestTree(SpeciesTreePtr);
    //SpeciesTreePtr->Update();
    //SpeciesTreePtr->Draw(cout);
    //cout<<ReturnFinalSpeciesTree(FormattedBestTrees.back());
    //FormattedBestTrees.pop_back();
    //cout<<endl;


    NodePtr currentnode = n.begin();
    //cout<<"    NodePtr currentnode = n.begin();\n";
    assert(currentnode!=NULL);
    //cout<<"beginning to update\n";
    SpeciesTreePtr->Update();
    //cout<<"finished updating\n";
    //cout<<endl;
    //SpeciesTreePtr->ReportTreeHealth();
    //cout<<endl<<"finished reporttreehealth"<<endl;
    int labelcount=1;
    while (currentnode)
    {
        //	cout<<"in while(currentnode)\n";
        //	cout<<"currentnode="<<currentnode<<endl;
        //	cout<<"currentnode anc = "<<currentnode->GetAnc()<<endl;
        //	cout<<"current node index = "<<currentnode->GetIndex()<<endl;
        currentnode->SetIndex(labelcount);
        //	cout<<"currentnode->SetIndex(labelcount)\n";
        labelcount++;
        //	cout<<"labelcount++\n";
        currentnode = n.next();
        //	cout<<"currentnode = n.next();\n";
    }
    //cout<<"out of while(currentnode)\n";
    int originalchosentree=chosentree;
    double weightednumDup=0;
    for (int selectedtree = 0; selectedtree < intrees.GetNumTrees(); selectedtree++) {
		if (gtptoohigh) {
			break;
		}
		bool usethistree=true;
		if (jackknifesearch) {
			if (jackknifevector[selectedtree]==0) {
				usethistree=false;
			}
		}
		if (usethistree) {
			int numDup=0;
			Tree t=intrees.GetIthTree(selectedtree);
			
        //t.Write(cout);
        //   cout<<endl;
        //cout<<"Gene Tree = "<<endl;
			t.Update();
        // t.Draw(cout);
        //  cout<<endl;
			NodeIterator <Node> m (t.GetRoot());
			NodePtr genetreenode=m.begin();
			NodePtr a;
			NodePtr b;
			NodePtr g1;
			NodePtr g2;
			NodePtr Mg1;
			NodePtr Mg2;
			map<NodePtr,NodePtr> GtoSmap;
			map<NodePtr,NodePtr>::iterator GtoSmapiter;
        //map<NodePtr,map<int, int> > setsM;
        //   map<int,int> g1map;
        //   map<int,int> g2map;
        //   map<int,int> a;
        // map<int,int> b;
        //map<NodePtr,map<int, int> >::iterator setsMiterator;
			while (genetreenode!=NULL) {
				if (gtptoohigh) {
					break;
				}				
				genetreenode->SetMarked(false); //in this context, whether it predates a speciation event and so could be a strong duplication
                                            //  cout<<"\ngenetreenode = "<<genetreenode<<" is leaf? "<<genetreenode->IsLeaf()<<" "<<genetreenode->GetLabel()<<endl;
				if (!(genetreenode->IsLeaf())) {
					g1=genetreenode->GetChild();
					g2=g1->GetSibling();
					if(g1->IsLeaf()) {
						int SampleNumber=taxa->FindTaxon(g1->GetLabel());
                    //cout<<"samplenum = "<<SampleNumber<<" convertsamplestospecies.size() = "<<convertsamplestospecies.size()<<endl;
                    //   tempmap.insert(make_pair(convertsamplestospecies[SampleNumber],1));
                    //  setsM.insert(make_pair(g1,tempmap));
						if (SampleNumber>=convertsamplestospecies.size()) {
							cout<<"SampleNumber "<<SampleNumber<<" but...";
							cout<<" convertsamplestospecies.size() = "<<convertsamplestospecies.size()<<endl;
							for(int i=0;i<convertsamplestospecies.size();i++) {
								cout<<convertsamplestospecies[i]<<" ";
							}
							cout<<endl;
						}
						assert(SampleNumber<convertsamplestospecies.size());
						GtoSmap.insert(make_pair(g1,SpeciesTreePtr->GetLeafWithNumber(convertsamplestospecies[SampleNumber])));
                    //  cout<<"Found leaf "<<SpeciesTreePtr->GetLeafWithNumber(convertsamplestospecies[SampleNumber]);
                    //  cout<<" "<<(SpeciesTreePtr->GetLeafWithNumber(convertsamplestospecies[SampleNumber]))->GetLabel()<<endl;
					}
					if(g2->IsLeaf()) {
						int SampleNumber=taxa->FindTaxon(g2->GetLabel());
						GtoSmap.insert(make_pair(g2,SpeciesTreePtr->GetLeafWithNumber(convertsamplestospecies[SampleNumber])));
                    //  cout<<"Found leaf "<<SpeciesTreePtr->GetLeafWithNumber(convertsamplestospecies[SampleNumber]);
                    //  cout<<" "<<(SpeciesTreePtr->GetLeafWithNumber(convertsamplestospecies[SampleNumber]))->GetLabel()<<endl;
						
					}
					if(g1->IsMarked() || g2->IsMarked()) {
						genetreenode->SetMarked(true); // it predates a speciation event, so can be a strong dup (idea from Sanderson)
					}
					GtoSmapiter=GtoSmap.find(g1);
					a=(*GtoSmapiter).second;
					Mg1=a;
					GtoSmapiter=GtoSmap.find(g2);
					b=(*GtoSmapiter).second;
					Mg2=b;
                //cout<<"a = "<<a<<" b = "<<b<<endl;
					while (a->GetIndex() != b->GetIndex()) {
						if (a->GetIndex() > b->GetIndex()) {
							a=a->GetAnc();
						}
						else {
							b=b->GetAnc();
						}
					}
					GtoSmap.insert(make_pair(genetreenode,a));
                //  cout<<"Mgenetreenode,Mg1,Mg2,genetreenode,g1,g2\n"<<a<<"\t"<<Mg1<<"\t"<<Mg2<<"\t"<<genetreenode<<"\t"<<g1<<"\t"<<g2<<endl;
					if((a==Mg1) || (a==Mg2)) {
                    //   cout<<"Dup";
						if(genetreenode->IsMarked()) {
                        //    cout<<"--Strong, will count"<<endl;
							numDup++;
                        //    cout<<"Match, numDup = "<<numDup<<endl;
						}
					}
					else {
						genetreenode->SetMarked(true); //since it's a speciation event
                                                   //  cout<<"SpeciationEvent"<<endl;
					}
				}
				genetreenode=m.next();
			}
        //     cout<<"numDup="<<numDup<<" weight = "<<trees->GetTreeWeight(chosentree)<<endl;
			weightednumDup+=(trees->GetTreeWeight(chosentree))*(1.0*numDup);
			if (weightednumDup*(1.0-structwt)>bestscorelocal) {
				weightednumDup=(0.0001+bestscorelocal)/(1.0-structwt);
				gtptoohigh=true;
				break;
			}
		}
    }
    chosentree=originalchosentree;
    //cout<<"weightednumDup="<<weightednumDup<<endl;
    return weightednumDup;
}

//took this function out as no longer depend on external gtp
//double BROWNIE::GetGTPScore(ContainingTree *SpeciesTreePtr)
//{
//    double totalscore=0;
//    badgtpcount=0;
//    if (SpeciesTreePtr->GetNumLeaves()>1) { //don't bother getting the GTP score if there's only one species
//        int originalchosentree=chosentree;
//        if (outputallatonce) {
//            totalscore=GSL_NEGINF;
//            badgtpcount=-1;
//            while (totalscore<0) {
//                totalscore=ReturnScore(OutputForGTP(SpeciesTreePtr),unrooted);
//                badgtpcount++;
//                if (badgtpcount==500) {
//                    cout<<"Potential error getting tree score, trying to recover...";
//                }
//                if (badgtpcount>10000) {
//                    cout<<"Has failed, aborting. Email brownie@brianomeara.info with this message, also include the info below:\n\n";
//                    cout<<OutputForGTP(SpeciesTreePtr)<<endl;
//                    errormsg="Aborted";
//                    throw XNexus( errormsg);
//                }
//            }
//            if (badgtpcount>=500) {
//                cout<<"... fixed! Required "<<badgtpcount<<" tries\n";
//            }
//        }
//       else {
//            for (chosentree = 0; chosentree < trees->GetNumTrees(); chosentree++) {
//                double thisscore=GSL_NEGINF;
//                while (thisscore<0) {
//                    nxsstring TreeFile="#nexus\nbegin trees;";
//                    TreeFile+="\ntree speciestree = ";
//                    TreeFile+=PipeSpeciesTree(SpeciesTreePtr);
//                    Tree t = intrees.GetIthTree(chosentree);
//                    TreeFile+="\ntree ";
//                    TreeFile+="genetree_";
//                    TreeFile+=chosentree+1;
//                    TreeFile+=" = ";
//                    TreeFile+=PipeGTP(t);
//                    TreeFile+="\nend;";
//                    thisscore=ReturnScore((char*)TreeFile.c_str(),unrooted);
//                    if (thisscore<0) {
//                        badgtpcount++;
//                    }
//                 }
//                totalscore+=(trees->GetTreeWeight(chosentree))*thisscore; //Get properly weighted score
//            }
//        }
//        chosentree=originalchosentree;
//    }
//return totalscore;
//}


/**
* @method OutputComment [virtual void:public]
 * @param msg [nxsstring&] the output comment to be displayed
 *
 * This function is called whenever an output comment (i.e., a comment
                                                       * beginning with an exclamation point) is found in the data file.
 */
void BROWNIE::OutputComment( nxsstring& msg )
{
    message=msg;
    PrintMessage();
}


/**
* @method Read [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XNexus
 *
 * This function provides the ability to read everything following
 * the block name (which is read by the Nexus object) to the end or
 * endblock statement. Characters are read from the input stream
 * in. Overrides the pure virtual function in the base class.
 */
void BROWNIE::Read( NexusToken& token )
{
    isEmpty = false;
	//need to  have these called here so that any brownie commands using characters work properly
	if(!characters->IsEmpty() ) {
		if (characters->GetDataType()==6) {
			continuouscharacters=characters;
			//cout<<"Found continuous characters\n";
		}
		else {
			discretecharacters=characters;
			numbercharstates=discretecharacters->GetMaxObsNumStates();
			localnumbercharstates=numbercharstates;
			//cout<<"Found discrete characters\n";
		}
	}
	
	if(!characters2->IsEmpty() ) {
		if (characters2->GetDataType()==6) {
			continuouscharacters=characters2;
			//cout<<"Found continuous characters\n";
		}
		else {
			discretecharacters=characters2;
			numbercharstates=discretecharacters->GetMaxObsNumStates();	
			localnumbercharstates=numbercharstates;
			//cout<<"Found discrete characters\n";
		}
	}
	
    // this should be the semicolon after the block name
    //
    token.GetNextToken();

    if( !token.Equals(";") ) {
        errormsg = "Expecting ';' after ";
        errormsg += id;
        errormsg += " block name, but found ";
        errormsg += token.GetToken();
        errormsg += " instead";
        throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
    }

    for(;;)
    {
        token.GetNextToken();

        if( token.Abbreviation("ENdblock") ) {
            HandleEndblock( token );
            break;
        }
        else if( token.Abbreviation("Help") ) {
            HandleHelp( token );
        }
        else if( token.Abbreviation("Log") ) {
            HandleLog( token );
        }
        else if( token.Abbreviation("ECho") ) {
            HandleEcho( token );
        }
        else if( token.Abbreviation("EXecute") ) {
            HandleExecute( token );
        }
        else if( token.Abbreviation("Gettrees") ) {
            HandleGettrees( token );
        }
        else if( token.Abbreviation("CHoose") ) {
            HandleChoose( token );
        }
        else if( token.Abbreviation("SEt") ) {
            HandleSet( token );
        }
        else if( token.Abbreviation("COMpare") ) {
            HandleCompareRandomTrees( token );
        }
        else if( token.Abbreviation("HSearch") ) {
            HandleHeuristicSearch( token );
        }
        else if( token.Abbreviation("HEUristicsearch") ) {
            HandleHeuristicSearch( token );
        }
		else if( token.Abbreviation("DETTman") ) {
            HandleDettmanCollapse( token );
        }
		else if( token.Abbreviation("JAckknife") ) {
			jackknifesearch=true;
			HandleHeuristicSearch( token ) ;  //a jackknife search is just like a heuristic search, just with jackknife options
		}
		else if( token.Abbreviation("EXhaustivesearch") ) {
			DoExhaustiveSearch();
		}
		else if( token.Abbreviation("NAst") ) {
			HandleNast(token);
		}
        else if( token.Abbreviation("ASsign") ) {
            Assign( token );
        }
        else if( token.Abbreviation("CItation") ) {
            HandleCitation( token );
        }
		else if( token.Abbreviation("ACcuracy") ) {
            HandleAccuracy( token );
        }
		else if( token.Abbreviation("PARtitionededgesupport") ) {
            HandlePartitionedEdgeSupport( token );
        }
        else if( token.Abbreviation("SHowtree") ) {
            HandleShowtree( token );
        }
        else if( token.Abbreviation("Blocks") ) {
            HandleBlocks( token );
        }
        else if( token.Abbreviation("VCV") ) {
            HandleVCV( token );
        }
        else if( token.Abbreviation("TIPVALues") ) {
            HandleTipValues( token );
        }
        else if( token.Abbreviation("TIPVARiance") ) {
            HandleTipVariance( token );
        }
        else if( token.Abbreviation("MOdel") ) {
            HandleModel( token );
        }
		else if( token.Abbreviation("DIscrete") ) {
			HandleDiscrete( token );
		}
        else if( token.Abbreviation("PAGELcorrelation") ) {
			HandlePagelDiscrete( token );
		}
		else if( token.Abbreviation("SIMulate") ) {
			HandleSimulateCharacters( token );
		}
		else if( token.Abbreviation("LOss") ) {
			HandleLoss( token );
		}
        else if( token.Abbreviation("TIMeslice") ) {
            HandleTimeSlice( token );
        }
        else if( token.Abbreviation("Debug") ) {
            HandleDebug( token );
        }
		else if( token.Abbreviation("NOQuitonerror") ) {
            HandleNoQuitOnErr( token );
        }
        else if (token.Abbreviation("PREorder") ) {
            PreOrderTraversal(token);
        }
        else if (token.Abbreviation("MRCA") ) {
            HandleMRCA( token );
        }
        else if (token.Abbreviation("GREP") ) {
            HandleGrepCount( token );
        }
        else if (token.Abbreviation("EXport") ) {
            HandleExport(token);
        }
        else if ( token.Abbreviation("OPtimization") || token.Abbreviation("CONtinuous") ) {
            HandleDebugOptimization( token );
        }
        else if ( token.Abbreviation("NUmopt") ) {
            NumOpt( token );
        }
        else if (token.Abbreviation("GARLAND") ) {
            HandleGarland(token);
        }
        else if( token.Abbreviation("RateTest") ) {
            HandleRateTest( token );
        }
		else if( token.Abbreviation("ORderbytree") ) {
			HandleOrderByTree( token );
		}
        else if(token.Abbreviation("Taxset") ) {
            //errormsg = "Sorry, enter taxsets in the ASSUMPTIONS block";
            //throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            assumptions->HandleTaxset( token );
        }
        else if( token.Abbreviation("Printedgelengths") ) {
            HandlePrintEdgeLengths( token);
        }
        else if( token.Abbreviation("Quit") ) {
            quit_now = true;
			PrintCitations();
            message = "\nPlease remember to send bug reports to omeara.brian@gmail.com.\n";
            PrintMessage();
            break;
        }
        else
        {
			if ((token.GetToken()).length()>1) {
				SkippingCommand( token.GetToken() );
				do {
					token.GetNextToken();
				} while( !token.AtEOF() && !token.Equals(";") );
				
				if( token.AtEOF() ) {
					errormsg = "Unexpected end of file encountered";
					throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
				}
			}
			else {
				do {
					token.GetNextToken();
				} while( !token.AtEOF() && !token.Equals(";") );
				break;
			}
        }
    }
}

/**
* @method Reset [void:protected]
 *
 * Overrides the pure virtual function in the base class.
 */
void BROWNIE::Reset()
{
    isEmpty = true;

    inf_open = false;
    quit_now = false;
    message = "";
    next_command[0] = '\0';
}

/**
* @method Report [virtual void:public]
 * @param out [ostream&] the output stream to which to write the report
 *
 * This function outputs a brief report of the contents of this BROWNIE block.
 * Overrides the pure virtual function in the base class.
 */
void BROWNIE::Report( ostream& /* out */ )
{
    message = "";
    PrintMessage();

    message = id;
    message += " block contains...";
    PrintMessage();
}

/**
* @method RunCmdLine [void:public]
 *
 * Runs the command line interpreter, allowing BROWNIE to interact with user.
 */
void BROWNIE::RunCmdLine(bool inputfilegiven, nxsstring fn)
{
    taxa = new TaxaBlock();
    trees = new TreesBlock(*taxa);
    assumptions = new AssumptionsBlock( *taxa );
    characters = new CharactersBlock( *taxa, *assumptions );
	characters2 = new CharactersBlock2 (*taxa, *assumptions );
    Add( taxa );
    Add( trees );
    Add( assumptions );
    Add( characters );
	Add( characters2 );
    Add( this );
    cout<<endl<<endl<<endl;
    cout<<"                               Brownie V2.1.1 (r"<<$SVN_VERSION<<")"<<endl;
	cout<<"                          Character evolution models,"<<endl;
	cout<<"                             species delimitation,"<<endl; 
	cout<<"                               and tree search"<<endl<<endl;
    cout<<"                                Brian O'Meara"<<endl;
    cout<<"                     http://www.brianomeara.info/brownie"<<endl<<endl;
    cout<<"                     Uses Paul Lewis' Nexus Class Library"<<endl;
    cout<<"                (modified to deal with continuous characters)"<<endl;
    cout<<"                        Rod Page's TreeLib & supertree,"<<endl;
    cout<<"                             Mike Sanderson's GTP,"<<endl;
    cout<<"                            and the GNU Scientific Library"<<endl<<endl;
    cout<<"                  Type \"help\" [no quotes] for a list of commands"<<endl;
    cout<<"                  and \"commandname ?\" for help for each command."<<endl<<endl;
    cout<<"                       Compiled on "<<__DATE__<<" at "<<__TIME__<<endl;
	cout<<"                         Using GSL version "<<GSL_VERSION<<endl<<endl;
	
    quit_now = false;
    if (inputfilegiven) {
        HandleExecuteCmdLine(fn);
    }
    while( !quit_now )
    {
        cerr << endl;
        cerr << "brownie> ";
        //cin.getline( next_command, 256 );
        fgets(next_command, COMMAND_MAXLEN,stdin); //use fgets since GCC on OS 10.2 has trouble with cin.
        PreprocessNextCommand();
        HandleNextCommand();
    }
}

/**
* @method Run [void:public]
 *
 * Runs the command line interpreter, allowing BROWNIE to interact with user.
 * Typically, this is the only function called in main after a BROWNIE object
 * is created.
 */
//void BROWNIE::Run()
//{
//	taxa = new TaxaBlock();
//	trees = new TreesBlock(*taxa);
//	assumptions = new AssumptionsBlock( *taxa );
//	characters = new CharactersBlock( *taxa, *assumptions );

    //	Add( taxa );
    //	Add( trees );
    //	Add( assumptions );
    //	Add( characters );
    //	Add( this );
    //      cout<<endl<<endl<<endl<<"                               Brownie V2.0b2"<<endl;
    //        cout<<"                 Testing rates of continuous character evolution"<<endl<<endl;
    //        cout<<"                                Brian O'Meara"<<endl;
    //        cout<<"                     http://www.brianomeara.info/brownie"<<endl<<endl;
    //        cout<<"                     Uses Paul Lewis' Nexus Class Library"<<endl;
    //        cout<<"                (modified to deal with continuous characters)"<<endl;
    //        cout<<"                         Rod Page's TreeLib, and code"<<endl;
    //        cout<<"                      from Ligia Mateiu & John Burkardt."<<endl<<endl;
    //        cout<<"                  Type \"help\" [no quotes] for a list of commands"<<endl;
    //        cout<<"                  and \"commandname ?\" for help for each command."<<endl<<endl;
    //	quit_now = false;
    //	while( !quit_now )
    //	{
    //		cerr << endl;
    //		cerr << "brownie> ";
    //		//cin.getline( next_command, 256 );
    //              fgets(next_command, COMMAND_MAXLEN,stdin); //use fgets since GCC on OS 10.2 has trouble with cin.
    //		PreprocessNextCommand();
    //		HandleNextCommand();
    //	}
    //}

/**
    * @method SkippingBlock [virtual void:public]
 * @param blockName [nxsstring] the unrecognized block name
 *
 * Called when program does not recognize a block name encountered in a
 * NEXUS file.  Virtual function that overrides the pure virtual function
 * in the base class Nexus.
 */
    void BROWNIE::SkippingBlock( nxsstring blockName )
    {
        message = "Skipping unknown block (";
        message += blockName;
        message += ")";
        PrintMessage();
    }

    /**
    * @method SkippingCommand [virtual void:public]
     * @param commandName [nxsstring] the name of the command being skipped
     *
     * This function is called when an unknown command named commandName is
     * about to be skipped.  This version of the function (which is identical
                                                           * to the base class version) does nothing (i.e., no warning is issued
                                                                                                      * that a command was unrecognized).  Modify this virtual function to
     * provide such warnings to the user (or eliminate it altogether since
                                          * the base class version already does what this does).
     */
    void BROWNIE::SkippingCommand( nxsstring commandName )
    {
        message = "Skipping unknown command (";
        message += commandName;
        message += ")\nType help for a list of available commands.";
        PrintMessage();
    }

    /**
    * @method SkippingDisabledBlock [virtual void:public]
     * @param blockName [nxsstring] the name of the block just exited
     *
     * Called by the Nexus object when skipping a block named blockName
     * that has been disabled. Allows program to notify user of progress
     * in parsing the NEXUS file. Virtual function that overrides the
     * pure virtual function in the base class Nexus.
     */
    void BROWNIE::SkippingDisabledBlock( nxsstring /*blockName*/ )
    {
    }

    /**
    * @method TaxonLabelToNumber [int:protected]
     * @param s [nxsstring] the taxon label to be translated to a taxon number
     *
     * The code here is modified from charactersblock.cpp
     */
    int BROWNIE::TaxonLabelToNumber( nxsstring s )
    {
        int i;
        try {
            i = 1 + taxa->FindTaxon(s);
        }
        catch( TaxaBlock::nosuchtaxon ) {
            i = 0;
        }
        return i;
    }

    /**
    * Asks user if "something" is ok, where "something" is expressed
     * in the title and message displayed.  This is a virtual function
     * so it can be overridden in a derived class to use a different
     * (perhaps graphical) means for displaying the message.
     * Note: mb_message should terminate with a quesiton mark; none
     * will be provided by this function.
     */
    bool BROWNIE::UserSaysOk( nxsstring mb_message, nxsstring mb_title )
    {
        cerr << endl;
        cerr << mb_title << endl;
        cerr << "  " << mb_message;
        cerr << " (y/n) ";
        fgets(next_command, COMMAND_MAXLEN,stdin);
        //cin.getline( next_command, COMMAND_MAXLEN ); //GCC problem
        //bool yep  = ( next_command[0] == 'y' && next_command[1] == '\0' );
        //bool nope = ( next_command[0] == 'n' && next_command[1] == '\0' );
        bool yep  = ( next_command[0] == 'y');
        bool nope = ( next_command[0] == 'n');
        int patience=5;
        int loopnumber=0;

        while( !yep && !nope )
        {
            loopnumber++;
            if (loopnumber>patience) {
                errormsg = "ERROR: You're too indecisive.";
                throw XNexus( errormsg);
            }
            cerr << endl;
            cerr << "Must answer by typing either y or n and then pressing the Enter key" << endl;
            cerr << endl;
            cerr << mb_title << endl;
            cerr << "  " << mb_message;
            cerr << " (y/n) ";
            fgets(next_command, COMMAND_MAXLEN,stdin);
            //cin.getline( next_command, COMMAND_MAXLEN ); //OS 10.2 GCC problem
            //yep  = ( next_command[0] == 'y' && next_command[1] == '\0' );
            //nope = ( next_command[0] == 'n' && next_command[1] == '\0' );
            yep  = ( next_command[0] == 'y');
            nope = ( next_command[0] == 'n');
        }

        return yep;
    }

    void BROWNIE::HandleMRCA( NexusToken& token )
    {
        nxsstring mrca_name;
        IntSet taxonnumbers;
        token.GetNextToken();

        if( token.Abbreviation("?") ) {
            message="Usage: MRCA <name> = <taxon 1> <taxon 2> >...>\n\n";
            message+="Loads the name of the most recent common ancestor of a given set of taxa (using names or numbers, as in a taxset definition).\nNote that the MRCA name is case-sensitive.";
            PrintMessage();
            token.GetNextToken();

        }
        else {
            // Token now stored should be the name of a mrca
            mrca_name = token.GetToken();

            // Now grab the equals sign
            token.GetNextToken();
            if( !token.Equals("=") ) {
                errormsg = "Expecting '=' in MRCA definition but found ";
                errormsg += token.GetToken();
                errormsg += " instead";
                throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            int totalTaxa = taxa->GetNumTaxonLabels();
            SetReader( token, totalTaxa, taxonnumbers, *this, SetReader::taxset ).Run();
            MrcaMap[mrca_name]=taxonnumbers;
        }
    }
    
    void BROWNIE::HandleGrepCount( NexusToken& token )
    {
        bool finishexecuting=true;
        bool returnmatches=false;
        bool usefileregex=false;
        bool usestepwise=true;
      	citationarray[3]=true;
     	nxsstring observedtreefile;
     	nxsstring simtreefile;
     	nxsstring assignmentfile;
     	multimap<nxsstring, nxsstring> speciestosamples;
     	multimap<nxsstring, nxsstring> samplestospecies; 
     	map<nxsstring, int> specieslist;
     	map<nxsstring, nxsstring> sampleslist;
     	vector<nxsstring> speciesvector;
     	vector<nxsstring> samplesvector;

     	nxsstring speciesnamenxs;
     	nxsstring samplenamenxs;
     	for(;;)
    		{
			token.GetNextToken();
			if( token.Abbreviation("?") ) {
				message="Usage: GREP  observedtrees=file_of_observedtrees_trees.txt simulatedtrees=file_of_simulated_trees.txt assignments=file_of_assignments.txt returnmatches=no usefileregex=no usestepwise=no\n\n";
			message+="Multiple gene trees may be consistent with a given species tree (i.e, if the species tree is ((A,B),C), gene trees ((((A1,A2,),A3),B1),C1) and ((((A3,A2,),A1),B1),C1) both match. This function takes one or more observed gene trees (in a file), an assignment (tab-delimited text, with the species name followed by a tab and then the gene sample name (i.e, SpeciesA<tab>A1), and makes a string to use with grep that will tell you how many trees in the gene trees file (set of newick trees, one line per tree) are consistent with each observed gene tree, given the assignment of samples to species. Note that passing Brownie's output to 'grep Match' will return only the relevant results. For very big trees, it may be better to set usefileregex=yes so that the grep line is stored in a file (good if you get the error message 'regular expression too big'). Usestepwise=y passes multiple short arguments to a chain of greps instead of one huge regex.";
				PrintMessage();
				finishexecuting=false;
			}
			else if (token.Abbreviation("Returnmatches") ) {
				nxsstring yesnomatches=GetFileName(token);
				if (yesnomatches[0] == 'n') {
					returnmatches=false;
				}
				else {
					returnmatches=true;
				}
			}		
			else if (token.Abbreviation("USEFileregex") ) {
				nxsstring yesnomatches=GetFileName(token);
				if (yesnomatches[0] == 'n') {
					usefileregex=false;
				}
				else {
					usefileregex=true;
				}
			}		
			else if (token.Abbreviation("USEStepwise") ) {
				nxsstring yesnomatches=GetFileName(token);
				if (yesnomatches[0] == 'n') {
					usestepwise=false;
				}
				else {
					usestepwise=true;
				}
			}		
			else if( token.Abbreviation("Simulatedtrees") ) {
				simtreefile = GetFileName(token);
			}
			else if( token.Abbreviation("Observedtrees") ) {
				observedtreefile = GetFileName(token);
			}
			else if( token.Abbreviation("Assignments") ) {
				assignmentfile = GetFileName(token);
			}        
			else if( token.Equals(";") ) {
				if (finishexecuting) {
					Profile<Tree> inObservedTrees;
					 if( FileExists( observedtreefile.c_str() ) )
					 {
						 ifstream inf( observedtreefile.c_str(), ios::binary | ios::in );
						 ifstream observedtreefile_stream;
						 observedtreefile_stream.open(observedtreefile.c_str(),ios::in);
						 if (!inObservedTrees.ReadTrees(observedtreefile_stream))
						 {
							 errormsg="No species trees read from file\n";
							throw XNexus( errormsg);
						 }
						 observedtreefile_stream.close();
					}
					else {
						errormsg="File ";
						errormsg+=observedtreefile;
						errormsg+=" does not exist, at least where Brownie is looking for it";
						throw XNexus(errormsg);
					}
					
					 if( FileExists( assignmentfile.c_str() ) )
					 {
						 //ifstream inf( assignmentfile.c_str(), ios::binary | ios::in );
						 //ifstream assignmentfile_stream;
						 //assignmentfile_stream.open(assignmentfile.c_str(),ios::in);
						 int count=0;
						 FILE *assignmentfile_open = fopen(assignmentfile.c_str(), "r");
						 char speciesname [100];
						 char samplename [100];
						 message="Assignments\n\tNo.\tSpecies\tSample\n";
						 while (fscanf(assignmentfile_open, "%s\t%s", speciesname, samplename) == 2) {
						  speciesnamenxs="";
						  speciesnamenxs+=speciesname;
						  samplenamenxs="";
						  samplenamenxs+=samplename;
						  count=count+1;
						  message+="\t";
						  message+=count;
						  message+="\t";
						  message+=speciesnamenxs;
						  message+="\t";
						  message+=samplenamenxs;
						  message+="\n";
						  speciestosamples.insert(pair<nxsstring, nxsstring>(speciesnamenxs,samplenamenxs));
						  samplestospecies.insert(pair<nxsstring, nxsstring>(samplenamenxs, speciesnamenxs));
						  ++specieslist[speciesnamenxs]; //keeps track of number of samples per species and keys are unique set of species labels
						  sampleslist[samplenamenxs]=speciesnamenxs;
						  samplesvector.push_back(samplenamenxs);
						}
						fclose(assignmentfile_open);
						PrintMessage();
						
						map<nxsstring, int>::const_iterator itr;

						message="Samples per species:";
						for (itr = specieslist.begin(); itr != specieslist.end(); ++itr){
							 speciesvector.push_back(itr->first);
							 message+="\n\t";
							 message+=itr->first;
							 message+="\t";
							 message+=itr->second;
						}
						PrintMessage();
						
						//assignmentfile_stream.close();
					}
					else {
						errormsg="File ";
						errormsg+=assignmentfile;
						errormsg+=" does not exist, at least where Brownie is looking for it";
						throw XNexus(errormsg);
					}				
					
					
					
//						ContainingTree CurrentTree;
//						CurrentTree.SetRoot((intrees.GetIthTree(inputchosentree)).CopyOfSubtree((intrees.GetIthTree(inputchosentree)).GetRoot()));
//						CurrentTree.Update();
//						CurrentTree.GetNodeDepths();
//						int numspecies=CurrentTree.GetNumLeaves();		
//						int ntax=taxa->GetNumTaxonLabels();
						vector <int> samplesperspecies;
						vector <nxsstring> labelswithinspecies;
						int numsamples=speciestosamples.size();
						int numspecies=specieslist.size();

						double numberofpermutations=1;
						int previoustotal=0;
						for (int currentspecies=0; currentspecies<numspecies; currentspecies++) {
							int samplecountthisspecies=0;
							nxsstring labelregex="\\(";
							for (int currentsample=0; currentsample<numsamples; currentsample++) {
								if (sampleslist[ samplesvector[currentsample] ]==speciesvector[currentspecies]) {
									samplecountthisspecies++;
									if (samplecountthisspecies>1) {
									   labelregex+="\\|";
									}
									labelregex+=samplesvector[currentsample];
								}
							}
							labelregex+="\\):[0-9]*.[0-9]*"; //ms exports brlen, too
							labelswithinspecies.push_back(labelregex);
							previoustotal+=samplecountthisspecies;
							numberofpermutations*=gsl_sf_fact(samplecountthisspecies);
						}
				
					
					
						//rather than looking for exact match, get probabilities of given topology with all possible permutations of labels of samples from a given species onto ms taxon numbers, then divide by number of such permutations
						
						//loop over all the gene trees
						message="\nGrep match results\n\tTree\tNumber of matches";
						PrintMessage();
						for (int chosentreenum=0; chosentreenum<inObservedTrees.GetNumTrees(); chosentreenum++) { //loop over all the observed gene trees
							nxsstring grepstring="";
							nxsstring grepfilestring="";
							nxsstring grepstringreturnmatch="";
							nxsstring grepstepwise="";
							int numberintraspecificcherries=0;
							Tree CurrentGeneTreeTreeFmt=inObservedTrees.GetIthTree(chosentreenum);
							NodeIterator <Node> n (CurrentGeneTreeTreeFmt.GetRoot());
							//traverse once to count intraspecific sister cherries
							cur = n.begin();
							while (cur) {
								if (cur->IsLeaf()) {
									NodePtr sister=cur->GetSibling();
									if (sister!=NULL) {
										if (sister->IsLeaf()) {
											//cout<<cur->GetLabel()<<" with "<<sister->GetLabel()<<endl;
											
											if (sampleslist[cur->GetLabel()] == sampleslist[sister->GetLabel()]) { //they have the same species (though node labels are changed, we do the one that's a child before its sib)
											//if (convertsamplestospecies[taxa->FindTaxon(cur->GetLabel())]==convertsamplestospecies[taxa->FindTaxon(sister->GetLabel())]) { //they have the same species (though node labels are changed, we do the one that's a child before its sib)
												numberintraspecificcherries++;
											}
										}
									}
								}
								cur = n.next();
							}
							
							//now get the grep lines
							cur = n.begin();
							while (cur) {
								if (cur->IsLeaf()) {
									nxsstring specieslabel=sampleslist[cur->GetLabel()];
									nxsstring greplabel="ERROR";
									for (int currentspecies=0; currentspecies<numspecies; currentspecies++) {
										if (speciesvector[currentspecies]==specieslabel) {
											greplabel=labelswithinspecies[currentspecies];
										}
									}
									cur->SetLabel(greplabel); //gets regex for tip labels
								}
								else {
									nxsstring combinedstring="";
									nxsstring leftchild=(cur->GetChild())->GetLabel();
									nxsstring rightchild=((cur->GetChild())->GetSibling())->GetLabel();
									nxsstring additionalregex="(\\(";
									additionalregex+=leftchild;
									additionalregex+=",";
									additionalregex+=rightchild;
									additionalregex+="\\|";
									additionalregex+=rightchild;
									additionalregex+=",";
									additionalregex+=leftchild;
									additionalregex+="\\))";
									if (cur->GetAnc()!=NULL) { // not root
										additionalregex+=":*[0-9]*.[0-9]*";
										grepstepwise+= " | grep \"";
										grepstepwise+=additionalregex;
										grepstepwise+="\"";
									}
									else { //originally, did this at each node, with the idea that from MS, you first filter out the lines that don't match a cherry, then filter from that filtered set, etc., thinking it would be faster than just using the regex at the root. Actually, just using the regex at the root was faster, so I do that.
										grepstring+=" | grep -c \"";
										grepstringreturnmatch+=" | grep \"";
										grepstring+=additionalregex;
										grepstringreturnmatch+=additionalregex;
										grepfilestring+=additionalregex;
										grepstring+="\"";
										grepstringreturnmatch+="\"";
										grepstepwise+=" | grep -c ';'";
									}
									cur->SetLabel(additionalregex);
								}
								cur = n.next();
							}
							/*ofstream newickf;
							newickf.open("tmp_newick.nwk");
							CurrentGeneTreeTreeFmt.Write(newickf);
							newickf.close();
							nxsstring finalsystemcall="cat tmp_newick.nwk";*/
							nxsstring finalsystemcall="cat ";
							finalsystemcall+=simtreefile;
							if (!usefileregex) {
								if (!usestepwise) {
									finalsystemcall+=grepstring;
								}
								else {
									finalsystemcall+=grepstepwise;
								}
							}
							else {
								nxsstring grepregex="tmp_grepregex.txt";
								finalsystemcall+=" | grep -c -f tmp_grepregex.txt";
								ofstream grepregexf;
								grepregexf.open(grepregex.c_str(), ios::out | ios::app );
								grepregexf<<grepfilestring;
								grepregexf.close();
							}
							nxsstring msinputfile="tmp_mscount.txt";
							finalsystemcall+=" > ";
							finalsystemcall+=msinputfile;
							//system("rm mscount.txt");
							int returncode=system(finalsystemcall.c_str());
							if (debugmode) {
								if (!usestepwise) {
									cout<<"grep line is "<<endl<<grepstring<<endl;
								}
								else {
									cout<<"grep line is "<<endl<<grepstepwise<<endl;

								}
								cout<<"return code is "<<returncode<<" (divided by 256, is "<<returncode/256<<")"<<endl;
							}
							ifstream msin;
							msin.open( msinputfile.c_str(), ios::binary | ios::in );
							if (msin) {
								char inputitem [COMMAND_MAXLEN];
								msin>>inputitem;
								double numbermatches=atoi(inputitem);
								message="Match\t";
								message+=chosentreenum+1;
								message+="\t";
								message+=numbermatches;
								PrintMessage();
								msin.close();
								if (returnmatches) {
									nxsstring finalsystemcallmatch="cat ";
									finalsystemcallmatch+=simtreefile;
									finalsystemcallmatch+=grepstringreturnmatch;
									finalsystemcallmatch+=" > matching_observed_tree_";
									finalsystemcallmatch+=chosentreenum+1;
									finalsystemcallmatch+=".tre";
									int returncodematch=system(finalsystemcallmatch.c_str());
									message="\t\tSaved matches to matching_observed_tree_";
									message+=chosentreenum+1;
									message+=".tre";
									PrintMessage();
								}
							}
							else {
								message="\nWarning: got return code of ";
								message+=returncode;
								message+="\n\nfor grep string of \n\n";
								message+=grepstring;
								message+="\n\n";
							
							}
							if (!debugmode) {
								system("rm tmp_mscount.txt");
								if (usefileregex) {
									system("rm tmp_grepregex.txt");
								}
							}
						}
				
					
					

				}
				break;
				
			
			
				}
			}
    }



    void BROWNIE::HandleExport( NexusToken& token)
    {
        nxsstring taxatoexport="";
        nxsstring internalstoexport="";
        double rate=1;
        double trend=0;
        nxsstring exportfname;
        int format=0; //0 means pagel
        bool replacing = false;
        bool nameprovided=true;
        bool hadfirstinternalnode=false;
        int source=0; //0 means use default character matrix, other sources involve simulation
        map<string,double> simcharmatrix;
        for(;;)
        {
            token.GetNextToken();

            if( token.Equals(";") ) {
                break;
            }
            else if( token.Abbreviation("FIle") ) {
                exportfname = GetFileName(token);
                nameprovided = true;
                //cout<<"Read file name"<<endl;
            }
            else if( token.Abbreviation("Replace") ) {
                nxsstring yesnoreplace=GetFileName(token);
                if (yesnoreplace[0] == 'n' || yesnoreplace[0] == 'N') {
                    replacing=false;
                }
                else {
                    replacing=true;
                }
            }
            else if( token.Abbreviation("Format") ) {
                nxsstring sourcestring=GetFileName(token);
                if (sourcestring[0] == 'p' || sourcestring[0] == 'P') {
                    format=0;
                }
            }

            else if( token.Abbreviation("Source") ) {
                nxsstring sourcestring=GetFileName(token);
                if (sourcestring[0] == 'c' || sourcestring[0] == 'C') {
                    source=0;
                }
                else if (sourcestring[0] == 't' || sourcestring[0] == 'T') {
                    source=1;
                }
                else {
                    //Simulate
                }
                // cout<<"Read source string"<<endl;
            }
            else if( token.Abbreviation("Trend") ) {
                nxsstring numbernexus = GetNumber(token);
                trend=atof( numbernexus.c_str() );
            }
            else if( token.Abbreviation("?") ) {
                message="Usage: Export [options...]\n\n";
                message+="Saves tree and data in a new format.\n\n";
                message+="Available options:\n\n";
                message+="Keyword ---- Option type ------------------------ Current setting --";
                message+="\nFile         <export-file-name>                   [None]";
                message+="\nReplace      yes | no                             false";
                //message+="\nFormat       Pagel | [other]                      Pagel";
                //message+="\nSource       CharMatrix | Trends                  CharMatrix";
                //message+="\nTrend        <double>                             ";
                //message+=trend;
                message+="\n\nCurrently, this only outputs the character matrix in Pagel format to a specified file.\nMake sure the tree has no polytomies.";
                PrintMessage();

            }
            else {
                errormsg = "Unexpected keyword (";
                errormsg += token.GetToken();
                errormsg += ") encountered reading PrintEdgeLengths command";
                throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        if (!nameprovided) {
            errormsg = "Must provide an output file name\n";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
        else {
            exportf.open( exportfname.c_str(), ios::out | ios::app );
            Tree t = intrees.GetIthTree(chosentree-1);
            //First do uppass to get Pagel like node numbers
            PreorderIterator <Node> m (t.GetRoot());
            Node *p = m.begin();
            int count=1+(continuouscharacters->GetNTax());
            while (p)
            {
                if(!(p->IsLeaf())) {
                    ostringstream NewLabelStream;
                    NewLabelStream<<count;
                    p->SetLabel(NewLabelStream.str());
                    count--;
                }
                p=m.next();
            }
            //Now do downpass to print out results
            // int count = 0;
            NodeIterator <Node> n (t.GetRoot());
            Node *q = n.begin();
            Node *anc = n.begin();
            // Node *a=n.begin();
            // Node *p=n.begin();
            //Node *b=n.begin();
            //Node *mrca=n.begin();
            //exportf<<"#PAG"<<endl;
            //cout<<"Successfully saved\n";
            if (source==0) {
                exportf<<continuouscharacters->GetNTax()<<" "<< continuouscharacters->GetNChar()<<"\n";
            }
            else if (source==1) {
                int ntax=continuouscharacters->GetNTax();
                exportf<<ntax;
                exportf<<" 1\n";
                simcharmatrix=SimulateBrownian(trend,1,0);
            }
            while (q->GetAnc())
            {
                string qlabel;
                qlabel=q->GetLabel();
                if (debugmode) {
                    cout<<qlabel<<endl;
                }
                // exportf << setiosflags(ios::right) << setw(4) << count					// arbitrary counter
                //      << " " << setw(8) << q->GetEdgeLength() 							// edge length
                //      << " " << setw(32) << setiosflags(ios::left) << q->GetLabel()	// node label
                //      ;
                //cout <<setw(32) << setiosflags(ios::left) << q->GetLabel()<<", ";
                if (q->IsLeaf()) {
                    taxatoexport+=(q->GetLabel()).c_str();
                    taxatoexport+=", ";
                    //exportf  << setiosflags(ios::left) << q->GetLabel()<<", ";
                }
                else if (q->GetAnc()) {
                    internalstoexport+=(q->GetLabel()).c_str();
                    internalstoexport+=", ";
                }
                if (q->GetAnc()) {
                    anc=q->GetAnc();
                    // if (strlen((anc->GetLabel()).c_str())==0) {
                    //     count++;
                    //     ostringstream NewLabelStream;
                    //     NewLabelStream<<count;
                    //     anc->SetLabel(NewLabelStream.str());
                    // }
                    if (q->IsLeaf()) {
                        taxatoexport+=(anc->GetLabel()).c_str();
                        taxatoexport+=", ";
                    }
                    else {
                        internalstoexport+=(anc->GetLabel()).c_str();
                        internalstoexport+=", ";
                    }
                    //exportf<<anc->GetLabel()<<", ";
                }
                else {
                    if (q->IsLeaf()) {
                        taxatoexport+=", ";
                    }
                    else if (q->GetAnc()) {
                        internalstoexport+=", ";
                    }
                    //exportf<<count<<", ";
                }
                double outputvalue=q->GetEdgeLength();
                if (q->IsLeaf()) {
                    taxatoexport+=outputvalue;
                }
                else if (q->GetAnc()) {
                    internalstoexport+=outputvalue;
                }
                //exportf << q->GetEdgeLength();
                if (q->IsLeaf())
                {
                    if (source==0) {
                        nxsstring currenttaxonlabel=(q->GetLabel()).c_str();
                        int taxonnumber=TaxonLabelToNumber( currenttaxonlabel );
                        for (int charnumber=0; charnumber<(continuouscharacters->GetNChar()); charnumber++) {
                            taxatoexport+=", ";
                            taxatoexport+=continuouscharacters->GetValue( taxonnumber, charnumber, true);
                        }
                    }
                    else if (source==1) {
                        map<string,double>::iterator pos;
                        pos=simcharmatrix.find(qlabel);
                        taxatoexport+=", ";
                        taxatoexport+=pos->second;
                        //exportf<<", "<<pos->second;
                    }

                }
                if (q->IsLeaf()) {
                    taxatoexport+="\n";
                }
                else if (q->GetAnc()) {
                    internalstoexport+="\n";
                }

                // exportf << endl;
                q = n.next();
            }
            exportf<<taxatoexport;
            exportf<<internalstoexport;
            exportf.close();
        }
    }

    nxsstring BROWNIE::PipeLeafGTP ()
    {
        int ntax=taxa->GetNumTaxonLabels();
        nxsstring output="";
        nxsstring TaxonLabel="";
        TaxonLabel+=cur->GetLabel();
        int TaxonNumber=taxa->FindTaxon(TaxonLabel);
        nxsstring NewLabel="taxon";
        if(TaxonNumber>=convertsamplestospecies.size()) {
            cout<<"Current taxon label is "<<TaxonLabel<<" which has number "<<TaxonNumber<<" but the convertsamplestospecies vector is of length "<<convertsamplestospecies.size()<<":"<<endl;
            for (int i=0;i<convertsamplestospecies.size();i++) {
                cout<<convertsamplestospecies[i]<<" ";
            }
            cout<<endl;
        }
        assert(TaxonNumber<convertsamplestospecies.size());
        NewLabel+=convertsamplestospecies[TaxonNumber];
        output=NEXUSString (NewLabel);
        return output;
    }

    nxsstring BROWNIE::PipeLeafSpeciesTree ()
    {
        nxsstring output="";
        output+=(cur->GetLabel());
		if (useCOAL || useMS)
		{
			output+=":";
			output+=cur->GetEdgeLength();
		}		
        return output;
    }

    nxsstring BROWNIE::PipeGTP (Tree t)
    {
        nxsstring TreeDescription="";
        cur = t.GetRoot();

        while (cur)
        {
            if (cur->GetChild())
            {
                TreeDescription+=PipeLeftParenthesis ();
                stk.push (cur);
                cur = cur->GetChild();
            }
            else
            {
                TreeDescription+=PipeLeafGTP ();
                while (!stk.empty() && (cur->GetSibling() == NULL))
                {
                    TreeDescription+=PipeRightParenthesis ();
                    cur = stk.top();
                    PipeInternal ();
                    stk.pop();
                }
                if (stk.empty())
                    cur = NULL;
                else
                {
                    TreeDescription+=PipeSiblingSymbol ();
                    cur = cur->GetSibling();
                }
            }
        }
        TreeDescription+=PipeEndOfTreeGTP ();
        return TreeDescription;
    }


    nxsstring BROWNIE::PipeSpeciesTree (ContainingTree *SpeciesTreePtr)
    {
        nxsstring TreeDescription="";
        cur = SpeciesTreePtr->GetRoot();

        while (cur)
        {
            if (cur->GetChild())
            {
                TreeDescription+=PipeLeftParenthesis ();
                stk.push (cur);
                cur = cur->GetChild();
            }
            else
            {
                TreeDescription+=PipeLeafSpeciesTree ();
                while (!stk.empty() && (cur->GetSibling() == NULL))
                {
                    TreeDescription+=PipeRightParenthesis ();
                    cur = stk.top();
                    TreeDescription+=PipeInternal ();
                    stk.pop();
                }
                if (stk.empty())
                    cur = NULL;
                else
                {
                    TreeDescription+=PipeSiblingSymbol ();
                    cur = cur->GetSibling();
                }
            }
        }
        TreeDescription+=PipeEndOfTreeSpeciesTree ();
        return TreeDescription;
    }



    nxsstring BROWNIE::PipeEndOfTreeGTP ()
    {
        nxsstring output="";
        if (outputallatonce) {
            output+=":";
            output+=trees->GetTreeWeight(chosentree);
        }
        else {
            output+=":1.00"; //Since we're piping one tree at a time, we don't need weights
        }

        output+=';';
        return output;
    }

    nxsstring BROWNIE::PipeEndOfTreeSpeciesTree ()
    {
        nxsstring output="";
        output+=';';
        return output;
    }


    //------------------------------------------------------------------------------
    nxsstring BROWNIE::PipeLeftParenthesis ()
    {
        nxsstring output="";
        output+='(';
        return output;
    }

    //------------------------------------------------------------------------------
    nxsstring BROWNIE::PipeRightParenthesis ()
    {
        nxsstring output="";
        output+=')';
        return output;
    }

    //------------------------------------------------------------------------------
    nxsstring BROWNIE::PipeSiblingSymbol ()
    {
        nxsstring output="";
        output+=',';
        return output;
    }


    //------------------------------------------------------------------------------
    nxsstring BROWNIE::PipeInternal ()
    {
        nxsstring output="";
        if (cur->GetLabel() != "")
            output+=cur->GetLabel();
		if (useCOAL || useMS)
		{
			output+=":";
			output+=cur->GetEdgeLength();
		}
		
        return output;
    }


    nxsstring BROWNIE::ReturnFinalSpeciesTree (Tree t)
    {
        nxsstring TreeDescription="";
        cur = t.GetRoot();

        while (cur)
        {
            if (cur->GetChild())
            {
                TreeDescription+=PipeLeftParenthesis ();
                stk.push (cur);
                cur = cur->GetChild();
            }
            else
            {
                TreeDescription+=PipeLeafSpeciesTree ();
                while (!stk.empty() && (cur->GetSibling() == NULL))
                {
                    TreeDescription+=PipeRightParenthesis ();
                    cur = stk.top();
                    PipeInternal ();
                    stk.pop();
                }
                if (stk.empty())
                    cur = NULL;
                else
                {
                    TreeDescription+=PipeSiblingSymbol ();
                    cur = cur->GetSibling();
                }
            }
        }
        TreeDescription+=PipeEndOfTreeSpeciesTree();
        return TreeDescription;
    }

    nxsstring BROWNIE::PipeLeafFinalSpeciesTree ()
    {
        int speciesnumber;
        nxsstring specieslabel=cur->GetLabel();
        string speciesstring=specieslabel.c_str();
        size_t index = speciesstring.find('t');
        speciesstring.erase(index,5); //erase "taxon"
                                      //cout<<"species string now "<<speciesstring<<endl;
        speciesnumber=atoi(speciesstring.c_str());
        nxsstring finallabel="(";
        int ntax=taxa->GetNumTaxonLabels();
        int timesprinted=0;
        for (int taxon=0; taxon<ntax; taxon++) {
            if (convertsamplestospecies[taxon]==speciesnumber) {
                if (timesprinted>0) {
                    finallabel+=",";
                }
                finallabel+=GetTaxonLabel(taxon);
                timesprinted++;
            }
        }
        finallabel+=")";
        // finallabel+=specieslabel;
        return finallabel;
    }


    map<string,double> BROWNIE::SimulateBrownian(double trend,double rate,double rootstate)
    {
        //variables
        map<string, double> newcharmatrix;
        message="Now simulating, trend = ";
        message+=trend;
        message+=" rate = ";
        message+=rate;
        message+=" state at root = ";
        message+=rootstate;
        PrintMessage();
        Tree t = intrees.GetIthTree(chosentree-1);
        PreorderIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        Node *ancestor=n.begin();
        while (q)
        {
            //cout <<setw(32) << setiosflags(ios::left) << q->GetLabel()<<", ";
            string qlabel;
            double newstate;
            qlabel=q->GetLabel();
            if (!q->GetAnc()) {
                ostringstream NewLabelStream;
                NewLabelStream<<rootstate;
                q->SetLabel(NewLabelStream.str());
            }
            else {
                //cout<<"Not doing root node\n";
                ancestor=q->GetAnc();
                double ancestralstate=atof((ancestor->GetLabel()).c_str());
                double brlen=q->GetEdgeLength();
                newstate=ancestralstate+(trend*brlen)+(gsl_ran_gaussian ( r, sqrt(rate*brlen)));
                ostringstream NewLabelStream;
                NewLabelStream<<newstate;
                q->SetLabel(NewLabelStream.str());
                if (debugmode) {
                    cout<<"Ancestral state = "<<ancestralstate<<", brlen = "<<brlen<<", newstate = "<<newstate<<endl;
                }
            }
            if (q->IsLeaf())
            {
                newcharmatrix[qlabel]=newstate;
            }
            //cout<<"Moving on to next node\n";
            q = n.next();
        }
        message="Done simulating";
        PrintMessage();
        return newcharmatrix;
    }
	
	double BROWNIE::browniesafe_gsl_sf_exp(double x) //Gets an exponential, but returns zero in case of underflow error
	{
		gsl_error_handler_t * old_handler =gsl_set_error_handler_off ();
		double result=gsl_sf_exp(x);
		if (gsl_isnan(result)) {
			result=0; //had some error, generally underflow
		}
		gsl_set_error_handler (old_handler);
		return result;
	}

    void BROWNIE::PreOrderTraversal(NexusToken& token)
    {
        //for(;;)
        // {
        //     token.GetNextToken();
        //      if( token.Equals(";") ) {
        //          break;
        //       }
        //   }
        Tree t = intrees.GetIthTree(chosentree-1);
        //t.SetInternalLabels(true);
        PreorderIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        Node *ancestor=n.begin();
        bool notfirstleaf=false;
        int count=0;
        while (q)
        {
            count++;
            cout << setiosflags(ios::right) << setw(4) << count					// arbitrary counter
                << " " << setw(8) << q->GetEdgeLength() 							// edge length
                << " " << setw(32) << setiosflags(ios::left) << q->GetLabel()	// node label
                ;
            if (q->IsLeaf())
            {
                cout << " [LEAF]";
                cout<<" length from root="<<q->GetHeight();
                notfirstleaf=true;
                ancestor=q;
                ancestor=ancestor->GetAnc();
                cout<<" label of ancestor = "<<ancestor->GetLabel();
            }
            else
            {
                cout << " [INTERNAL]";
                cout<<" length from root="<<q->GetHeight();
                char * NewLabel;
                sprintf(NewLabel,"%u",count);
                cout<<" did sprintf ";
                q->SetLabel(NewLabel);
                cout<<" new label = "<<q->GetLabel();
                cout<<" twice new label is "<<2*atoi((q->GetLabel()).c_str());
            }
            cout << endl;
            q = n.next();
        }
    }


    void BROWNIE::HandlePrintEdgeLengths(NexusToken& token)
    {
        for(;;)
        {
            token.GetNextToken();

            if( token.Equals(";") ) {
                break;
            }
            else {
                errormsg = "Unexpected keyword (";
                errormsg += token.GetToken();
                errormsg += ") encountered reading PrintEdgeLengths command";
                throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        Tree t = intrees.GetIthTree(chosentree-1);
		t.Update();
		t.SetPathLengths();
        cout << endl << "Postorder traversal of tree" << endl;
        cout << "Node Length     Label                           Type" << endl;
        cout << "----------------------------------------------------" << endl;
        int count = 0;
        NodeIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        Node *a=n.begin();
        Node *p=n.begin();
        Node *b=n.begin();
        Node *mrca=n.begin();
		nxsstring listofchanges="\n\nList of changes";
        bool notfirstleaf=false;
        while (q)
        {
            count++;
            cout << setiosflags(ios::right) << setw(4) << count					// arbitrary counter
                << " " << setw(8) << q->GetEdgeLength() 							// edge length
                << " " << setw(32) << setiosflags(ios::left) << q->GetLabel()	// node label
                ;
            if (q->IsLeaf())
            {
                cout << " [LEAF]";
                if (q->IsMarked()) {
                    cout <<" Node is marked\n\tModelCategorySize: ";
                }
                else {
                    cout <<" Node not marked\n\tModelCategorySize: ";
                }
                vector<double> modelcatoutput(q->GetModelCategory());
				cout<<"( ";
				for (int i=0;i<modelcatoutput.size();i++) {
					cout<<modelcatoutput[i]<<" ";
				}
				cout<<") "<<endl;
				vector<int> stateordervector(q->GetStateOrder()); 
				vector<double> statetimesvector(q->GetStateTimes());
				cout<<"\tStateOrder: ( ";
				for (int i=0;i<stateordervector.size();i++) {
					cout<<stateordervector[i]<<" ";
					if (i>0) {
						listofchanges+="\nchange ";
						listofchanges+=stateordervector[i-1];
						listofchanges+=" -> ";
						listofchanges+=stateordervector[i];
						double timeofchange=(q->GetPathLength())-(q->GetEdgeLength());
						for (int j=0; j<i; j++) {
							timeofchange+=statetimesvector[j];
						}
						listofchanges+=" at time above root of ";
						listofchanges+=timeofchange;
					}
				}
				cout<<") "<<endl;			
				cout<<"\tStateTimes: ( ";
				for (int i=0;i<statetimesvector.size();i++) {
					cout<<statetimesvector[i]<<" ";
				}
				cout<<") "<<endl;			
				
				
                //gsl_vector modelcatoutput(q->GetModelCategory());
                //cout<<modelcatoutput->size;
                cout<<"\t length from root="<<q->GetPathLength();
				cout<<endl<<" length of edge="<<q->GetEdgeLength();
                float pathlength=0;
                a=q;
                if(notfirstleaf) {
                    //  cout<<"Get MRCA of "<< q->GetLabel()<<" and "<<p->GetLabel()<<endl;
                    bool mrcanotfound=true;
                    while (a->GetAnc() && mrcanotfound) {
                        a=a->GetAnc();
                        b=p;
                        while (b->GetAnc() && mrcanotfound) {
                            b=b->GetAnc();
                            if (a==b) {
                                mrca=a;
                                mrcanotfound=false;
                                while(mrca->GetAnc()) {
                                    pathlength+=mrca->GetEdgeLength();
                                    //         cout<<"Edge length: "<<mrca->GetEdgeLength()<<" Total: "<<pathlength<<endl;
                                    mrca=mrca->GetAnc();
                                }
                            }
                        }
                    }
                }
                notfirstleaf=true;
                p=q;
                //  cout<<"Root to MRCA length="<<pathlength<<endl;
                //    while (aanc != Root && mrcanotfound) {
                //        banc=b;
                //        aanc=aanc->GetAnc();
                //       while (banc != Root && mrcanotfound) {
                //           banc=banc->GetAnc();
                //          if (aanc == banc) {
                //              mrcaptr=aanc;
                //               mrcanotfound=false;
                //           }
                //        }
                //    }
                //    return mrcaptr;


                //a=q;
                //cout<<endl<<"interim path length 1: "<<pathlength<<endl;
                //while (a->GetAnc()) {
                //   r=a->GetAnc();
                //  pathlength+=r->GetEdgeLength();
                // cout<<"interim: "<<pathlength<<" and edge "<<r->GetEdgeLength()<<endl;
                // a=r;
                //}
                //cout<<endl<<"Total length="<<pathlength<<endl;
            }
            else
            {
                cout << " [INTERNAL]";
                if (q->IsMarked()) {
                    cout <<" Node is marked  ModelCategorySize: ";
                }
                else {
                    cout <<" Node not marked  ModelCategorySize: ";
                }
				vector<double> modelcatoutput(q->GetModelCategory());
				cout<<"( ";
				for (int i=0;i<modelcatoutput.size();i++) {
					cout<<modelcatoutput[i]<<" ";
				}
				cout<<") "<<endl;
				vector<int> stateordervector(q->GetStateOrder()); 
				vector<double> statetimesvector(q->GetStateTimes());
				cout<<"\tStateOrder: ( ";
				for (int i=0;i<stateordervector.size();i++) {
					cout<<stateordervector[i]<<" ";
					if (i>0) {
						listofchanges+="\nchange ";
						listofchanges+=stateordervector[i-1];
						listofchanges+=" -> ";
						listofchanges+=stateordervector[i];
						double timeofchange=(q->GetPathLength())-(q->GetEdgeLength());
						for (int j=0; j<i; j++) {
							timeofchange+=statetimesvector[j];
						}
						listofchanges+=" at time above root of ";
						listofchanges+=timeofchange;
					}					
				}
				cout<<") "<<endl;			
				cout<<"\tStateTimes: ( ";
				for (int i=0;i<statetimesvector.size();i++) {
					cout<<statetimesvector[i]<<" ";
				}
				cout<<") "<<endl;	
                //gsl_vector *modelcatoutput;
                //modelcatoutput=gsl_vector_calloc(q->GetModelCategory());
                //gsl_vector modelcatoutput(q->GetModelCategory());
                //cout<<modelcatoutput->size;
                cout<<" length from root="<<q->GetPathLength();
				cout<<endl<<" length of edge="<<q->GetEdgeLength();
            }
            cout << endl;
            q = n.next();
        }
        cout << "----------------------------------------------------" << endl;
		message=listofchanges;
		PrintMessage();
    }

/**
* Copied from charactersblock.cpp
 * @method GetTaxonLabel [nxsstring:public]
 * @param i [int] the taxon's position in the taxa block
 *
 * Returns label for taxon number i (i ranges from 0 to ntax-1).
 */
nxsstring BROWNIE::GetTaxonLabel( int i )
{
    nxsstring s = taxa->GetTaxonLabel(i);
    return s;
}

//From TreeLib
nxsstring& BROWNIE::blanks_to_underscores( nxsstring& s )
{
    int len = s.length();
    for( int k = 0; k < len; k++ ) {
        if( s[k] == ' ' )
            s[k] = '_';
    }
    return s;
}

nxsstring& BROWNIE::underscores_to_blanks( nxsstring& s )
{
    int len = s.length();
    for( int k = 0; k < len; k++ ) {
        if( s[k] == '_' )
            s[k] = ' ';
    }
    return s;
}



void BROWNIE::HandleVCV(NexusToken& token)
{
    bool donenothing=true;
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") ) {
            if (donenothing) {
                message = "Usage: VCV taxset=<taxset-name>\n\n";
                PrintMessage();
            }
            break;
        }
        else if (token.Equals("?") ) {
            donenothing=false;
            message = "Usage: VCV taxset=<taxset-name>\n\n";
            PrintMessage();
        }
        else if (token.Abbreviation("Taxset") ) {
            donenothing=false;
            if (trees->GetNumTrees()<1) {
                errormsg = "Error: No valid trees are loaded yet\nYou can use Gettrees to load the trees.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            nxsstring currenttaxonlabel;
            nxsstring currenttaxonlabeltaxset;
std::string currenttaxonlabelnodeptr;
            nxsstring chosentaxset;
            int ntaxintaxset=0;
            chosentaxset=GetFileName(token);
            message="Taxset ";
            message+=chosentaxset.c_str();
            message+=": \n";
            IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
            if (taxonlist.empty()) {
                errormsg= "Error: Taxset ";
                errormsg+=chosentaxset.c_str();
                errormsg+=" does not exist.\nYou can define it using the taxset command.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
IntSet::const_iterator xi;
            for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
                ntaxintaxset++;
                currenttaxonlabel=taxa->GetTaxonLabel(*xi);
                message+="\t";
                message+=(*xi + 1);
                message+="=";
                message+=currenttaxonlabel.c_str();
                message+=" \n";
            }
            Tree t = intrees.GetIthTree(chosentree-1);
            gsl_matrix* VCVmatrix=gsl_matrix_calloc(ntaxintaxset,ntaxintaxset);
            VCVmatrix=GetVCV(chosentaxset);
			double 	TotalVCV=0.0;
            message+="With stem\n";
            for (int currentrow=0;currentrow<ntaxintaxset;currentrow++) {
                for (int currentcol=0;currentcol<ntaxintaxset;currentcol++) {
                    message+=gsl_matrix_get(VCVmatrix,currentrow,currentcol);
					TotalVCV+=gsl_matrix_get(VCVmatrix,currentrow,currentcol);
                    message+="\t";
                }
                message+="\n";
            }
            message+="\n";
			message+="\nTotal of all entries in VCV with stem: ";
			message+=TotalVCV;
			message+="\nAverage of all entries in VCV with stem: ";
			message+=TotalVCV/(ntaxintaxset*ntaxintaxset);
			PrintMessage();
			
            VCVmatrix=DeleteStem(GetVCV(chosentaxset));
			TotalVCV=0.0;
            message+="Without stem\n";
            for (int currentrow=0;currentrow<ntaxintaxset;currentrow++) {
                for (int currentcol=0;currentcol<ntaxintaxset;currentcol++) {
                    message+=gsl_matrix_get(VCVmatrix,currentrow,currentcol);
					TotalVCV+=gsl_matrix_get(VCVmatrix,currentrow,currentcol);
                    message+="\t";
                }
                message+="\n";
            }
			message+="\nTotal of all entries in VCV without stem: ";
			message+=TotalVCV;
			message+="\nAverage of all entries in VCV without stem: ";
			message+=TotalVCV/(ntaxintaxset*ntaxintaxset);
			PrintMessage();
			gsl_matrix_free(VCVmatrix);
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading VCV command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}

void BROWNIE::HandleTipValues(NexusToken& token)
{
    nxsstring chosentaxset;
    bool taxsetloaded=false;
    for(;;)
    {
        token.GetNextToken();

        if( token.Equals(";") && taxsetloaded!=true) {
            message="Error: Insufficient arguments. At least include taxset=taxset-name.\n\n";
            message += "Usage: TipValues taxset=taxset-name\n";
            PrintMessage();

            break;
        }
        else if (token.Equals(";") && taxsetloaded) {
            //Do tipvalue thing once everything else has been loaded
            int ntaxintaxset=0;
            IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
            if (taxonlist.empty()) {
                errormsg= "Error: Taxset ";
                errormsg+=chosentaxset.c_str();
                errormsg+=" does not exist.\nYou can define it using the taxset command.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
IntSet::const_iterator xi;
            for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
                ntaxintaxset++;
            }
            gsl_vector *tipvector;
            tipvector=gsl_vector_calloc(ntaxintaxset);
            //gsl_vector tipvector(ntaxintaxset,0);
            int rowcount=-1;
IntSet::const_iterator ri;
            for( ri = taxonlist.begin(); ri != taxonlist.end(); ri++ ) {
                rowcount++;
                gsl_vector_set(tipvector,rowcount,continuouscharacters->GetValue( *ri, chosenchar-1 ,true));
                message = "  ";
                message+=gsl_vector_get(tipvector,rowcount);
                message+="\n";
                PrintMessage();
            }
			gsl_vector_free(tipvector);
            break;
        }
        else if (token.Equals("?") ) {
            message = "Usage: TipValues taxset=taxset-name\nNote that char defaults to 1.";
            PrintMessage();
        }
        else if (token.Abbreviation("Taxset") ) {
            nxsstring currenttaxonlabel;
            nxsstring currenttaxonlabeltaxset;
            chosentaxset=GetFileName(token);
            message="Taxset " ;
            message+=chosentaxset.c_str();
            message+=", character ";
            message+=chosenchar;
            message+=":\n";
            PrintMessage();
            taxsetloaded=true;
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading TipValues command";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
    }
}

//Returns a vector containing tip values IN THE SAME ORDER AS THE TAXA IN THE TAXSET
gsl_vector* BROWNIE::GetTipValues(nxsstring chosentaxset, int chosenchar)
{
    int ntaxintaxset=0;
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg );
    }

IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntaxintaxset++;

    }
    gsl_vector *tipvalues;
    tipvalues=gsl_vector_calloc(ntaxintaxset);
    //gsl_vector tipvalues(ntaxintaxset,0);
    int rowcount=-1;
IntSet::const_iterator ri;
    for( ri = taxonlist.begin(); ri != taxonlist.end(); ri++ ) {
        rowcount++;
        gsl_vector_set(tipvalues,rowcount,continuouscharacters->GetValue( *ri, chosenchar-1, true));
    }

    return tipvalues;
}

//Returns a VCV matrix with columns and rows IN THE SAME ORDER AS THE TAXA IN THE TAXSET
gsl_matrix* BROWNIE::GetVCV(nxsstring chosentaxset)
{
    nxsstring currenttaxonlabel;
    nxsstring currenttaxonlabeltaxset;
	std::string currenttaxonlabelnodeptr;
    int ntaxintaxset=0;
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistr = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistc = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg );
    }
	IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntaxintaxset++;
        currenttaxonlabel=taxa->GetTaxonLabel(*xi);
    }
    Tree t = intrees.GetIthTree(chosentree-1);
    gsl_matrix *VCV=gsl_matrix_calloc(ntaxintaxset,ntaxintaxset);
    int rowcount=-1;
    nxsstring rtaxon;
    nxsstring ctaxon;
	IntSet::const_iterator ri;
	
    for( ri = taxonlistr.begin(); ri != taxonlistr.end(); ri++ ) {
        int colcount=-1;
        rowcount++;
        rtaxon=taxa->GetTaxonLabel(*ri);
        nxsstring rtaxonUnderscores=blanks_to_underscores(rtaxon);
        nxsstring rtaxonBlanks=underscores_to_blanks(rtaxon);
        Node *rleafptr;
        NodeIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        bool goodmatch=false;
        while (q)
        {
            if (q->IsLeaf())
            {
                string qlabel;
                qlabel=q->GetLabel();
                //nxsstring qlabel2<<qlabel;
                //qlabel2=blanks_to_underscores(qlabel2);
                if (qlabel==rtaxon.c_str() || qlabel==rtaxonUnderscores.c_str() || qlabel==rtaxonBlanks.c_str()) {
                    rleafptr=q;
                    goodmatch=true;
                }
            }
            q = n.next();
        }
        if (goodmatch==false) {
            errormsg= "Error: there was trouble identifying taxon ";
            errormsg+=rtaxon.c_str();
            errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
            throw XNexus (errormsg );
        }
		IntSet::const_iterator ci;
        for( ci = taxonlistc.begin(); ci != taxonlistc.end(); ci++ ) {
            colcount++;
            ctaxon=taxa->GetTaxonLabel(*ci);
            nxsstring ctaxonUnderscores=blanks_to_underscores(ctaxon);
            nxsstring ctaxonBlanks=underscores_to_blanks(ctaxon);
            Node *cleafptr;
            NodeIterator <Node> n (t.GetRoot());
            Node *q = n.begin();
            bool goodmatch=false;
            while (q)
            {
                if (q->IsLeaf())
                {
                    string qlabel;
                    qlabel=q->GetLabel();
					// nxsstring qlabel2<<qlabel;
     //qlabel2=blanks_to_underscores(qlabel2);
                    if (qlabel==ctaxon.c_str() || qlabel==ctaxonUnderscores.c_str() || qlabel==ctaxonBlanks.c_str()) {
                        cleafptr=q;
                        goodmatch=true;
                    }
                }
                q = n.next();
            }
            if (goodmatch==false) {
                errormsg= "Error: there was trouble identifying taxon ";
                errormsg+=ctaxon.c_str();
                errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
                throw XNexus (errormsg );
            }
            float pathlength=0;
            Node *a;
            a=rleafptr;
   /*         if (debugmode) {
                cout<<"Taxa "<<rleafptr->GetLabel()<<" (rtaxon="<<rtaxon<<") and "<<cleafptr->GetLabel()<<" (ctaxon="<<ctaxon<<"): ";
            } */
            if (cleafptr==rleafptr) { //We need the pendant edge length
                pathlength+=cleafptr->GetEdgeLength();
            }
            bool mrcanotfound=true;
 /*           if (debugmode) {
                cout<<pathlength;
            } */
            while (a->GetAnc() && mrcanotfound) { //General strategy here:
                a=a->GetAnc();
                Node *b;
                b=cleafptr;
                while (b->GetAnc() && mrcanotfound) {
                    b=b->GetAnc();
                    if (a==b) {
                        Node* mrca;
                        mrca=a;
                        mrcanotfound=false;
                        while(mrca->GetAnc()) {
                            pathlength+=mrca->GetEdgeLength();
/*                            if (debugmode) {
                                cout<<" + "<<mrca->GetEdgeLength();
                            } */
                            mrca=mrca->GetAnc();
                        }
                    }
                }
            }
            gsl_matrix_set(VCV,rowcount,colcount,pathlength);
  /*          if (debugmode) {
                cout<<" = "<<pathlength<<endl;
            } */
        }
    }
    //matrixsingular=TestSingularity(VCV);
    return VCV;
}

//Returns a VCV matrix with columns and rows IN THE SAME ORDER AS THE TAXA IN THE TAXSET, with edge lengths all raised to kappa power
gsl_matrix* BROWNIE::GetVCVwithKappa(nxsstring chosentaxset,double kappa)
{
    nxsstring currenttaxonlabel;
    nxsstring currenttaxonlabeltaxset;
	std::string currenttaxonlabelnodeptr;
    int ntaxintaxset=0;
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistr = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistc = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg );
    }
	IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntaxintaxset++;
        currenttaxonlabel=taxa->GetTaxonLabel(*xi);
    }
    Tree t = intrees.GetIthTree(chosentree-1);
    gsl_matrix *VCV=gsl_matrix_calloc(ntaxintaxset,ntaxintaxset);
    int rowcount=-1;
    nxsstring rtaxon;
    nxsstring ctaxon;
	IntSet::const_iterator ri;
	
    for( ri = taxonlistr.begin(); ri != taxonlistr.end(); ri++ ) {
        int colcount=-1;
        rowcount++;
        rtaxon=taxa->GetTaxonLabel(*ri);
        nxsstring rtaxonUnderscores=blanks_to_underscores(rtaxon);
        nxsstring rtaxonBlanks=underscores_to_blanks(rtaxon);
        Node *rleafptr;
        NodeIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        bool goodmatch=false;
        while (q)
        {
            if (q->IsLeaf())
            {
                string qlabel;
                qlabel=q->GetLabel();
                //nxsstring qlabel2<<qlabel;
                //qlabel2=blanks_to_underscores(qlabel2);
                if (qlabel==rtaxon.c_str() || qlabel==rtaxonUnderscores.c_str() || qlabel==rtaxonBlanks.c_str()) {
                    rleafptr=q;
                    goodmatch=true;
                }
            }
            q = n.next();
        }
        if (goodmatch==false) {
            errormsg= "Error: there was trouble identifying taxon ";
            errormsg+=rtaxon.c_str();
            errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
            throw XNexus (errormsg );
        }
		IntSet::const_iterator ci;
        for( ci = taxonlistc.begin(); ci != taxonlistc.end(); ci++ ) {
            colcount++;
            ctaxon=taxa->GetTaxonLabel(*ci);
            nxsstring ctaxonUnderscores=blanks_to_underscores(ctaxon);
            nxsstring ctaxonBlanks=underscores_to_blanks(ctaxon);
            Node *cleafptr;
            NodeIterator <Node> n (t.GetRoot());
            Node *q = n.begin();
            bool goodmatch=false;
            while (q)
            {
                if (q->IsLeaf())
                {
                    string qlabel;
                    qlabel=q->GetLabel();
					// nxsstring qlabel2<<qlabel;
     //qlabel2=blanks_to_underscores(qlabel2);
                    if (qlabel==ctaxon.c_str() || qlabel==ctaxonUnderscores.c_str() || qlabel==ctaxonBlanks.c_str()) {
                        cleafptr=q;
                        goodmatch=true;
                    }
                }
                q = n.next();
            }
            if (goodmatch==false) {
                errormsg= "Error: there was trouble identifying taxon ";
                errormsg+=ctaxon.c_str();
                errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
                throw XNexus (errormsg );
            }
            float pathlength=0;
            Node *a;
            a=rleafptr;
/*            if (debugmode) {
                cout<<"Taxa "<<rleafptr->GetLabel()<<" (rtaxon="<<rtaxon<<") and "<<cleafptr->GetLabel()<<" (ctaxon="<<ctaxon<<"): ";
            } */
            if (cleafptr==rleafptr) { //We need the pendant edge length
                pathlength+=pow(cleafptr->GetEdgeLength(),kappa);
            }
            bool mrcanotfound=true;
/*            if (debugmode) {
                cout<<pathlength;
            } */
            while (a->GetAnc() && mrcanotfound) { //General strategy here:
                a=a->GetAnc();
                Node *b;
                b=cleafptr;
                while (b->GetAnc() && mrcanotfound) {
                    b=b->GetAnc();
                    if (a==b) {
                        Node* mrca;
                        mrca=a;
                        mrcanotfound=false;
                        while(mrca->GetAnc()) {
                            pathlength+=pow(mrca->GetEdgeLength(),kappa);
/*                            if (debugmode) {
                                cout<<" + "<<pow(mrca->GetEdgeLength(),kappa);
                            } */
                            mrca=mrca->GetAnc();
                        }
                    }
                }
            }
            gsl_matrix_set(VCV,rowcount,colcount,pathlength);
/*            if (debugmode) {
                cout<<" = "<<pathlength<<endl;
            } */
        }
    }
    //matrixsingular=TestSingularity(VCV);
    return VCV;
}


//Returns a VCV matrix with columns and rows IN THE SAME ORDER AS THE TAXA IN THE TAXSET
gsl_matrix* BROWNIE::GetVCVwithTree(nxsstring chosentaxset,Tree t)
{
    //SimpleLCAQuery lca;
    //lca.SetTree (&t);
    nxsstring currenttaxonlabel;
    nxsstring currenttaxonlabeltaxset;
	std::string currenttaxonlabelnodeptr;
    int ntaxintaxset=0;
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistr = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistc = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg );
    }
	IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntaxintaxset++;
        currenttaxonlabel=taxa->GetTaxonLabel(*xi);
    }
    gsl_matrix *VCV=gsl_matrix_calloc(ntaxintaxset,ntaxintaxset);
    int rowcount=-1;
    nxsstring rtaxon;
    nxsstring ctaxon;
	IntSet::const_iterator ri;
	
    for( ri = taxonlistr.begin(); ri != taxonlistr.end(); ri++ ) {
        int colcount=-1;
        rowcount++;
        rtaxon=taxa->GetTaxonLabel(*ri);
        nxsstring rtaxonUnderscores=blanks_to_underscores(rtaxon);
        nxsstring rtaxonBlanks=underscores_to_blanks(rtaxon);
        Node *rleafptr;
        NodeIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        bool goodmatch=false;
        while (q)
        {
            if (q->IsLeaf())
            {
                string qlabel;
                qlabel=q->GetLabel();
                //nxsstring qlabel2<<qlabel;
                //qlabel2=blanks_to_underscores(qlabel2);
                if (qlabel==rtaxon.c_str() || qlabel==rtaxonUnderscores.c_str() || qlabel==rtaxonBlanks.c_str()) {
                    rleafptr=q;
                    goodmatch=true;
                }
            }
            q = n.next();
        }
        if (goodmatch==false) {
            errormsg= "Error: there was trouble identifying taxon ";
            errormsg+=rtaxon.c_str();
            errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
            throw XNexus (errormsg );
        }
		IntSet::const_iterator ci;
        for( ci = taxonlistc.begin(); ci != taxonlistc.end(); ci++ ) {
            colcount++;
            ctaxon=taxa->GetTaxonLabel(*ci);
            nxsstring ctaxonUnderscores=blanks_to_underscores(ctaxon);
            nxsstring ctaxonBlanks=underscores_to_blanks(ctaxon);
            Node *cleafptr;
            NodeIterator <Node> n (t.GetRoot());
            Node *q = n.begin();
            bool goodmatch=false;
            while (q)
            {
                if (q->IsLeaf())
                {
                    string qlabel;
                    qlabel=q->GetLabel();
                    // nxsstring qlabel2<<qlabel;
                    //qlabel2=blanks_to_underscores(qlabel2);
                    if (qlabel==ctaxon.c_str() || qlabel==ctaxonUnderscores.c_str() || qlabel==ctaxonBlanks.c_str()) {
                        cleafptr=q;
                        goodmatch=true;
                    }
                }
                q = n.next();
            }
            if (goodmatch==false) {
                errormsg= "Error: there was trouble identifying taxon ";
                errormsg+=ctaxon.c_str();
                errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
                throw XNexus (errormsg );
            }
            float pathlength=0;
            Node *a;
            a=rleafptr;
/*            if (debugmode) {
                cout<<"Taxa "<<rleafptr->GetLabel()<<" (rtaxon="<<rtaxon<<") and "<<cleafptr->GetLabel()<<" (ctaxon="<<ctaxon<<"): ";
            }*/
            if (cleafptr==rleafptr) { //We need the pendant edge length
                pathlength+=cleafptr->GetEdgeLength();
            }
            bool mrcanotfound=true;
/*            if (debugmode) {
                cout<<pathlength;
            }*/
            while (a->GetAnc() && mrcanotfound) { //General strategy here:
                a=a->GetAnc();
                Node *b;
                b=cleafptr;
                while (b->GetAnc() && mrcanotfound) {
                    b=b->GetAnc();
                    if (a==b) {
                        Node* mrca;
                        mrca=a;
                        mrcanotfound=false;
                        while(mrca->GetAnc()) {
                            pathlength+=mrca->GetEdgeLength();
/*                            if (debugmode) {
                                cout<<" + "<<mrca->GetEdgeLength();
                            }*/
                            mrca=mrca->GetAnc();
                        }
                    }
                }
            }
            gsl_matrix_set(VCV,rowcount,colcount,pathlength);
/*            if (debugmode) {
                cout<<" = "<<pathlength<<endl;
            }*/
        }
    }
    //matrixsingular=TestSingularity(VCV);
    return VCV;
}

//Returns for each taxon a table of start and stop times in the selected state. Assumes no more than maxstartstops/2 changes occur root to tip along tree per state
//Actually lists stop, then corresponding start, then next stop, then next corresponding start
gsl_matrix* BROWNIE::GetStartStopTimesforOneState(nxsstring chosentaxset, int selectedstate) {
	
    nxsstring currenttaxonlabel;
    nxsstring currenttaxonlabeltaxset;
	std::string currenttaxonlabelnodeptr;
    int ntaxintaxset=0;
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistr = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistc = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg );
    }
	IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntaxintaxset++;
        currenttaxonlabel=taxa->GetTaxonLabel(*xi);
    }
    Tree t = intrees.GetIthTree(chosentree-1);
    gsl_matrix * StartStopTimes=gsl_matrix_calloc(ntaxintaxset,maxstartstops*ntaxintaxset);
    gsl_matrix * VCV=gsl_matrix_calloc(ntaxintaxset,ntaxintaxset);
    VCV=DeleteStem(GetVCV(chosentaxset));
    double TotalTime=gsl_matrix_max (VCV);
    int rowcount=-1;
    nxsstring rtaxon;
	IntSet::const_iterator ri;
	
    for( ri = taxonlistr.begin(); ri != taxonlistr.end(); ri++ ) {
        int colcount=0;
        rowcount++;
        rtaxon=taxa->GetTaxonLabel(*ri);
        nxsstring rtaxonUnderscores=blanks_to_underscores(rtaxon);
        nxsstring rtaxonBlanks=underscores_to_blanks(rtaxon);
        Node *rleafptr;
        NodeIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        bool goodmatch=false;
        while (q)
        {
            if (q->IsLeaf())
            {
                string qlabel;
                qlabel=q->GetLabel();
                //nxsstring qlabel2<<qlabel;
                //qlabel2=blanks_to_underscores(qlabel2);
                if (qlabel==rtaxon.c_str() || qlabel==rtaxonUnderscores.c_str() || qlabel==rtaxonBlanks.c_str()) {
                    rleafptr=q;
                    goodmatch=true;
                }
            }
            q = n.next();
        }
        if (goodmatch==false) {
            errormsg= "Error: there was trouble identifying taxon ";
            errormsg+=rtaxon.c_str();
            errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
            throw XNexus (errormsg );
        }
        double pathlength=0;
        Node *a;
        a=rleafptr;
        while (a->GetAnc()) {
            vector<int> stateordervector(a->GetStateOrder());
            vector<double>  statetimesvector(a->GetStateTimes());
            double edgelength=a->GetEdgeLength();
            pathlength+=edgelength;
            double elapsedlength=0;
            for (int element=0;element<stateordervector.size();element++) {
                double segmentlength=statetimesvector[element];
              //  for (int position=0; position<staterestrictionvector.size(); position++) {
                  //  if (staterestrictionvector[position]==selectedstate) { //position is the called state, staterestrictionvector.at(position) gives the rate category to which that will be assigned
				if (staterestrictionvector[(stateordervector[element])]==selectedstate) {
					gsl_matrix_set(StartStopTimes,rowcount,colcount,fabs(TotalTime-pathlength+elapsedlength+segmentlength)); //do fabs because sometimes rounding error introduces -0.0...
					colcount++;
					gsl_matrix_set(StartStopTimes,rowcount,colcount,fabs(TotalTime-pathlength+elapsedlength)); 
					colcount++;
				}
				if (colcount>=StartStopTimes->size2) {
					errormsg= "Error: Too many appearances of state ";
					errormsg+=selectedstate;
					errormsg+=" on a path from the root to the tip of a tree\nUse \"model change=X\" and enter a larger number for changes than the current number (";
					errormsg+=maxstartstops/2;
					errormsg+=")";
					cout<<"StartStopTimes matrix\n";
					PrintMatrix(StartStopTimes);
					throw XNexus (errormsg );
				}
				elapsedlength+=segmentlength;
                   // }
               // }
            }
            a=a->GetAnc();
        }
    }
	gsl_matrix_free(VCV);
    return StartStopTimes;
}


/*Returns a VCV matrix for the selected taxset for one model
* with columns and rows IN THE SAME ORDER AS THE TAXA
* For each edge, there is a ModelCategory vector.
* For each ModelCategory (such as a morphological state), the corresponding
* entry in the vector is the amount of time spent in that Category on that edge.
*/
gsl_matrix* BROWNIE::GetVCVforOneModel(nxsstring chosentaxset, int selectedmodel)
{
    nxsstring currenttaxonlabel;
    nxsstring currenttaxonlabeltaxset;
std::string currenttaxonlabelnodeptr;
    int ntaxintaxset=0;
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistr = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistc = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg );
    }
IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntaxintaxset++;
        currenttaxonlabel=taxa->GetTaxonLabel(*xi);
    }
    Tree t = intrees.GetIthTree(chosentree-1);
    gsl_matrix *VCV=gsl_matrix_calloc(ntaxintaxset,ntaxintaxset);
    int rowcount=-1;
    nxsstring rtaxon;
    nxsstring ctaxon;
IntSet::const_iterator ri;

    for( ri = taxonlistr.begin(); ri != taxonlistr.end(); ri++ ) {
        int colcount=-1;
        rowcount++;
        rtaxon=taxa->GetTaxonLabel(*ri);
        nxsstring rtaxonUnderscores=blanks_to_underscores(rtaxon);
        nxsstring rtaxonBlanks=underscores_to_blanks(rtaxon);
        Node *rleafptr;
        NodeIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        bool goodmatch=false;
        while (q)
        {
            if (q->IsLeaf())
            {
                string qlabel;
                qlabel=q->GetLabel();
                //nxsstring qlabel2<<qlabel;
                //qlabel2=blanks_to_underscores(qlabel2);
                if (qlabel==rtaxon.c_str() || qlabel==rtaxonUnderscores.c_str() || qlabel==rtaxonBlanks.c_str()) {
                    rleafptr=q;
                    goodmatch=true;
                }
            }
            q = n.next();
        }
        if (goodmatch==false) {
            errormsg= "Error: there was trouble identifying taxon ";
            errormsg+=rtaxon.c_str();
            errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
            throw XNexus (errormsg );
        }
IntSet::const_iterator ci;
        for( ci = taxonlistc.begin(); ci != taxonlistc.end(); ci++ ) {
            colcount++;
            ctaxon=taxa->GetTaxonLabel(*ci);
            nxsstring ctaxonUnderscores=blanks_to_underscores(ctaxon);
            nxsstring ctaxonBlanks=underscores_to_blanks(ctaxon);
            Node *cleafptr;
            NodeIterator <Node> n (t.GetRoot());
            Node *q = n.begin();
            bool goodmatch=false;
            while (q)
            {
                if (q->IsLeaf())
                {
                    string qlabel;
                    qlabel=q->GetLabel();
                    // nxsstring qlabel2<<qlabel;
                    //qlabel2=blanks_to_underscores(qlabel2);
                    if (qlabel==ctaxon.c_str() || qlabel==ctaxonUnderscores.c_str() || qlabel==ctaxonBlanks.c_str()) {
                        cleafptr=q;
                        goodmatch=true;
                    }
                }
                q = n.next();
            }
            if (goodmatch==false) {
                errormsg= "Error: there was trouble identifying taxon ";
                errormsg+=ctaxon.c_str();
                errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
                throw XNexus (errormsg );
            }
            double pathlength=0;
            Node *a;
            a=rleafptr;
/*            if (debugmode) {
                cout<<"Taxa "<<rleafptr->GetLabel()<<" (rtaxon="<<rtaxon<<") and "<<cleafptr->GetLabel()<<" (ctaxon="<<ctaxon<<"): ";
            }*/
            if (cleafptr==rleafptr) { //We need the pendant edge length
                                      // float unweightedpathlength=cleafptr->GetEdgeLength();
                                      //gsl_vector *modelcategoryvector;
                                      //modelcategoryvector=gsl_vector_calloc(cleafptr->GetModelCategory());
                vector<double> modelcategoryvector(cleafptr->GetModelCategory());
                // if (modelcategoryvector->size==0) {
                //     pathlength+=0;
                //}
                //else if(modelcategoryvector->size<(selectedmodel-1)) {
                //			errormsg= "Error: ModelCategory vector too short for the number of submodels ";
                //			throw XNexus (errormsg );
                //		}
                //		else {
                for (int position=0;position<staterestrictionvector.size();position++) {
                    if (staterestrictionvector[position]==selectedmodel) { //position is the called state, staterestrictionvector.at(position) gives the rate category to which that state will be assigned
                        pathlength+=modelcategoryvector[position];
                    }
                }
                }
            bool mrcanotfound=true;
/*            if (debugmode) {
                cout<<pathlength;
            }*/
            while (a->GetAnc() && mrcanotfound) { //General strategy here:
                a=a->GetAnc();
                Node *b;
                b=cleafptr;
                while (b->GetAnc() && mrcanotfound) {
                    b=b->GetAnc();
                    if (a==b) {
                        Node* mrca;
                        mrca=a;
                        mrcanotfound=false;
                        while(mrca->GetAnc()) {
                            vector<double> modelcategoryvector(mrca->GetModelCategory());
							for (int position=0;position<staterestrictionvector.size();position++) {
								if (staterestrictionvector[position]==selectedmodel) { //position is the called state, staterestrictionvector.at(position) gives the rate category to which that state will be assigned
									pathlength+=modelcategoryvector[position];
								}
							}
                            mrca=mrca->GetAnc();
                        }
                    }
                }
            }
            gsl_matrix_set(VCV,rowcount,colcount,pathlength);
/*            if (debugmode) {
                cout<<" = "<<pathlength<<endl;
            }*/
            }
        }
//matrixsingular=TestSingularity(VCV);
return VCV;
    }


/*Returns a VCV matrix for the selected taxset for one model
* with columns and rows IN THE SAME ORDER AS THE TAXA
* For each edge, there is a ModelCategory vector.
* For each ModelCategory (such as a morphological state), the corresponding
* entry in the vector is the amount of time spent in that Category on that edge.
* This returns a VCV consisting only of branches with changes on them (wantchangeedges=true)
* Or a VCV consisting only of branches with no changes on them (wantchangeedges=false)
*/
gsl_matrix* BROWNIE::GetVCVforChangeNoChange(nxsstring chosentaxset, bool wantchangeedges)
{
    nxsstring currenttaxonlabel;
    nxsstring currenttaxonlabeltaxset;
std::string currenttaxonlabelnodeptr;
    int ntaxintaxset=0;
    IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistr = assumptions->GetTaxSet( chosentaxset );
    IntSet& taxonlistc = assumptions->GetTaxSet( chosentaxset );
    if (taxonlist.empty()) {
        errormsg= "Error: Taxset ";
        errormsg+=chosentaxset.c_str();
        errormsg+=" does not exist.\nYou can define it using the taxset command.";
        throw XNexus (errormsg );
    }
IntSet::const_iterator xi;
    for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
        ntaxintaxset++;
        currenttaxonlabel=taxa->GetTaxonLabel(*xi);
    }
    Tree t = intrees.GetIthTree(chosentree-1);
    gsl_matrix *VCV=gsl_matrix_calloc(ntaxintaxset,ntaxintaxset);
    int rowcount=-1;
    nxsstring rtaxon;
    nxsstring ctaxon;
IntSet::const_iterator ri;

    for( ri = taxonlistr.begin(); ri != taxonlistr.end(); ri++ ) {
        int colcount=-1;
        rowcount++;
        rtaxon=taxa->GetTaxonLabel(*ri);
        nxsstring rtaxonUnderscores=blanks_to_underscores(rtaxon);
        nxsstring rtaxonBlanks=underscores_to_blanks(rtaxon);
        Node *rleafptr;
        NodeIterator <Node> n (t.GetRoot());
        Node *q = n.begin();
        bool goodmatch=false;
        while (q)
        {
            if (q->IsLeaf())
            {
                string qlabel;
                qlabel=q->GetLabel();
                //nxsstring qlabel2<<qlabel;
                //qlabel2=blanks_to_underscores(qlabel2);
                if (qlabel==rtaxon.c_str() || qlabel==rtaxonUnderscores.c_str() || qlabel==rtaxonBlanks.c_str()) {
                    rleafptr=q;
                    goodmatch=true;
                }
            }
            q = n.next();
        }
        if (goodmatch==false) {
            errormsg= "Error: there was trouble identifying taxon ";
            errormsg+=rtaxon.c_str();
            errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
            throw XNexus (errormsg );
        }
IntSet::const_iterator ci;
        for( ci = taxonlistc.begin(); ci != taxonlistc.end(); ci++ ) {
            colcount++;
            ctaxon=taxa->GetTaxonLabel(*ci);
            nxsstring ctaxonUnderscores=blanks_to_underscores(ctaxon);
            nxsstring ctaxonBlanks=underscores_to_blanks(ctaxon);
            Node *cleafptr;
            NodeIterator <Node> n (t.GetRoot());
            Node *q = n.begin();
            bool goodmatch=false;
            while (q)
            {
                if (q->IsLeaf())
                {
                    string qlabel;
                    qlabel=q->GetLabel();
                    // nxsstring qlabel2<<qlabel;
                    //qlabel2=blanks_to_underscores(qlabel2);
                    if (qlabel==ctaxon.c_str() || qlabel==ctaxonUnderscores.c_str() || qlabel==ctaxonBlanks.c_str()) {
                        cleafptr=q;
                        goodmatch=true;
                    }
                }
                q = n.next();
            }
            if (goodmatch==false) {
                errormsg= "Error: there was trouble identifying taxon ";
                errormsg+=ctaxon.c_str();
                errormsg+=".\nTry removing strange characters (underscores, dashes,\nperiods, spaces, etc.) in its name. Sorry.\nPlease let me know about this error.";
                throw XNexus (errormsg );
            }
            double pathlength=0;
            Node *a;
            a=rleafptr;
/*            if (debugmode) {
                cout<<"Taxa "<<rleafptr->GetLabel()<<" (rtaxon="<<rtaxon<<") and "<<cleafptr->GetLabel()<<" (ctaxon="<<ctaxon<<"): ";
            }*/
            if (cleafptr==rleafptr) { //We need the pendant edge length
                                      // float unweightedpathlength=cleafptr->GetEdgeLength();
                                      //gsl_vector *modelcategoryvector;
                                      //modelcategoryvector=gsl_vector_calloc(cleafptr->GetModelCategory());
                vector<double> modelcategoryvector(cleafptr->GetModelCategory());
                // if (modelcategoryvector->size==0) {
                //     pathlength+=0;
                //}
                //else if(modelcategoryvector->size<(selectedmodel-1)) {
                //			errormsg= "Error: ModelCategory vector too short for the number of submodels ";
                //			throw XNexus (errormsg );
                //		}
                //		else {
                int numberofnonzeroentries=0;
                double temppathlength=0;
                for (int modelcat=0;modelcat<modelcategoryvector.size();modelcat++) {
                    if (modelcategoryvector[modelcat]>0) {
                        temppathlength+=modelcategoryvector[modelcat];
                        numberofnonzeroentries++;
                    }
                }
                if (wantchangeedges && (numberofnonzeroentries>1)) {
                    pathlength+=temppathlength;
                }
                if (!wantchangeedges && (numberofnonzeroentries<2)) {
                    pathlength+=temppathlength;
                }

                //		}
            }
            bool mrcanotfound=true;
/*            if (debugmode) {
                cout<<pathlength;
            }*/
            while (a->GetAnc() && mrcanotfound) { //General strategy here:
                a=a->GetAnc();
                Node *b;
                b=cleafptr;
                while (b->GetAnc() && mrcanotfound) {
                    b=b->GetAnc();
                    if (a==b) {
                        Node* mrca;
                        mrca=a;
                        mrcanotfound=false;
                        while(mrca->GetAnc()) {
                            vector<double> modelcategoryvector(mrca->GetModelCategory());
                            int numberofnonzeroentries=0;
                            double temppathlength=0;
                            for (int modelcat=0;modelcat<modelcategoryvector.size();modelcat++) {
                                if (modelcategoryvector[modelcat]>0) {
                                    temppathlength+=modelcategoryvector[modelcat];
                                    numberofnonzeroentries++;
                                }
                            }
                            if (wantchangeedges && (numberofnonzeroentries>1)) {
                                pathlength+=temppathlength;
/*                                if (debugmode) {
                                    cout<<" + "<<temppathlength;
                                }*/
                            }
                            if (!wantchangeedges && (numberofnonzeroentries<2)) {
                                pathlength+=temppathlength;
/*                                if (debugmode) {
                                    cout<<" + "<<temppathlength;
                                }*/
                            }
                            mrca=mrca->GetAnc();
                        }
                    }
                }
            }
            gsl_matrix_set(VCV,rowcount,colcount,pathlength);
/*            if (debugmode) {
                cout<<" = "<<pathlength<<endl;
            }*/
        }
    }
//matrixsingular=TestSingularity(VCV);
return VCV;
}

void BROWNIE::PrintMatrix(gsl_matrix *VCV)
{
	message="";
    for (int r=0;r<VCV->size1;r++) {
        for (int c=0;c<VCV->size2;c++) {
            message+=gsl_matrix_get(VCV,r,c);
			message+="\t";
        }
		message+="\n";
    }
    PrintMessage();
	
}

void BROWNIE::PrintVector(gsl_vector *somevector)
{
	message="( ";
    for (int r=0;r<somevector->size;r++) {
			message+=gsl_vector_get(somevector,r);
			message+=" ";
    }
	message+=")";
    PrintMessage();
	
}


//Use's Pagel Delta transform for continuous characters. First, run DeleteStem on input matrix
gsl_matrix* BROWNIE::ConvertVCVwithDelta(gsl_matrix *VCVorig,double delta)
{
    int ntax=VCVorig->size1;
    gsl_matrix* VCVfinal=gsl_matrix_calloc(ntax,ntax);
    for (int r=0;r<ntax;r++) {
        for (int c=0;c<ntax;c++) {
            gsl_matrix_set(VCVfinal,r,c,(pow(gsl_matrix_get(VCVorig,r,c),delta)));
        }
    }
	if (debugmode) {
		PrintMatrix(VCVorig);
		cout<<"\nAfter transformation with delta of "<<delta<<endl<<endl;
		PrintMatrix(VCVfinal);
	}
	//cout<<gsl_matrix_get(VCVorig,0,0)<<"^"<<delta<<" = (using pow) "<<pow(gsl_matrix_get(VCVorig,0,0),delta)<<" and from output "<<gsl_matrix_get(VCVfinal,0,0)<<endl;
    return VCVfinal;
}

//Use's Pagel Lambda (continuous char) transform for continuous characters. First, run DeleteStem on input matrix
gsl_matrix* BROWNIE::ConvertVCVwithLambda(gsl_matrix *VCVorig,double lambda)
{
    int ntax=VCVorig->size1;
    gsl_matrix* VCVfinal=gsl_matrix_calloc(ntax,ntax);
    for (int r=0;r<ntax;r++) {
        for (int c=0;c<ntax;c++) {
			if (r!=c) {
            gsl_matrix_set(VCVfinal,r,c,(lambda*(gsl_matrix_get(VCVorig,r,c))));
			}
        }
    }
    return VCVfinal;
}

gsl_matrix* BROWNIE::DeleteStem(gsl_matrix *VCVorig)
{
    int ntax=VCVorig->size1;
    double minvalue=gsl_matrix_min(VCVorig);
    gsl_matrix* VCVfinal=gsl_matrix_calloc(ntax,ntax);
    for (int r=0;r<ntax;r++) {
        for (int c=0;c<ntax;c++) {
            gsl_matrix_set(VCVfinal,r,c,((gsl_matrix_get(VCVorig,r,c))-minvalue));
        }
    }
    return VCVfinal;
}

//This takes a matrix (VCV times rate) and adds a vector of known tip variances to it (by, basically, multiplying the identity matrix by the tip variance vector and adding the resulting matrix element by element to the VCV times rate matrix)
gsl_matrix* BROWNIE::AddTipVarianceVectorToRateTimesVCV(gsl_matrix * VCVorig,gsl_vector * TipVariance)
{
    int ntax=VCVorig->size1;
    gsl_matrix* VCVfinal=gsl_matrix_calloc(ntax,ntax);
    gsl_matrix_memcpy (VCVfinal, VCVorig);
    for (int r=0;r<ntax;r++) {
        gsl_matrix_set(VCVfinal,r,r,(gsl_matrix_get(VCVfinal,r,r)+(gsl_vector_get(TipVariance,r))));
    }
    return VCVfinal;
}






void BROWNIE::HandleDebugOptimization( NexusToken& token )
{
	bool justdohelp=false;
    ofstream tablef;
    nxsstring tablefname;
    bool tablef_open=false;
    bool name_provided=false;
    nxsstring chosentaxset=globalchosentaxset;
    bool treeloop=false;
	bool charloop=false;
    bool adequateinput=true;
    int ntax=0;
    int nbest=1;
    int cstart=0;
    int cend=-1;
	bool appending=true;
	bool replacing=false;

    nxsstring treefilename="";
    bool definedfilename=false;
    nxsstring tmessage;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (adequateinput==false) {
                message="Insufficient input: type \"continuous ?\" for help";
                PrintMessage();
            }
            break;
        }
        else if (token.Abbreviation("TReeloop") ) {
            nxsstring yesnotreeloop=GetFileName(token);
            if (yesnotreeloop[0] == 'n') {
                treeloop=false;
            }
            else {
                treeloop=true;
            }
        }
        else if (token.Abbreviation("CHarloop") ) {
            nxsstring yesnotreeloop=GetFileName(token);
            if (yesnotreeloop[0] == 'n') {
                charloop=false;
            }
            else {
                charloop=true;
            }
        }
		else if( token.Abbreviation("Replace") ) {
            nxsstring yesnoreplace=GetFileName(token);
            if (yesnoreplace[0] == 'n') {
                replacing=false;
            }
            else {
                replacing=true;
            }
        }
        else if( token.Abbreviation("APpend") ) {
            nxsstring yesnoappend=GetFileName(token);
            if (yesnoappend[0] == 'n') {
                appending=false;
            }
            else {
                appending=true;
            }
        }		
        else if( token.Abbreviation("FIle") ) {
            tablefname = GetFileName(token);
            name_provided = true;			
        }
        else if( token.Abbreviation("Nbest") ) {
            
            nxsstring numbernexus = GetNumber(token);
            nbest=atoi( numbernexus.c_str() ); //convert to int
            message="You have chosen to save ";
            message+=nbest;
            message+=" trees";
            PrintMessage();
            if (nbest<1) {
                errormsg = "Error: must select a number greater than zero";
                nbest=1;
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
        }
        else if( token.Abbreviation("CStart") ) {
			
            nxsstring numbernexus = GetNumber(token);
            cstart=atoi( numbernexus.c_str() ); //convert to int
        }
        else if( token.Abbreviation("CEnd") ) {
			
            nxsstring numbernexus = GetNumber(token);
            cend=atoi( numbernexus.c_str() ); //convert to int
        }
        
		
		
        else if( token.Equals("?") ) {
            message="Usage: Continuous taxset=<chosen taxset> treeloop=[yes|no] append=[yes|no] replace=[yes|no]\n\n";
            message+="Returns the likelihood and AICc under the current model.\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --";
            message+="\nTaxset       <taxset name>                        *";
			message+=globalchosentaxset;
            message+="\nTreeloop     No|Yes                               *No";
            message+="\nCharloop     No|Yes                               *No";
            message+="\nFile         <file name>                          *None";
			message+="\nAppend       No|Yes                               *Yes";
			message+="\nReplace      No|Yes                               *No";

			//Next three commands are for a project with Justen
          //  message+="\nNBest        <integer>                            *1";
          //  message+="\nCStart       <integer>                            *none";
          //  message+="\nCEnd         <integer>                            *none";
            message+="\n                                                 *Option is nonpersistent\n\n";
            PrintMessage();
			justdohelp=true;
        }
        else if( token.Abbreviation("TReefile") ) {
            treefilename=GetFileName(token);
            definedfilename=true;
            break;
        }
        else if (token.Abbreviation("TAxset") ) {
            //adequateinput=true;
            if (trees->GetNumTrees()<1) {
                errormsg = "Error: No valid trees are loaded.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
            chosentaxset=GetFileName(token);
            IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
            if (taxonlist.empty()) {
                errormsg= "Error: Taxset ";
                errormsg+=chosentaxset.c_str();
                errormsg+=" does not exist.\nYou can define it using the taxset command.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
			IntSet::const_iterator xi;
            for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
                ntax++;
            }
        }
    }
    if (adequateinput && !justdohelp) {
		if( appending && replacing ) {
			errormsg = "Cannot specify APPEND and REPLACE at the same time";
			throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
		}		
		bool exists = FileExists( tablefname.c_str() );
		bool userok = true;
		if (appending && name_provided) {
			tablef_open = true;
			tablef.open( tablefname.c_str(), ios::out | ios::app );
			message = "\nAppending to continuous model output file (creating it if need be) ";
			message += tablefname;
			PrintMessage();
		}
		else if (name_provided) {
			if( exists && !replacing && !UserSaysOk( "Ok to replace?", "Continuous model output file specified already exists" ) )
				userok = false;
			if( userok && !tablef_open) {
				tablef_open = true;
				tablef.open( tablefname.c_str() );
			}
			if( exists && userok ) {
				message = "\nReplacing continuous model output file ";
				message += tablefname;
			}
			else if( userok ) {
				message = "\nContinuous model output file ";
				message += tablefname;
				message += " opened";
			}
			else {
				errormsg = "Aborting the continuous optimization so as not to overwrite the file.\n";
				throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
			}
			PrintMessage();		
		}

		message="Now optimizing with taxset ";
		message+=chosentaxset;
		PrintMessage();
        int starttree, stoptree, startchar, stopchar;
        int originalchosentree=chosentree;
        int originalchosenchar=chosenchar;
        if (treeloop) {
            starttree=1;
            stoptree=trees->GetNumTrees();
        }
        else {
            starttree=chosentree;
            stoptree=chosentree;
        }
		if (charloop) {
            startchar=1;
            stopchar=continuouscharacters->GetNChar();
        }
        else {
            startchar=chosenchar;
            stopchar=chosenchar;
        }
		if(ntax==0) {
			IntSet& taxonlist = assumptions->GetTaxSet( chosentaxset );
            if (taxonlist.empty()) {
                errormsg= "Error: Taxset ";
                errormsg+=chosentaxset.c_str();
                errormsg+=" does not exist.\nYou can define it using the taxset command.";
                throw XNexus (errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
            }
			IntSet::const_iterator xi;
            for( xi = taxonlist.begin(); xi != taxonlist.end(); xi++ ) {
                ntax++;
            }
			
		}
        for (chosentree=starttree;chosentree<=stoptree;chosentree++) {
			double treeweight=trees->GetTreeWeight(chosentree-1);
			nxsstring treename=trees->GetTreeName(chosentree-1);
			for (chosenchar=startchar;chosenchar<=stopchar;chosenchar++) {
				if (tablef_open) {
				//	tmessage="\n";
				//	tmessage+=chosentree;
				//	tmessage+="\t";
				//	tmessage+=chosenchar;
					//tablef<<tmessage;
					message="Now working on tree number ";
					message+=chosentree;
					message+=" char number ";
					message+=chosenchar;
					PrintMessage();
				}
				gsl_matrix * VCV=gsl_matrix_calloc(ntax,ntax);
				gsl_vector * tips=gsl_vector_calloc(ntax);
				gsl_vector * variance=gsl_vector_calloc(ntax);
				tips=GetTipValues(chosentaxset,chosenchar);
				if (tipvariancetype==2) {
					variance=GetTipValues(chosentaxset,chosenchar+1);
				}
				if (chosenmodel<5 || chosenmodel==21 || chosenmodel==22) {
					VCV=DeleteStem(GetVCV(chosentaxset));
					OptimizationFn my_fn(VCV,tips,variance,maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput);
					
					if (tipvariancetype==1 && chosenmodel==2) {
						gsl_vector *optimalrate=gsl_vector_calloc(5);
					//gsl_vector_memcpy(optimalrate,my_fn.OptimizeRateWithOptimizedTipVariance());
						gsl_vector_memcpy(optimalrate,my_fn.GeneralOptimization(2));
						optimalvaluescontinuouschar=gsl_vector_calloc(3);
						//optimalvalueslabels.clear();
						//optimalvalueslabels.push_back("rate");
						//optimalvalueslabels.push_back("tip variance");
						//optimalvalueslabels.push_back("lnL");
						gsl_vector_set(optimalvaluescontinuouschar,0,gsl_vector_get(optimalrate,0));
						gsl_vector_set(optimalvaluescontinuouschar,1,gsl_vector_get(optimalrate,1));
						gsl_vector_set(optimalvaluescontinuouschar,2,gsl_vector_get(optimalrate,4));
						message="\nOptimal rate = ";
						message+=gsl_vector_get(optimalrate,0);
						message+=" +/- ";
						message+=gsl_vector_get(optimalrate,2);
						message+="\nTip variance = ";
						message+=gsl_vector_get(optimalrate,1);
						message+=" +/- ";
						message+=gsl_vector_get(optimalrate,3);
						message+="\n-lnL = ";
						message+=gsl_vector_get(optimalrate,4);
						if (tablef_open) {
							tmessage="\t";
							tmessage+=gsl_vector_get(optimalrate,0);
							tmessage+="\t";
							tmessage+=gsl_vector_get(optimalrate,2);
							tmessage+="\t";
							tmessage+=gsl_vector_get(optimalrate,1);
							tmessage+="\t";
							tmessage+=gsl_vector_get(optimalrate,3);
							tmessage+="\t";
							tmessage+=gsl_vector_get(optimalrate,4);
							tablef<<tmessage;
						}
						gsl_vector_free(optimalrate);
					}
					else if (chosenmodel==1) {
						gsl_vector *optimalrate=gsl_vector_calloc(3);
					//gsl_vector_memcpy(optimalrate,my_fn.OptimizeRateWithGivenTipVariance());
						gsl_vector_memcpy(optimalrate,my_fn.GeneralOptimization(1));
						optimalvaluescontinuouschar=gsl_vector_calloc(2);
						//optimalvalueslabels.clear();
						//optimalvalueslabels.push_back("rate");
						//optimalvalueslabels.push_back("lnL");
						gsl_vector_set(optimalvaluescontinuouschar,0,gsl_vector_get(optimalrate,0));
						gsl_vector_set(optimalvaluescontinuouschar,1,gsl_vector_get(optimalrate,2));
						message+="\n-lnL = ";
						message+=gsl_vector_get(optimalrate,2);
						message+="\nAIC = ";
						message+=2*(1.0*gsl_vector_get(optimalrate,2) + 2);
						message+="\nAICc = ";
						message+=(2*1.0*gsl_vector_get(optimalrate,2))+4.0+12.0/(ntax-3);
						tips=GetTipValues(chosentaxset,chosenchar);
						message+="\nAncestral state = ";
						gsl_matrix *currentVCVmat=gsl_matrix_calloc(ntax,ntax);
						currentVCVmat=DeleteStem(GetVCV(chosentaxset));
						gsl_vector *tips=gsl_vector_calloc(ntax);
						tips=GetTipValues(chosentaxset,chosenchar);
						double ancstate=GetAncestralState(currentVCVmat,tips);	
						message+=ancstate;
						gsl_matrix_free(currentVCVmat);
						gsl_vector_free(tips);
						message+="\nOptimal rate = ";
						message+=gsl_vector_get(optimalrate,0);
						message+=" +/- ";
						message+=gsl_vector_get(optimalrate,1);
						if (tablef_open) {
							tmessage="Tree\tTree weight\tTree name\tChar\tModel\t-LnL\tAIC\tAICc\tAncState\tBMrate\n";
							tmessage+=chosentree;
							tmessage+="\t";
							tmessage+=treeweight;
							tmessage+="\t";
							tmessage+=treename;
							tmessage+="\t";
							tmessage+=chosenchar;
							tmessage+="\tBM1\t";
							tmessage+=gsl_vector_get(optimalrate,2);
							tmessage+="\t";
							tmessage+=2*(1.0*gsl_vector_get(optimalrate,2) + 2);
							tmessage+="\t";
							tmessage+=(2*1.0*gsl_vector_get(optimalrate,2))+2*2+2*2.0*(2+1)/(ntax-2-1);
							tmessage+="\t";
							tmessage+=ancstate;
							tmessage+="\t";
							tmessage+=gsl_vector_get(optimalrate,0);
							tmessage+="\n";
							tablef<<tmessage;
						}
						gsl_vector_free(optimalrate);
					}
					else if (chosenmodel==3) {
						gsl_vector *optimalrate=gsl_vector_calloc(7);
						gsl_vector_memcpy(optimalrate,my_fn.GeneralOptimization(3));
						optimalvaluescontinuouschar=gsl_vector_calloc(4);
						//optimalvalueslabels.clear();
						//optimalvalueslabels.push_back("rate");
						//optimalvalueslabels.push_back("ancestral state");
						//optimalvalueslabels.push_back("d");
						//optimalvalueslabels.push_back("lnL");
						gsl_vector_set(optimalvaluescontinuouschar,0,gsl_vector_get(optimalrate,0));
						gsl_vector_set(optimalvaluescontinuouschar,1,gsl_vector_get(optimalrate,1));
						gsl_vector_set(optimalvaluescontinuouschar,2,gsl_vector_get(optimalrate,2));
						gsl_vector_set(optimalvaluescontinuouschar,3,gsl_vector_get(optimalrate,6));
						message="\nOptimal rate = ";
						message+=gsl_vector_get(optimalrate,0);
						message+=" +/- ";
						message+=gsl_vector_get(optimalrate,3);
						message+="\nAncestral state = ";
						message+=gsl_vector_get(optimalrate,1);
						message+=" +/- ";
						message+=gsl_vector_get(optimalrate,4);
						message+="\nd = ";
						message+=gsl_vector_get(optimalrate,2);
						message+=" +/- ";
						message+=gsl_vector_get(optimalrate,5);
						message+="\n-lnL = ";
						message+=gsl_vector_get(optimalrate,6);
						gsl_vector_free(optimalrate);
					}
					else if (chosenmodel==4) {
						gsl_vector *optimalrate=gsl_vector_calloc(7);
						gsl_vector_memcpy(optimalrate,my_fn.GeneralOptimization(4));
						optimalvaluescontinuouschar=gsl_vector_calloc(4);
						//optimalvalueslabels.clear();
						//optimalvalueslabels.push_back("rate");
						//optimalvalueslabels.push_back("ancestral state");
						//optimalvalueslabels.push_back("g");
						//optimalvalueslabels.push_back("lnL");
						gsl_vector_set(optimalvaluescontinuouschar,0,gsl_vector_get(optimalrate,0));
						gsl_vector_set(optimalvaluescontinuouschar,1,gsl_vector_get(optimalrate,1));
						gsl_vector_set(optimalvaluescontinuouschar,2,gsl_vector_get(optimalrate,2));
						gsl_vector_set(optimalvaluescontinuouschar,3,gsl_vector_get(optimalrate,6));
						message="\nOptimal rate = ";
						message+=gsl_vector_get(optimalrate,0);
						message+=" +/- ";
						message+=gsl_vector_get(optimalrate,3);
						message+="\nAncestral state = ";
						message+=gsl_vector_get(optimalrate,1);
						message+=" +/- ";
						message+=gsl_vector_get(optimalrate,4);
						message+="\ng = ";
						message+=gsl_vector_get(optimalrate,2);
						message+=" +/- ";
						message+=gsl_vector_get(optimalrate,5);
						message+="\n-lnL = ";
						message+=gsl_vector_get(optimalrate,6);
						gsl_vector_free(optimalrate);
					}
					else if (chosenmodel==21) {
						gsl_vector *optimalrate=gsl_vector_calloc(7);
						gsl_vector_memcpy(optimalrate,my_fn.GeneralOptimization(21));
						optimalvaluescontinuouschar=gsl_vector_calloc(4);
						//optimalvalueslabels.clear();
						//optimalvalueslabels.push_back("rate");
						//optimalvalueslabels.push_back("ancestral state");
						//optimalvalueslabels.push_back("delta");
						//optimalvalueslabels.push_back("lnL");
						gsl_vector_set(optimalvaluescontinuouschar,0,gsl_vector_get(optimalrate,0));
						gsl_vector_set(optimalvaluescontinuouschar,1,gsl_vector_get(optimalrate,1));
						gsl_vector_set(optimalvaluescontinuouschar,2,gsl_vector_get(optimalrate,2));
						gsl_vector_set(optimalvaluescontinuouschar,3,gsl_vector_get(optimalrate,6));
						message="\nOptimal rate = ";
						message+=gsl_vector_get(optimalrate,0);
						message+=" +/- ";
						message+=gsl_vector_get(optimalrate,3);
						message+="\nAncestral state = ";
						message+=gsl_vector_get(optimalrate,1);
						message+=" +/- ";
						message+=gsl_vector_get(optimalrate,4);
						message+="\ndelta = ";
						message+=gsl_vector_get(optimalrate,2);
						message+=" +/- ";
						message+=gsl_vector_get(optimalrate,5);
						message+="\n-lnL = ";
						message+=gsl_vector_get(optimalrate,6);
						gsl_vector_free(optimalrate);
					}
					else if (chosenmodel==22) {
						gsl_vector *optimalrate=gsl_vector_calloc(7);
						gsl_vector_memcpy(optimalrate,my_fn.GeneralOptimization(22));
						optimalvaluescontinuouschar=gsl_vector_calloc(4);
						//optimalvalueslabels.clear();
						//optimalvalueslabels.push_back("rate");
						//optimalvalueslabels.push_back("ancestral state");
						//optimalvalueslabels.push_back("lambda");
						//optimalvalueslabels.push_back("lnL");
						gsl_vector_set(optimalvaluescontinuouschar,0,gsl_vector_get(optimalrate,0));
						gsl_vector_set(optimalvaluescontinuouschar,1,gsl_vector_get(optimalrate,1));
						gsl_vector_set(optimalvaluescontinuouschar,2,gsl_vector_get(optimalrate,2));
						gsl_vector_set(optimalvaluescontinuouschar,3,gsl_vector_get(optimalrate,6));
						message="\nOptimal rate = ";
						message+=gsl_vector_get(optimalrate,0);
						message+=" +/- ";
						message+=gsl_vector_get(optimalrate,3);
						message+="\nAncestral state = ";
						message+=gsl_vector_get(optimalrate,1);
						message+=" +/- ";
						message+=gsl_vector_get(optimalrate,4);
						message+="\nlambda = ";
						message+=gsl_vector_get(optimalrate,2);
						message+=" +/- ";
						message+=gsl_vector_get(optimalrate,5);
						message+="\n-lnL = ";
						message+=gsl_vector_get(optimalrate,6);
						gsl_vector_free(optimalrate);
					}
					message+="\n\nNote that +/- reflects imprecision due to numerical optimization";
					PrintMessage();
				}
				else { //we must have a model that deals with multiple VCVs
					gsl_matrix * VCV0=gsl_matrix_calloc(ntax,ntax);
					gsl_matrix * VCV1=gsl_matrix_calloc(ntax,ntax);
					gsl_matrix * VCV2=gsl_matrix_calloc(ntax,ntax);
					gsl_matrix * VCV3=gsl_matrix_calloc(ntax,ntax);
					gsl_matrix * VCV4=gsl_matrix_calloc(ntax,ntax);
					gsl_matrix * VCV5=gsl_matrix_calloc(ntax,ntax);
					gsl_matrix * VCV6=gsl_matrix_calloc(ntax,ntax);
					gsl_matrix * VCV7=gsl_matrix_calloc(ntax,ntax);
					gsl_matrix * VCV8=gsl_matrix_calloc(ntax,ntax);
					gsl_matrix * VCV9=gsl_matrix_calloc(ntax,ntax);
					if (chosenmodel==5) {
						VCV0=GetVCVforOneModel(chosentaxset,0);
                    //cout<<"VCV0 "<<gsl_matrix_get(VCV0,0,0);
						VCV1=GetVCVforOneModel(chosentaxset,1);
                   // cout<<"VCV1 "<<gsl_matrix_get(VCV1,0,0);
						VCV2=GetVCVforOneModel(chosentaxset,2);
                   // cout<<"VCV2 "<<gsl_matrix_get(VCV2,0,0);
						VCV3=GetVCVforOneModel(chosentaxset,3);
						VCV4=GetVCVforOneModel(chosentaxset,4);
						VCV5=GetVCVforOneModel(chosentaxset,5);
						VCV6=GetVCVforOneModel(chosentaxset,6);
						VCV7=GetVCVforOneModel(chosentaxset,7);
						VCV8=GetVCVforOneModel(chosentaxset,8);
						VCV9=GetVCVforOneModel(chosentaxset,9);
					}
					else if (chosenmodel==6) {
						VCV0=GetVCVforChangeNoChange(chosentaxset,false); //branches with no changes
						VCV1=GetVCVforChangeNoChange(chosentaxset,true); //branches with changes
					}
					if (debugmode) {
						cout<<"\n\nVCV0\n";
						for (int rt=0;rt<ntax;rt++) {
							cout<<"\n";
							for (int ct=0;ct<ntax;ct++) {
								cout<<gsl_matrix_get(VCV0,rt,ct)<<"\t";
							}
						}
						cout<<"\n\nVCV1\n";
						for (int rt=0;rt<ntax;rt++) {
							cout<<"\n";
							for (int ct=0;ct<ntax;ct++) {
								cout<<gsl_matrix_get(VCV1,rt,ct)<<"\t";
							}
						}
						cout<<"\n\nVCV2\n";
						for (int rt=0;rt<ntax;rt++) {
							cout<<"\n";
							for (int ct=0;ct<ntax;ct++) {
								cout<<gsl_matrix_get(VCV2,rt,ct)<<"\t";
							}
						}
						cout<<"\n\nVCV3\n";
						for (int rt=0;rt<ntax;rt++) {
							cout<<"\n";
							for (int ct=0;ct<ntax;ct++) {
								cout<<gsl_matrix_get(VCV3,rt,ct)<<"\t";
							}
						}
						cout<<"\n\nVCV4\n";
						for (int rt=0;rt<ntax;rt++) {
							cout<<"\n";
							for (int ct=0;ct<ntax;ct++) {
								cout<<gsl_matrix_get(VCV4,rt,ct)<<"\t";
							}
						}
						cout<<"\n\nVCV5\n";
						for (int rt=0;rt<ntax;rt++) {
							cout<<"\n";
							for (int ct=0;ct<ntax;ct++) {
								cout<<gsl_matrix_get(VCV5,rt,ct)<<"\t";
							}
						}
						cout<<"\n\nVCV6\n";
						for (int rt=0;rt<ntax;rt++) {
							cout<<"\n";
							for (int ct=0;ct<ntax;ct++) {
								cout<<gsl_matrix_get(VCV6,rt,ct)<<"\t";
							}
						}
						cout<<"\n\nVCV7\n";
						for (int rt=0;rt<ntax;rt++) {
							cout<<"\n";
							for (int ct=0;ct<ntax;ct++) {
								cout<<gsl_matrix_get(VCV7,rt,ct)<<"\t";
							}
						}
						cout<<"\n\nVCV8\n";
						for (int rt=0;rt<ntax;rt++) {
							cout<<"\n";
							for (int ct=0;ct<ntax;ct++) {
								cout<<gsl_matrix_get(VCV8,rt,ct)<<"\t";
							}
						}
						cout<<"\n\nVCV9\n";
						for (int rt=0;rt<ntax;rt++) {
							cout<<"\n";
							for (int ct=0;ct<ntax;ct++) {
								cout<<gsl_matrix_get(VCV9,rt,ct)<<"\t";
							}
						}
						
					}
					if (chosenmodel==5) {
						OptimizationFnMultiModel my_fn(VCV0,VCV1,VCV2,VCV3,VCV4,VCV5,VCV6,VCV7,VCV8,VCV9,tips,variance,maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput);
						gsl_vector *optimalvalues=gsl_vector_calloc(24);
						gsl_vector_memcpy(optimalvalues,my_fn.GeneralOptimization(5));
						int np=int(gsl_vector_get(optimalvalues,1));
						optimalvaluescontinuouschar=gsl_vector_calloc(np+1);
						//optimalvalueslabels.clear();
						//optimalvalueslabels.push_back("ancestral state");
						gsl_vector_set(optimalvaluescontinuouschar,0,gsl_vector_get(optimalvalues,2));
						for (int modelstate=0;modelstate<(np-1);modelstate++) {
							gsl_vector_set(optimalvaluescontinuouschar,modelstate+1,gsl_vector_get(optimalvalues,3+modelstate));
							nxsstring labelstring="rate";
							labelstring+=modelstate;
							//optimalvalueslabels.push_back(labelstring);
						}
						gsl_vector_set(optimalvaluescontinuouschar,np,gsl_vector_get(optimalvalues,0));
						//optimalvalueslabels.push_back("lnL");
						message="\n-lnL = ";
						message+=gsl_vector_get(optimalvalues,0);
						message+="\nAIC = ";
						message+=2*(1.0*gsl_vector_get(optimalvalues,0) + gsl_vector_get(optimalvalues,1));
						message+="\nAICc = ";
						message+=(2*1.0*gsl_vector_get(optimalvalues,0))+2*gsl_vector_get(optimalvalues,1)+2.0*gsl_vector_get(optimalvalues,1)*(gsl_vector_get(optimalvalues,1)+1)/(ntax-gsl_vector_get(optimalvalues,1)-1);
						message+="\nAncestral state = ";
						message+=gsl_vector_get(optimalvalues,2);
						message+=" +/- ";
						message+=gsl_vector_get(optimalvalues,13);
						tmessage="Tree\tTree weight\tTree name\tChar\tModel\t-LnL\tAIC\tAICc\tAncState";
						for (int modelstate=0;modelstate<(np-1);modelstate++) {
							message+="\nRate in state ";
							message+=modelstate;
							tmessage+="\tRate_in_state_";
							tmessage+=modelstate;
							message+=" = ";
							message+=gsl_vector_get(optimalvalues,3+modelstate);
							message+=" +/- ";
							message+=gsl_vector_get(optimalvalues,14+modelstate);
						}
						if (tablef_open) {
							tmessage+="\n";
							tmessage+=chosentree;
							tmessage+="\t";
							tmessage+=treeweight;
							tmessage+="\t";
							tmessage+=treename;
							tmessage+="\t";
							tmessage+=chosenchar;
							tmessage+="\tBMS\t";
							tmessage+=gsl_vector_get(optimalvalues,0);
							tmessage+="\t";
							tmessage+=2*(1.0*gsl_vector_get(optimalvalues,0) + gsl_vector_get(optimalvalues,1));
							tmessage+="\t";
							tmessage+=(2*1.0*gsl_vector_get(optimalvalues,0))+2*gsl_vector_get(optimalvalues,1)+2.0*gsl_vector_get(optimalvalues,1)*(gsl_vector_get(optimalvalues,1)+1)/(ntax-gsl_vector_get(optimalvalues,1)-1);
							tmessage+="\t";
							tmessage+=gsl_vector_get(optimalvalues,2);
							for (int modelstate=0;modelstate<(np-1);modelstate++) {
								tmessage+="\t";
								tmessage+=gsl_vector_get(optimalvalues,3+modelstate);
							}
							tmessage+="\n";
							tablef<<tmessage;
						}
						
						PrintMessage();
						gsl_vector_free(optimalvalues);
					}
					if (chosenmodel==12) {
						VCV0=DeleteStem(GetVCV(chosentaxset));
						gsl_matrix_free (VCV1);
						VCV1=gsl_matrix_calloc(ntax,maxstartstops*ntax); //assumes that states change fewer than three times per branch on pectinate tree
						VCV1=GetStartStopTimesforOneState(chosentaxset,0);
						gsl_matrix_free (VCV2);
						VCV2=gsl_matrix_calloc(ntax,maxstartstops*ntax);
						VCV2=GetStartStopTimesforOneState(chosentaxset,1);
						gsl_matrix_free (VCV3);
						VCV3=gsl_matrix_calloc(ntax,maxstartstops*ntax);
						VCV3=GetStartStopTimesforOneState(chosentaxset,2);
						gsl_matrix_free (VCV4);
						VCV4=gsl_matrix_calloc(ntax,maxstartstops*ntax);
						VCV4=GetStartStopTimesforOneState(chosentaxset,3);
						gsl_matrix_free (VCV5);
						VCV5=gsl_matrix_calloc(ntax,maxstartstops*ntax);
						VCV5=GetStartStopTimesforOneState(chosentaxset,4);
						gsl_matrix_free (VCV6);
						VCV6=gsl_matrix_calloc(ntax,maxstartstops*ntax);
						VCV6=GetStartStopTimesforOneState(chosentaxset,5);
						gsl_matrix_free (VCV7);
						VCV7=gsl_matrix_calloc(ntax,maxstartstops*ntax);
						VCV7=GetStartStopTimesforOneState(chosentaxset,6);
						gsl_matrix_free (VCV8);
						VCV8=gsl_matrix_calloc(ntax,maxstartstops*ntax);
						VCV8=GetStartStopTimesforOneState(chosentaxset,7);
						gsl_matrix_free (VCV9);
						VCV9=gsl_matrix_calloc(ntax,maxstartstops*ntax);
						VCV9=GetStartStopTimesforOneState(chosentaxset,8);
						OptimizationFnMultiModel my_fn(VCV0,VCV1,VCV2,VCV3,VCV4,VCV5,VCV6,VCV7,VCV8,VCV9,tips,variance,maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput);
						gsl_vector *optimalvalues=gsl_vector_calloc(28);
						gsl_vector_memcpy(optimalvalues,my_fn.GeneralOptimization(12));
						int np=int(gsl_vector_get(optimalvalues,1));
						optimalvaluescontinuouschar=gsl_vector_calloc(np+1);
						//optimalvalueslabels.clear();
						//optimalvalueslabels.push_back("rate");
						//optimalvalueslabels.push_back("attraction");
						//optimalvalueslabels.push_back("ancestral state");
						gsl_vector_set(optimalvaluescontinuouschar,0,gsl_vector_get(optimalvalues,2));
						gsl_vector_set(optimalvaluescontinuouschar,1,gsl_vector_get(optimalvalues,3));
						gsl_vector_set(optimalvaluescontinuouschar,2,gsl_vector_get(optimalvalues,4));
						for (int modelstate=0;modelstate<(np-3);modelstate++) {
							gsl_vector_set(optimalvaluescontinuouschar,modelstate+3,gsl_vector_get(optimalvalues,5+modelstate));
							nxsstring labelstring="optimum_state";
							labelstring+=modelstate;
							//optimalvalueslabels.push_back(labelstring);
						}
						gsl_vector_set(optimalvaluescontinuouschar,np,gsl_vector_get(optimalvalues,0));
						//optimalvalueslabels.push_back("lnL");
						message="\n-lnL = ";
						message+=gsl_vector_get(optimalvalues,0);
						message+="\nAIC = ";
						message+=2*(1.0*gsl_vector_get(optimalvalues,0) + gsl_vector_get(optimalvalues,1));
						message+="\nAICc = ";
						message+=(2*1.0*gsl_vector_get(optimalvalues,0))+2*gsl_vector_get(optimalvalues,1)+2.0*gsl_vector_get(optimalvalues,1)*(gsl_vector_get(optimalvalues,1)+1)/(ntax-gsl_vector_get(optimalvalues,1)-1);						
						message+="\nBM rate (sigma-squared) = ";
						message+=gsl_vector_get(optimalvalues,2);
                   // message+=" +/- ";
                    //message+=gsl_vector_get(optimalvalues,15);
						message+="\nOU attraction = ";
						message+=gsl_vector_get(optimalvalues,3);
						if (gsl_vector_get(optimalvalues,3)<=0.0015) {
							message+="  *Warning: program minimum value is 0.001! Try simple Brownian motion, since there is so little attraction.";
						}
						if (gsl_vector_get(optimalvalues,3)>=19.9) {
							message+="  *Warning: program maximum value is 20! Do you have variation in the terminals?";
						}
						
					//message+=" +/- ";
					//message+=gsl_vector_get(optimalvalues,16);
						message+="\nRoot state = ";
						message+=gsl_vector_get(optimalvalues,4);
					//message+=" +/- ";
					//message+=gsl_vector_get(optimalvalues,15);
						for (int modelstate=0;modelstate<(np-3);modelstate++) {
							message+="\nMean value in state ";
							message+=modelstate;
							message+=" = ";
							message+=gsl_vector_get(optimalvalues,5+modelstate);
                        //message+=" +/- ";
                        //message+=gsl_vector_get(optimalvalues,18+modelstate);
						}
						PrintMessage();
						if (tablef_open) {
							tmessage="Tree\tTree weight\tTree name\tChar\tModel\t-LnL\tAIC\tAICc\tAncState\tBMrate\tAttraction";
							for (int modelstate=0;modelstate<(np-3);modelstate++) {
								tmessage+="\tMean_in_state_";
								tmessage+=modelstate;
							}
							tmessage+="\n";
							tmessage+=chosentree;
							tmessage+="\t";
							tmessage+=treeweight;
							tmessage+="\t";
							tmessage+=treename;
							tmessage+="\t";
							tmessage+=chosenchar;
							if ((np-3)==1) {
								tmessage+="\tOU1\t";
							}
							else {
								tmessage+="\tOUSM\t";
							}
							tmessage+=gsl_vector_get(optimalvalues,0);
							tmessage+="\t";
							tmessage+=2*(1.0*gsl_vector_get(optimalvalues,0) + gsl_vector_get(optimalvalues,1));
							tmessage+="\t";
							tmessage+=(2*1.0*gsl_vector_get(optimalvalues,0))+2*gsl_vector_get(optimalvalues,1)+2.0*gsl_vector_get(optimalvalues,1)*(gsl_vector_get(optimalvalues,1)+1)/(ntax-gsl_vector_get(optimalvalues,1)-1);
							tmessage+="\t";
							tmessage+=gsl_vector_get(optimalvalues,4);
							tmessage+="\t";							
							tmessage+=gsl_vector_get(optimalvalues,2);
							tmessage+="\t";
							tmessage+=gsl_vector_get(optimalvalues,3);
							for (int modelstate=0;modelstate<(np-3);modelstate++) {
								tmessage+="\t";
								tmessage+=gsl_vector_get(optimalvalues,5+modelstate);
							}
							tmessage+="\n";
							tablef<<tmessage;
						}
						
						
						/*		gsl_matrix_free(VCV1);
						gsl_matrix_free(VCV2);
						gsl_matrix_free(VCV3);
						gsl_matrix_free(VCV4);
						gsl_matrix_free(VCV5);
						gsl_matrix_free(VCV6);
						gsl_matrix_free(VCV7);
						gsl_matrix_free(VCV8);
						gsl_matrix_free(VCV9);*/
						//cout<<"Vector of all output:\n";
						//PrintVector(optimalvalues);
						
						gsl_vector_free(optimalvalues);
					}
					if (chosenmodel==14) {
						Tree t = intrees.GetIthTree(chosentree-1);
						double MaxRootTipLength=t.GetMaxPathLength();
						if (debugmode) {
							cout<<"MaxRootTipLength = "<<MaxRootTipLength<<endl;
						}
					}
					if (chosenmodel==16 || chosenmodel==17) {
						if (chosenmodel==16) {
							treefilename="OutputTreesBMAN.nex";
						}
						else {
							treefilename="OutputTreesBMAO.nex";
						}
						ofstream outtreef;
						outtreef.open(treefilename.c_str());
						outtreef<<"#nexus\nbegin trees;";
						
						gsl_combination *c;
						size_t combinationsize;
                    //int ntax=(intrees.GetIthTree(0)).GetNumLeaves();
						int totalnumberofcombinations=0;
						for (combinationsize=0;combinationsize<(2*ntax-2);combinationsize++) {
							totalnumberofcombinations+=int(gsl_sf_choose ((2*ntax-2), combinationsize));
						}
						message="There are ";
						message+=totalnumberofcombinations;
						message+=" total combinations to try.";
						PrintMessage();
						ProgressBar(totalnumberofcombinations);
						int combinationcount=0;
						for (combinationsize=0;combinationsize<(2*ntax-2);combinationsize++) { //we want to try all possible assignments of zero, except assigning all branches to zero
							c=gsl_combination_calloc((2*ntax-2),combinationsize);
							do
							{
								combinationcount++;
								vector<int> nodestomakezero;
								for (int i=0;i<combinationsize;i++) {
									nodestomakezero.push_back(int(gsl_combination_get(c,i)));
								}
								nodestomakezero.push_back(-1); //just to keep things from going awry when we're done with all the nodes
								int vectorstep=0;
								int nodecount=0;
								Tree t=intrees.GetIthTree(chosentree-1);
								t.SetEdgeLengths(true);
								PreorderIterator <Node> n (t.GetRoot());
								NodePtr RootNode=t.GetRoot();
								NodePtr currentnode = n.begin();
								while (currentnode) {
									if (currentnode!=RootNode) {
										if (nodecount==nodestomakezero[vectorstep]) {
											currentnode->SetEdgeLength(0.0);
											vectorstep++;
										}
										else if (chosenmodel==17) {
											currentnode->SetEdgeLength(1.0);
										}
										nodecount++; //don't iterate on the root node, as it's not in the combination
									}
									currentnode=n.next();    
								}
								gsl_matrix *currentVCVmat=gsl_matrix_calloc(ntax,ntax);
								currentVCVmat=DeleteStem(GetVCVwithTree(chosentaxset,t));
								gsl_vector *tips=gsl_vector_calloc(ntax);
								gsl_vector *tipsresid=gsl_vector_calloc(ntax);
								double ancstate;
								double rate;
								double likelihood;
								tips=GetTipValues(chosentaxset,chosenchar);
								ancstate=GetAncestralState(currentVCVmat,tips);
								tipsresid=GetTipResiduals(tips,ancstate);
								rate=EstimateRate(currentVCVmat,tipsresid);
								gsl_matrix *RateTimesVCVfortest=gsl_matrix_calloc(currentVCVmat->size1,currentVCVmat->size2);
								gsl_matrix_memcpy(RateTimesVCVfortest, currentVCVmat);
								gsl_matrix_scale(RateTimesVCVfortest,rate);
								likelihood=GetLScore(currentVCVmat,tipsresid,rate);
								nxsstring rootlabel="-lnL_";
								rootlabel+=likelihood;
								RootNode->SetLabel(rootlabel);
								t.SetInternalLabels(true);
								t.Update();
                           // outtreef<<"gethasedgelengths = "<<t.GetHasEdgeLengths()<<endl;
								char outputstring[14];
								sprintf(outputstring,"%14.6f",likelihood);
                            //outtreef<<endl<<"gethasedgelengths = "<<t.GetHasEdgeLengths();
								
								outtreef<<"\ntree tree"<<combinationcount<<" = [-ln likelihood = "<<outputstring<<" ] ";
								NewickTreeWriter tw (&t);
								tw.SetStream (&outtreef);
								tw.SetWriteEdgeLengths(true);
								tw.Write();
								ProgressBar(0);
								gsl_matrix_free(currentVCVmat);
								gsl_vector_free(tips);
								gsl_vector_free(tipsresid);
								gsl_matrix_free(RateTimesVCVfortest);
								
							}
							while (gsl_combination_next(c) == GSL_SUCCESS);
							gsl_combination_free (c);
						}
						outtreef<<"\nend;";
						outtreef.close();
					}
					if (chosenmodel==18 || chosenmodel==19) {
						
                     //vector<Tree> BestTrees;
                     //vector<double> ScoresOfBestTrees;
						double bestscore=GSL_POSINF;
                     //cout<<"cstart = "<<cstart<<" cend = "<<cend<<" gslposinf = "<<int(GSL_POSINF)<<endl;
						if (chosenmodel==18 && !definedfilename) {
							treefilename="BMPN";
							if ((cstart>0) && (cend>0)) {
								treefilename+="_";
								treefilename+=cstart;
								treefilename+="_";
								treefilename+=cend;
								treefilename+=".tre";
							}
							else {
								treefilename+="_ALL.tre";
							}
						}
						else if (chosenmodel==19 && !definedfilename) {
							treefilename="BMPO";
							if ((cstart>0) && (cend>0)) {
								treefilename+="_";
								treefilename+=cstart;
								treefilename+="_";
								treefilename+=cend;
								treefilename+=".tre";
							}
							else {
								treefilename+="_ALL.tre";
							}
							
						}
						ofstream outtreef;
						outtreef.open(treefilename.c_str());
                   // if (outtreef.fail()) {
                   //     cout<<"Has trouble writing to file!!!!";
                 //   }
						
						outtreef<<"#nexus\nbegin trees;";
						outtreef.close();
						gsl_combination *c;
						size_t combinationsize;
                    //int ntax=(intrees.GetIthTree(0)).GetNumLeaves();
						int totalnumberofcombinations=0;
						for (combinationsize=0;combinationsize<=(ntax-2);combinationsize++) { //ntax-1 internal nodes, minus the root (reml problem)
							totalnumberofcombinations+=int(gsl_sf_choose ((ntax-2), combinationsize));
						}
						message="There are ";
						message+=totalnumberofcombinations;
						message+=" total combinations to try, you are doing ";
						if (cend<0) {
							cend=totalnumberofcombinations+1;
						}
                    //cout<<"cstart = "<<cstart<<" cend = "<<cend<<" gslposinf = "<<int(GSL_POSINF)<<endl;
						int actualcombinations=GSL_MIN(totalnumberofcombinations,cend-cstart+1);
						message+=actualcombinations;
						PrintMessage();
						int loopcounter=0;
						ProgressBar(actualcombinations);
						int combinationcount=0;
                    //int cutoff=int(floor(combinationcount/10.0));
						int cutoff=500;
						message="\nComb.\t-lnL";
						PrintMessage();
						for (combinationsize=0;combinationsize<=(ntax-2);combinationsize++) { //if a node is selected, its child has brlen of zero; otherwise, its child's sibling does
							c=gsl_combination_calloc((ntax-2),combinationsize);
							do
							{
								combinationcount++;
								if (combinationcount>=cstart && combinationcount<=cend) {
									loopcounter++;
									ProgressBar(0);
									vector<int> nodestomakezero;
									for (int i=0;i<combinationsize;i++) {
										nodestomakezero.push_back(int(gsl_combination_get(c,i)));
									}
									nodestomakezero.push_back(-1); //just to keep things from going awry when we're done with all the nodes
									int vectorstep=0;
									int nodecount=0;
									Tree t=intrees.GetIthTree(chosentree-1);
									t.SetEdgeLengths(true);
									t.SetInternalLabels(true);
									PreorderIterator <Node> n (t.GetRoot());
									NodePtr RootNode=t.GetRoot();
									NodePtr currentnode = n.begin();
									while (currentnode) {
										if (!(currentnode->IsLeaf()) && (currentnode!=RootNode)) {
											if (nodecount==nodestomakezero[vectorstep]) {
												(currentnode->GetChild())->SetEdgeLength(0.0);
												if (chosenmodel==19) {
													((currentnode->GetChild())->GetSibling())->SetEdgeLength(1.0);
												}
												vectorstep++;
											}
											else {
												((currentnode->GetChild())->GetSibling())->SetEdgeLength(0.0);
												if (chosenmodel==19) {
													(currentnode->GetChild())->SetEdgeLength(1.0);
												}
											}
											nodecount++;
										}
										else if((currentnode==RootNode) && (chosenmodel==19)) {
											(currentnode->GetChild())->SetEdgeLength(1.0);
											((currentnode->GetChild())->GetSibling())->SetEdgeLength(1.0);
										}
										currentnode=n.next();
									}
									gsl_matrix *currentVCVmat=gsl_matrix_calloc(ntax,ntax);
									currentVCVmat=DeleteStem(GetVCVwithTree(chosentaxset,t));
									gsl_vector *tips=gsl_vector_calloc(ntax);
									gsl_vector *tipsresid=gsl_vector_calloc(ntax);
									double ancstate;
									double rate;
									double likelihood;
									tips=GetTipValues(chosentaxset,chosenchar);
									ancstate=GetAncestralState(currentVCVmat,tips);
									tipsresid=GetTipResiduals(tips,ancstate);
									rate=EstimateRate(currentVCVmat,tipsresid);
									gsl_matrix *RateTimesVCVfortest=gsl_matrix_calloc(currentVCVmat->size1,currentVCVmat->size2);
									gsl_matrix_memcpy(RateTimesVCVfortest, currentVCVmat);
									gsl_matrix_scale(RateTimesVCVfortest,rate);
									likelihood=GetLScore(currentVCVmat,tipsresid,rate);
                                //t.Draw(cout);
                                //cout<<"likelihood is "<<likelihood<<endl;
									if (likelihood==likelihood) { //test for NaN
                                                              // if (likelihood<=bestscore) {
                                                              // if(BestTrees.size()<nbest) {
                                                              //       BestTrees.push_back(t);
                                                              //      ScoresOfBestTrees.push_back(likelihood);
                                                              //  }
                                                              //   else {
                                                              //       for (int i=0;i<nbest-1;i++) {
                                                              //           BestTrees[i]=BestTrees[i+1];
                                                              //           ScoresOfBestTrees[i]=ScoresOfBestTrees[i+1]; //move everything over one, except for the last cell, which we'll replace
                                                              //       }
                                                              //       BestTrees.pop_back();
                                                              //       BestTrees.push_back(t);
                                                              //       ScoresOfBestTrees.pop_back();
                                                              //       ScoresOfBestTrees.push_back(likelihood);
                                                              //   }
                                                              // if (likelihood<bestscore) { //overwrite
                                                              //      outtreef.open(treefilename.c_str());
                                                              //      outtreef<<"#nexus\nbegin trees;";
                                                              // }
                                                              // else {
                                                              //save all trees
										outtreef.open(treefilename.c_str(), ios::out | ios::app );
                                    // }
										if (likelihood<=bestscore) {
											bestscore=likelihood;
										}
										nxsstring rootlabel="-lnL_";
										rootlabel+=likelihood;
										RootNode->SetLabel(rootlabel);
										outtreef<<"\ntree tree"<<combinationcount<<" = [&R] [-ln likelihood = "<<likelihood<<" ancstate = "<<ancstate<<" rate = "<<rate<<" ] ";
										t.Write(outtreef);
										outtreef.close();
										
															  }
									if (loopcounter==cutoff) {
										message="";
										message+=combinationcount;
										message+="\t";
										message+=bestscore;
										PrintMessage();
										loopcounter=0;
									}
									gsl_matrix_free(currentVCVmat);
									gsl_vector_free(tips);
									gsl_vector_free(tipsresid);
									gsl_matrix_free(RateTimesVCVfortest);
									}
								}
							while ((gsl_combination_next(c) == GSL_SUCCESS) && (combinationcount<=cend));
							gsl_combination_free (c);
							}
                   // cout<<"Saving "<<BestTrees.size()<<" best trees"<<endl;
                    //for (int i=0;i<BestTrees.size();i++) {
                    //    Tree goodtree=BestTrees.back();
                     //   goodtree.Update();
                     //   goodtree.Draw(cout);
                     //   BestTrees.pop_back();
                   //     double goodlikelihood=ScoresOfBestTrees.back();
                   //     ScoresOfBestTrees.pop_back();
                   //     cout<<"lnL = "<<goodlikelihood<<endl;
                   // outtreef<<"\ntree tree"<<i+1<<" = [&R] [likelihood = "<<goodlikelihood<<" ] ";
                   //  goodtree.Write(outtreef);
                   // }
						outtreef.open(treefilename.c_str(), ios::out | ios::app );
						outtreef<<"\nend;";
						outtreef.close();
						}
					
					if (chosenmodel==6) {
						OptimizationFnMultiModel my_fn(VCV0,VCV1,VCV2,VCV3,VCV4,VCV5,VCV6,VCV7,VCV8,VCV9,tips,variance,maxiterations, stoppingprecision, randomstarts, stepsize,detailedoutput);
						gsl_vector *optimalvalues=gsl_vector_calloc(7);
						gsl_vector_memcpy(optimalvalues,my_fn.GeneralOptimization(6));
						optimalvaluescontinuouschar=gsl_vector_calloc(4);
						//optimalvalueslabels.clear();
						//optimalvalueslabels.push_back("ancestral state");
						//optimalvalueslabels.push_back("rate on branches with no changes");
						//optimalvalueslabels.push_back("rate on branches with changes");
						//optimalvalueslabels.push_back("lnL");
						gsl_vector_set(optimalvaluescontinuouschar,0,gsl_vector_get(optimalvalues,1));
						gsl_vector_set(optimalvaluescontinuouschar,1,gsl_vector_get(optimalvalues,2));
						gsl_vector_set(optimalvaluescontinuouschar,2,gsl_vector_get(optimalvalues,3));
						gsl_vector_set(optimalvaluescontinuouschar,3,gsl_vector_get(optimalvalues,0));
						message="\n-lnL = ";
						message+=gsl_vector_get(optimalvalues,0);
						message+="\nAIC = ";
						message+=2*(1.0*gsl_vector_get(optimalvalues,0) + 3);
						message+="\nAICc = ";
						message+=(2*1.0*gsl_vector_get(optimalvalues,0))+2*3+2.0*3*(4.0)/(ntax-3-1);
						message+="\nAncestral state = ";
						message+=gsl_vector_get(optimalvalues,1);
						message+=" +/- ";
						message+=gsl_vector_get(optimalvalues,4);
						message+="\nRate on branches with no changes = ";
						message+=gsl_vector_get(optimalvalues,2);
						message+=" +/- ";
						message+=gsl_vector_get(optimalvalues,5);
						message+="\nRate on branches with changes = ";
						message+=gsl_vector_get(optimalvalues,3);
						message+=" +/- ";
						message+=gsl_vector_get(optimalvalues,6);
						PrintMessage();
						if (tablef_open) {
							tmessage="Tree\tTree weight\tTree name\tChar\tModel\t-LnL\tAIC\tAICc\tAncState\tRate_No_changes\tRate_Changes\n";
							tmessage+=chosentree;
							tmessage+="\t";
							tmessage+=treeweight;
							tmessage+="\t";
							tmessage+=treename;
							tmessage+="\t";
							tmessage+=chosenchar;
							tmessage+="\tBMC\t";
							tmessage+=gsl_vector_get(optimalvalues,0);
							tmessage+="\t";
							tmessage+=2*(1.0*gsl_vector_get(optimalvalues,0) + 3);
							tmessage+="\t";
							tmessage+=(2*1.0*gsl_vector_get(optimalvalues,0))+2*3+2.0*3*(4.0)/(ntax-3-1);
							tmessage+="\t";
							tmessage+=gsl_vector_get(optimalvalues,1);
							tmessage+="\t";							
							tmessage+=gsl_vector_get(optimalvalues,2);
							tmessage+="\t";
							tmessage+=gsl_vector_get(optimalvalues,3);
							tmessage+="\n";
							tablef<<tmessage;
						}
						gsl_vector_free(optimalvalues);
					}
					gsl_matrix_free(VCV0);
					gsl_matrix_free(VCV1);
					gsl_matrix_free(VCV2);
					gsl_matrix_free(VCV3);
					gsl_matrix_free(VCV4);
					gsl_matrix_free(VCV5);
					gsl_matrix_free(VCV6);
					gsl_matrix_free(VCV7);
					gsl_matrix_free(VCV8);
					gsl_matrix_free(VCV9);
					
					}
				gsl_matrix_free(VCV);
				gsl_vector_free(tips);
				gsl_vector_free(variance);
				}
			}
        chosentree=originalchosentree; //restore initial values.
        chosenchar=originalchosenchar;
        if (tablef_open) {
            tablef.close();
        }
		}
	}


//Returns the ancestral state, using formula from Martins and Lamont 1998, Animal Behavior 55: 1685-1706, bottom right of page 1689.
double BROWNIE::GetAncestralState(gsl_matrix *VCV, gsl_vector *tips)
{
    if (debugmode) {
        cout<<"\nNow in GetAncestralState\n";
    }
    int ntax=tips->size;
    if (debugmode) {
        message="\nNtax = ";
        message+=ntax;
        message+="\n";
        PrintMessage();
    }
    //cout<<"ntax="<<ntax<<endl;
    double ancestralstate;
    //matrixsingular=TestSingularity(VCV);
    // if (matrixsingular) {
    //     errormsg="Error: VCV matrix is singular in GetAncestralState ";
    //    throw XNexus( errormsg);
    // }
    if (debugmode) {
        message="\nNow calculating inverse VCV for estimating ancestral state.\n";
        message+="\n\nntax is ";
        message+=ntax;
        PrintMessage();
    }
    if(debugmode) {
        message="debugging line 3494\n";
        for (int currentrow=0;currentrow<VCV->size1;currentrow++) {
            for (int currentcol=0;currentcol<VCV->size2;currentcol++) {
                message+=gsl_matrix_get(VCV,currentrow,currentcol);
                message+="\t";
            }
            message+="\n";
        }
        PrintMessage();
    }

    gsl_permutation * p = gsl_permutation_alloc (ntax);
    int signum;
    gsl_matrix * VCVinverse=gsl_matrix_calloc(ntax,ntax);
    gsl_matrix * VCVinversestart=gsl_matrix_calloc(ntax,ntax);
    gsl_matrix_memcpy (VCVinversestart, VCV);
    if(debugmode) {
        message="debugging line 3511\n";
        for (int currentrow=0;currentrow<VCV->size1;currentrow++) {
            for (int currentcol=0;currentcol<VCV->size2;currentcol++) {
                message+=gsl_matrix_get(VCV,currentrow,currentcol);
                message+="\t";
            }
            message+="\n";
        }
        PrintMessage();
    }
    gsl_linalg_LU_decomp (VCVinversestart,p, &signum);
    gsl_linalg_LU_invert (VCVinversestart,p, VCVinverse);
    //VCVinverse=Inverse(VCV);
    if (debugmode) {
        message="\nFinished calculating inverse VCV for estimating ancestral state.\n";
        PrintMessage();
    }
    gsl_matrix * tipsasmatrix=gsl_matrix_calloc(ntax,1);
    for (int i=0; i<ntax; i++) {
        gsl_matrix_set(tipsasmatrix,i,0,gsl_vector_get(tips,i));
    }
    gsl_matrix * stepA=gsl_matrix_calloc(ntax,1);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, VCVinverse, tipsasmatrix,0, stepA);
    //stepA=VCVinverse*tipsasmatrix;
    double stepB=0.0;
    for (int i=0; i<ntax; i++) {
        stepB+=gsl_matrix_get(stepA,i,0);
        //cout<<stepA[i][0]<<" stepB "<<stepB<<endl;

    }
    double stepC=0.0;
    for (int i=0; i<ntax; i++) {
        for (int j=0; j<ntax; j++) {
            stepC+=gsl_matrix_get(VCVinverse,i,j);
            // cout<<"StepC  "<<stepC<<endl;
        }
    }
    if (stepC==0) {
        errormsg="Error: Division by zero in GetAncestralState routine";
        throw XNexus(errormsg);
    }
    else {
        ancestralstate=stepB/stepC;
    }

    //gsl_matrix J(ntax,1,1);
    //gsl_matrix Jprime(1,ntax,1);
    //gsl_matrix step1(ntax,1,0);
    //   step1=VCVinverse*J;
    //  gsl_vector step1vect(ntax,0);
    //  for (int i=0; i<ntax; i++) {
    //      step1vect[i]=step1[i][0];
    //  }
    //cout<<"VCVinverse*J"<<endl<<step1<<endl;
    //   gsl_vector step2vect(1,0);
    //   step2vect=MatrixTimesVector(Jprime,step1vect);
    //   gsl_matrix step3(1,1,0);
    //   gsl_matrix step2(1,1,0);
    //   step2[0][0]=step2vect[0];
    //cout<<step2;
    //   matrixsingular=TestSingularity(step2);
    //   if (matrixsingular) {
    //       errormsg="Singular matrix (step2) during get ancestral state";
    //      throw XNexus(errormsg);
    //  }
    //  step3=Inverse(step2);
    //   gsl_matrix tipsasmatrix(ntax,1,0);
    //   gsl_matrix tipsasmatrixprime(1,ntax,0);
    //   for (int i=0; i<ntax; i++) {
    //       tipsasmatrix[i][0]=tips[i];
    //       tipsasmatrixprime[0][i]=tips[i];
    //   }
    //   gsl_matrix step4(ntax,1,0);
    //   step4=VCVinverse*tipsasmatrix;
    //   gsl_matrix step5(1,1,0);
    //   for (int i=0; i<ntax; i++) {
    //       step5[0][0]+=step4[i][0];
    //   }
    //   gsl_matrix finalvalue(1,1,0);
    //   finalvalue=step3*step5;
    //   ancestralstate=finalvalue[0][0];

	gsl_matrix_free(VCVinverse);
	gsl_matrix_free(VCVinversestart);
	gsl_matrix_free(tipsasmatrix);
	gsl_matrix_free(stepA);
	gsl_permutation_free(p);
    return ancestralstate;	
}


gsl_vector * BROWNIE::GetTipResiduals(gsl_vector * tips, double ancestralstate)
{
    int ntax=tips->size;
    gsl_vector *tipresiduals;
    tipresiduals=gsl_vector_calloc(ntax);
    //gsl_vector tipresiduals(ntax,0);
    //cout<<"GetTipResiduals ntax="<<ntax<<endl;
    //cout<<"Tips / tip residuals"<<endl;
    for (int i=0; i<ntax; i++) {
        // cout<<tips[i];
        //cout<<" / ";
        gsl_vector_set(tipresiduals,i,((gsl_vector_get(tips,i))-ancestralstate));
        //cout<<tipresiduals[i]<<endl;
    }
    return tipresiduals;
}

//Estimate rate
//Make sure to use tip residuals (tips minus estimated ancestral state)
double BROWNIE::EstimateRate(gsl_matrix * VCV, gsl_vector * tipresiduals)
{
    //rate=tipresiduals'*(inv(currenttreematrix))*tipresiduals/ntax
    double rateparameter;
    int ntax=VCV->size1;
    //matrixsingular=TestSingularity(VCV);
    //if (matrixsingular) {
    //    errormsg="Singular matrix (input VCV) during estimate rate";
    //   throw XNexus(errormsg);
    //}
    if (debugmode) {
        message="\nNow calculating inverse VCV for estimating rate.\n";
        PrintMessage();
    }
gsl_permutation * p = gsl_permutation_alloc (ntax);
int signum;
gsl_matrix * VCVinverse=gsl_matrix_calloc(ntax,ntax);
gsl_matrix * VCVinversestart=gsl_matrix_calloc(ntax,ntax);
gsl_matrix_memcpy (VCVinversestart, VCV);
if(debugmode) {
    message="debugging\n";
    for (int currentrow=0;currentrow<VCV->size1;currentrow++) {
        for (int currentcol=0;currentcol<VCV->size2;currentcol++) {
            message+=gsl_matrix_get(VCV,currentrow,currentcol);
            message+="\t";
        }
        message+="\n";
    }
    message+="\n\nduplicatedmatrix\n";
    for (int currentrow=0;currentrow<VCVinversestart->size1;currentrow++) {
        for (int currentcol=0;currentcol<VCVinversestart->size2;currentcol++) {
            message+=gsl_matrix_get(VCVinversestart,currentrow,currentcol);
            message+="\t";
        }
        message+="\n";
    }

    PrintMessage();
}
gsl_linalg_LU_decomp (VCVinversestart,p, &signum);
if(debugmode) {
    message="\n\nLUmatrix\n";
    for (int currentrow=0;currentrow<VCVinversestart->size1;currentrow++) {
        for (int currentcol=0;currentcol<VCVinversestart->size2;currentcol++) {
            message+=gsl_matrix_get(VCVinversestart,currentrow,currentcol);
            message+="\t";
        }
        message+="\n";
    }

    PrintMessage();
}

gsl_linalg_LU_invert (VCVinversestart,p, VCVinverse);

if(debugmode) {
    message="\n\nInverseMatrix\n";
    for (int currentrow=0;currentrow<VCVinverse->size1;currentrow++) {
        for (int currentcol=0;currentcol<VCVinverse->size2;currentcol++) {
            message+=gsl_matrix_get(VCVinverse,currentrow,currentcol);
            message+="\t";
        }
        message+="\n";
    }

    PrintMessage();
}


//VCVinverse=Inverse(VCV);
if (debugmode) {
    message="\nFinished calculating inverse VCV for estimating rate.\n";
    PrintMessage();
}
gsl_vector *step1vect;
step1vect=gsl_vector_calloc(ntax);
//gsl_vector step1vect(ntax,0);
gsl_blas_dgemv (CblasNoTrans,1, VCVinverse, tipresiduals,0, step1vect); ///TEST THIS
                                                                        //step1vect=MatrixTimesVector(VCVinverse,tipresiduals);
gsl_vector *step2vect;
step2vect=gsl_vector_calloc(1);
//gsl_vector step2vect(1,0);
gsl_matrix *tipsasmatrixprime=gsl_matrix_calloc(1,ntax);
for (int i=0; i<ntax; i++) {
    gsl_matrix_set(tipsasmatrixprime,0,i,(gsl_vector_get(tipresiduals,i)));
}
gsl_blas_dgemv (CblasNoTrans,1, tipsasmatrixprime, step1vect,0, step2vect); ///TEST THIS
                                                                            //step2vect=MatrixTimesVector(tipsasmatrixprime,step1vect);
rateparameter=gsl_vector_get(step2vect,0)/ntax;
if (debugmode) {
    message="rate estimate is ";
    message+=rateparameter;
    PrintMessage();
}
gsl_matrix_free(VCVinverse);
gsl_matrix_free(VCVinversestart);
gsl_vector_free(step1vect);
gsl_vector_free(step2vect);
gsl_matrix_free(tipsasmatrixprime);
gsl_permutation_free (p);
return rateparameter;
}

double BROWNIE::GetDiscreteCharLnL_gsl( const gsl_vector * variables, void *obj) 
{
	double temp;
	temp= ((BROWNIE*)obj)->GetDiscreteCharLnL(variables);
	if((gsl_isinf (temp))<0) {
		temp=BROWNIE_MAXLIKELIHOOD;
	}
	return temp;
}

double BROWNIE::GetDiscreteCharLnL(const gsl_vector * variables)
{
	if (debugmode) {
		cout<<endl<<endl<<"------   Using GetDiscreteCharLnL -------"<<endl<<endl;
	}
	gsl_vector *localvariables=gsl_vector_calloc(variables->size);
	gsl_vector_memcpy(localvariables,variables);
	if(nonnegvariables) { //N-M can get negative values for parameters. This is fine usually, but not with rates and frequencies, which must be nonnegative. Solution? NM variable x=log(true variable); true variable Y=exp(NM variable)
		for (int i=0;i<variables->size;i++) {
			gsl_vector_set(localvariables,i,exp(gsl_vector_get(localvariables,i)));
		}
	}
	//double likelihood=GSL_POSINF;
	double likelihood=BROWNIE_MAXLIKELIHOOD; //rather than an infinte value, use the maximum possible value, so numerical optimization doesn't fail
	negbounceparam=-1;
	if (numberoffreeparameters>0) { //if the input vector has useful variables; this number is only zero in the case of some user models
		if (gsl_vector_min(localvariables)<0) { //means we have a negative rate or state frequency if <0, so leave the likelihood set at a really bad number
			negbounceparam=gsl_vector_min_index(localvariables);
			if(detailedoutput) {
				cout<<"Had negative input, variables vector is ( ";
				for (int i=0;i<numberoffreeparameters;i++) {
					cout<<gsl_vector_get(localvariables,i)<<" ";
				}
				cout<<")"<<endl;
			}
			if (debugmode) {
				cout<<"GetDiscreteCharLnL output (first return) is "<<likelihood<<endl;
			}
			return likelihood;
		}
	}
	//else {
		//int numbercharstates=(discretecharacters->GetObsNumStates(discretechosenchar));
		int numberofrates=(localnumbercharstates*localnumbercharstates)-localnumbercharstates;
		int ntax=taxa->GetNumTaxonLabels();
		
		gsl_matrix *RateMatrix=gsl_matrix_calloc(localnumbercharstates,localnumbercharstates);
		gsl_vector *ancestralstatevector=gsl_vector_calloc(localnumbercharstates);
		int vectorposition=0;
		int position=-1; //used in user-set matrix only
		for (int i=0; i<localnumbercharstates;i++) {
			for (int j=0; j<localnumbercharstates;j++) {
				if (i!=j) {
					if (discretechosenmodel==1) { //one rate
						gsl_matrix_set(RateMatrix,i,j,gsl_vector_get(localvariables,0));
						vectorposition=1;
					}
					if (discretechosenmodel==2) { //rev
						if (i<j) {
							gsl_matrix_set(RateMatrix,i,j,gsl_vector_get(localvariables,vectorposition));
							vectorposition++;
							gsl_matrix_set(RateMatrix,j,i,gsl_matrix_get(RateMatrix,i,j)); //since symmetric
						}
					}
					if (discretechosenmodel==3) { //nonrev
						gsl_matrix_set(RateMatrix,i,j,gsl_vector_get(localvariables,vectorposition)); //all off-diagonal elements get their own rates
						vectorposition++;
					}
					if (discretechosenmodel==4) { //user
						if (debugmode) {
							cout<<"vectorposition = "<<vectorposition<<"ratematassignvector["<<vectorposition<<"]="<<ratematassignvector[vectorposition]<<endl;
						}
						if (ratematassignvector[vectorposition]>=0) { //means there's an assigned rate
							gsl_matrix_set(RateMatrix,i,j,ratematfixedvector[(ratematassignvector[vectorposition])]);
						}
						else {
							position=-1*(1+ratematassignvector[vectorposition]);
							if (debugmode) {
								cout<<"calling position "<<position<<" for variables vector of size "<<localvariables->size<<endl;
							}
							gsl_matrix_set(RateMatrix,i,j,gsl_vector_get(localvariables,position)); //means it's a variable rate
						}
						vectorposition++;
					}
				}
			}
		}
		if (discretechosenmodel==4) {
			vectorposition=position+1;
		}
		for (int i=0; i<localnumbercharstates; i++) { //fill in diagonal entries
			double ratesum=0;
			for (int j=0; j<localnumbercharstates; j++) {
				ratesum+=gsl_matrix_get(RateMatrix,i,j);
			}
			gsl_matrix_set(RateMatrix,i,i,-1.0*ratesum);
		}
		if (discretechosenstatefreqmodel==3) { //Calculate Equilibrium state freqs by just getting a Pmatrix (P=exp(QT)) for a really big time
			if (debugmode) {
				cout<<"now computing equilibrium state frequencies"<<endl;
			}
			gsl_matrix* Pmatrix=gsl_matrix_calloc(localnumbercharstates,localnumbercharstates);
			gsl_matrix* StartFreqs=gsl_matrix_calloc(1,localnumbercharstates);
			gsl_matrix_set_all (StartFreqs, 1.0/localnumbercharstates); //start with equal freqs
			gsl_matrix* EndFreqs=gsl_matrix_calloc(1,localnumbercharstates);
			gsl_matrix* EndFreqs2=gsl_matrix_calloc(1,localnumbercharstates);
			gsl_matrix_set_all(EndFreqs2,1.0/localnumbercharstates);
			double maxdiff=1.0;
			double tolerlimit=0.0001/localnumbercharstates; //More char states, want more precision
			double simulatedtime=1000;
			while (maxdiff>tolerlimit) {
				Pmatrix=ComputeTransitionProb(RateMatrix,simulatedtime); //A really long time
				simulatedtime*=100; //increase it in case we need to do this over a longer time interval
				gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, StartFreqs, Pmatrix,0.0, EndFreqs);
				gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, EndFreqs, Pmatrix,0.0, EndFreqs2);
				gsl_matrix_sub(EndFreqs,EndFreqs2); //leaves EndFreqs2 unchanged, stores diff in EndFreqs
				maxdiff=0.0;
				for (int i=0;i<localnumbercharstates;i++) {
					maxdiff=GSL_MAX(fabs(gsl_matrix_get(EndFreqs,0,i)),maxdiff);
				}
			}
			for (int i=0;i<localnumbercharstates;i++) {
				gsl_vector_set(ancestralstatevector,i,gsl_matrix_get(EndFreqs2,0,i));
			}
			gsl_matrix_free(Pmatrix);
			gsl_matrix_free(StartFreqs);
			gsl_matrix_free(EndFreqs);
			gsl_matrix_free(EndFreqs2);
		}
		else {	
			for (int i=0; i<localnumbercharstates; i++) { //do ancestralstatevector for freqs
				if (discretechosenstatefreqmodel==1) {
					if (debugmode) {
						cout<<"now setting uniform state frequencies"<<endl;
					}	
					//Uniform
					//gsl_vector_set(ancestralstatevector,i,1.0/localnumbercharstates);
					if (i<(localnumbercharstates-1)) {
						gsl_vector_set(ancestralstatevector,i,1.0/localnumbercharstates);
					}
					else { //last number must be 1-sum(other states)
						double frequencysum=0.0;
						for (int j=0; j<i;j++) {
							frequencysum+=gsl_vector_get(ancestralstatevector,j);
						}
						gsl_vector_set(ancestralstatevector,i,1.0-frequencysum);
					}
				}
				else if (discretechosenstatefreqmodel==2) {
					if (debugmode) {
						cout<<"now setting state frequencies based on empirical frequencies for this character"<<endl;
					}	
					double frequency=0.0;
					for (int j=0;j<ntax;j++) {
						if (discretecharacters->GetInternalRepresentation(j,discretechosenchar)==i) {
							frequency+=1.0/ntax;
						}
					}
					gsl_vector_set(ancestralstatevector,i,frequency);
				}
				else if (discretechosenstatefreqmodel==4) {
					if (debugmode) {
						cout<<"now setting state frequencies based on optimal frequencies for this character"<<endl;
					}	

					if (debugmode) {
						cout<<"i="<<i<<" variables->size="<<localvariables->size<<" vectorposition="<<vectorposition<<endl;
					}
					if (i<(localnumbercharstates-1)) {
						if (debugmode) {
							cout<<"vectorposition = "<<vectorposition<<endl;
						}
						gsl_vector_set(ancestralstatevector,i,gsl_vector_get(localvariables,vectorposition));
						vectorposition++;
					}
					else { //last number must be 1-sum(other states)
						double frequencysum=0.0;
						for (int j=0; j<i;j++) {
							frequencysum+=gsl_vector_get(ancestralstatevector,j);
						}
						gsl_vector_set(ancestralstatevector,i,1.0-frequencysum);
					}
				}
				else if (discretechosenstatefreqmodel==5) {
					//user
					if (debugmode) {
						cout<<"now setting state frequencies based on user-specified frequencies for this character"<<endl;
					}	

					if (userstatefreqvector.size()!=localnumbercharstates) {
						errormsg="The current (possibly default) vector of user-specified state frequencies ";
						errormsg+="\nwith size ";
						int freqvectorsize=userstatefreqvector.size();
						errormsg+=freqvectorsize;
						errormsg+=" should have size ";
						errormsg+=localnumbercharstates;
						errormsg+=" for the selected character";
						throw XNexus( errormsg);
					}
					gsl_vector_set(ancestralstatevector,i,userstatefreqvector[i]);
				}
			}
		}
		/*cout<<"\n\nbbbbbbbbbbbbbbbbbbbbbbb"<<endl;
		*/
		if (debugmode) {
			PrintMatrix(RateMatrix);
			cout<<"\nancestral freq: ";
			for (int i=0; i<ancestralstatevector->size;i++) {
				cout<<gsl_vector_get(ancestralstatevector,i)<<"\t";
			}
			cout<<endl;
		}
		double frequencysum=0.0;
		for (int i=0; i<ancestralstatevector->size;i++) {
			frequencysum+=gsl_vector_get(ancestralstatevector,i);
		}
//		if ((fabs(frequencysum-1.0)>BROWNIE_EPSILON) || gsl_vector_min(ancestralstatevector)<0 || gsl_vector_max(ancestralstatevector)>1) { //a way of constraining the search
		if ((gsl_fcmp(frequencysum,1.0,BROWNIE_EPSILON)!=0) || gsl_vector_min(ancestralstatevector)<0 || gsl_vector_max(ancestralstatevector)>1) { //a way of constraining the search
			if(detailedoutput) {
				cout<<"\n\tHad wrong state frequency input, ancestralstatevector vector is ( ";
				for (int i=0;i<ancestralstatevector->size;i++) {
					cout<<gsl_vector_get(ancestralstatevector,i)<<" ";
				}
				cout<<") ";
				if (gsl_fcmp(frequencysum,1.0,BROWNIE_EPSILON)!=0) {
					cout<<"[SUM ("<<frequencysum<<") NOT 1: Sum-1 = "<<fabs(frequencysum-1.0)<<", gsl_fcmp(frequencysum,1.0,BROWNIE_EPSILON) = "<<gsl_fcmp(frequencysum,1.0,BROWNIE_EPSILON)<<"] ";
				}
				if (gsl_vector_min(ancestralstatevector)<0) {
					cout<<"[MIN ("<<gsl_vector_min(ancestralstatevector)<<") < 0] ";
				}
				if (gsl_vector_max(ancestralstatevector)>1 ) {
					cout<<"[MAX ("<<gsl_vector_max(ancestralstatevector)<<") > 1] ";
				}
				cout<<endl;
			}
			//likelihood=GSL_POSINF;
			likelihood=BROWNIE_MAXLIKELIHOOD; //rather than an infinte value, use the maximum possible value, so numerical optimization doesn't fail

			if (debugmode) {
				cout<<"GetDiscreteCharLnL output (second return) is "<<likelihood<<endl;
			}
			return likelihood;	
		}
		
		likelihood=(CalculateDiscreteCharLnL(RateMatrix,ancestralstatevector));
		//cout<<"likelihood is "<<likelihood<<endl;
		gsl_matrix_swap(currentdiscretecharQmatrix,RateMatrix);
		//cout<<"currentdiscretecharQmatrix\n"; 
		//PrintMatrix(currentdiscretecharQmatrix);
		gsl_vector_swap(currentdiscretecharstatefreq,ancestralstatevector);
		gsl_matrix_free(RateMatrix);
		gsl_vector_free(ancestralstatevector);
		gsl_vector_free(localvariables);
		//cout<<"\neeeeeeeeeeeeeeeeeeeeeeeeeee"<<endl;
		if (debugmode) {
			cout<<"GetDiscreteCharLnL output (third return) is "<<likelihood<<endl;
		}
		return likelihood;
	//}
}



////////////////////////////////////////////////
double BROWNIE::GetLikelihoodUnderLindy2_gsl( const gsl_vector * variables, void *obj) 
{
	double temp;
	temp= ((BROWNIE*)obj)->GetLikelihoodUnderLindy2(variables);
	return temp;
}

double BROWNIE::GetLikelihoodUnderLindy2(const gsl_vector * variables)
{
	double rateA=gsl_vector_get(variables,0);
	double rateB=gsl_vector_get(variables,1);
	double likelihood=(CalculateDiscreteLindy2(rateA,rateB));
	return likelihood;
}

double BROWNIE::GetLikelihoodUnderLindy1_gsl( const gsl_vector * variables, void *obj) 
{
	double temp;
	temp= ((BROWNIE*)obj)->GetLikelihoodUnderLindy1(variables);
	return temp;
}

double BROWNIE::GetLikelihoodUnderLindy1(const gsl_vector * variables)
{
	double rateA=gsl_vector_get(variables,0);
	double likelihood=(CalculateDiscreteLindy1(rateA));
	return likelihood;
}

gsl_vector * BROWNIE::LindyGeneralOptimization(int ChosenModel)
{
	//Model 1=One param
	//Model 2=Two param
	double bestlikelihood=GSL_POSINF;
	int hitlimitscount=0;
	size_t np;
	if (ChosenModel==1) {
		np = 1;
	}
	else if (ChosenModel==2) {
		np = 2;
	}
	gsl_vector * results=gsl_vector_calloc(np);
	double estimates[randomstarts][np];
	double startingvalues[randomstarts][np];
	double likelihoods[randomstarts][1];
	for (int startnum=0;startnum<randomstarts;startnum++) {
		const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
		gsl_multimin_fminimizer *s = NULL;
		gsl_vector *ss, *x;
		size_t iter = 0, i;
		int status;
		double size;
		bool hitlimits=false;
		/* Initial vertex size vector */
		ss = gsl_vector_alloc (np);
		gsl_vector_set_all (ss, stepsize);
		
		/* Starting point */
		x = gsl_vector_calloc (np);
		if (ChosenModel==1) {
			gsl_vector_set (x,0,gsl_ran_exponential (r,1.0)); //starting rate
			startingvalues[startnum][0]=gsl_vector_get(x,0);
		}
		else if (ChosenModel==2) {
			gsl_vector_set(x,0,gsl_ran_exponential (r,1.0)); //starting rate
			startingvalues[startnum][0]=gsl_vector_get(x,0);
			gsl_vector_set(x,1,gsl_ran_exponential (r,1.0)); 
			startingvalues[startnum][1]=gsl_vector_get(x,1);
		}
		BROWNIE *pt;
		pt=(this);
		double (*F)(const gsl_vector *, void *);
		if (ChosenModel==1) {
			F = &BROWNIE::GetLikelihoodUnderLindy1_gsl;
		}
		else if (ChosenModel==2) {
			F = &BROWNIE::GetLikelihoodUnderLindy2_gsl;
		}
		gsl_multimin_function minex_func;
		minex_func.f=*F;
		minex_func.params=pt;
		minex_func.n = np;
		s = gsl_multimin_fminimizer_alloc (T, np);
		gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
		do
		{
			iter++;
			status = gsl_multimin_fminimizer_iterate(s);
			if (status!=0) { //0 Means it's a success in c++, but not in C
				printf ("error: %s\n", gsl_strerror (status));
				break;
			}
			size = gsl_multimin_fminimizer_size (s);
			status = gsl_multimin_test_size (size, stoppingprecision); //since we want more precision
			if (status == GSL_SUCCESS)
			{
				//printf ("converged to minimum at\n");
			}
		}
		while (status == GSL_CONTINUE && iter < maxiterations);
		if (s->fval<bestlikelihood) {
			gsl_vector_memcpy(results,s->x);
			bestlikelihood=s->fval;
		}
		for (int parameternumber=0; parameternumber<np; parameternumber++) {
			estimates[startnum][parameternumber]=gsl_vector_get(s->x,parameternumber);
		}
		if (iter==maxiterations) {
			hitlimits=true;
			hitlimitscount++;
		}
		message="Replicate ";
		message+=startnum+1;
		if (hitlimits) {
			message+=" **WARNING**";
		}
		message+="\n   NM iterations needed = ";
		int iterationsrequired=iter;
		message+=iterationsrequired;
		if (hitlimits) {
			message+=" **Max iterations hit; see WARNING below**";
		}
		message+="\n   -LnL = ";
		char outputstring[60];
		sprintf(outputstring,"%60.45f",1.0*(s->fval));
		message+=outputstring;
		// message+="\n   Rate = ";
		// message+=gsl_vector_get(s->x,0);
		if (detailedoutput) {
			PrintMessage();
		}
		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free (s);
	//	if(detailedoutput==false) {
	//		ProgressBar(0);
	//	}
	}
	if (hitlimitscount>0) {
		message="\n----------------------------------------------------------------------------\n WARNING: Out of ";
		message+=randomstarts;
		message+=" optimization starts, ";
		message+=hitlimitscount;
		if (hitlimitscount==1) {
			message+=" was ";
		}
		else {
			message+=" were ";
		}
		message+="stopped by hitting\n  the maximum # of iterations. This means that those replicates\n  may not even have hit the local maximum.\n\n  You can increase the maximum number of iterations or decrease the\n  precision with the NumOpt command. You could also consider\n  increasing the number of random starts using that same command.\n\n  If this happened on a small proportion of replicates, though,\n  or if the precision (below) is good enough, don't worry about it.\n----------------------------------------------------------------------------";
		if (detailedoutput) {
			PrintMessage();
		}
		else if ((randomstarts-hitlimitscount)<10 && (hitlimitscount/randomstarts)>.1) {
			PrintMessage();
		}
	}
	gsl_vector * finalvector=gsl_vector_calloc((2*np)+1);
	for (int position=0; position<np; position++) {
		gsl_vector_set(finalvector,position,gsl_vector_get(results,position));
		double paramestimate[randomstarts];
		for (int startnumber=0;startnumber<randomstarts;startnumber++) {
			paramestimate[startnumber]=estimates[startnumber][position];
		}
		gsl_vector_set(finalvector,position+np,gsl_stats_sd(paramestimate,1,randomstarts));
	}
	gsl_vector_set(finalvector,(2*np),bestlikelihood);
	return finalvector;
}

gsl_vector * BROWNIE::DiscreteGeneralOptimization()
{
	bestdiscretelikelihood=BROWNIE_MAXLIKELIHOOD;
	bool globalbesthadfixedzerosorones=false;
	int hitlimitscount=0;
	//int numbercharstates=(discretecharacters->GetObsNumStates(discretechosenchar));
	int numberofrates=(localnumbercharstates*localnumbercharstates)-localnumbercharstates;
	int ntax=taxa->GetNumTaxonLabels();
	numberoffreerates=0;
	numberoffreefreqs=0;
	//gsl_matrix_free(optimaldiscretecharQmatrix);
	optimaldiscretecharQmatrix=gsl_matrix_calloc(localnumbercharstates,localnumbercharstates);
	currentdiscretecharQmatrix=gsl_matrix_calloc(localnumbercharstates,localnumbercharstates);
	//gsl_vector_free(optimaldiscretecharstatefreq);
	optimaldiscretecharstatefreq=gsl_vector_calloc(localnumbercharstates);
	currentdiscretecharstatefreq=gsl_vector_calloc(localnumbercharstates);
	double currentstepsize=stepsize;

	if (discretechosenmodel==1) {
		numberoffreerates=1;
	}
	else if (discretechosenmodel==2) {
		numberoffreerates=numberofrates/2;
	}
	else if (discretechosenmodel==3) {
		numberoffreerates=numberofrates;
	}
	else if (discretechosenmodel==4) {
		numberoffreerates=freerateletterstring.size();
	}
	if (discretechosenstatefreqmodel==4) {
		numberoffreefreqs=localnumbercharstates-1;
	}
	numberoffreeparameters=numberoffreerates+numberoffreefreqs; //stored globally for later calculation of AIC/AICc
	if (((1.0*ntax)/(1.0*numberoffreeparameters))<10 && !allchar) {
		message="\n-------------------------------------------------------------------------------\n WARNING: You are trying to estimate ";
		message+=numberoffreeparameters;
		if (numberoffreeparameters==1) {
			message+=" parameter, but only have ";
		}
		else {
			message+=" parameters, but only have ";
		}
		message+=ntax; 
		message+=" taxa.\n Make sure to try some simpler models, and expect quite imprecise estimates.\n-------------------------------------------------------------------------------";
		PrintMessage();
	}
	else if (((discretecharacters->GetNChar())*((1.0*ntax)/(1.0*numberoffreeparameters)))<10 && allchar) {
		message="\n-------------------------------------------------------------------------------\n WARNING: You are trying to estimate ";
		message+=numberoffreeparameters;
		if (numberoffreeparameters==1) {
			message+=" parameter, but only have ";
		}
		else {
			message+=" parameters, but only have ";
		}
		message+=ntax; 
		message+=" taxa and ";
		message+=discretecharacters->GetNChar();
		message+=" characters.\n Make sure to try some simpler models, and expect quite imprecise estimates.\n-------------------------------------------------------------------------------";
		PrintMessage();
	}
	
	size_t np=numberoffreeparameters;
	gsl_vector * finalvector=gsl_vector_calloc((2*np)+1);
	if (np>0) {
		gsl_vector * results=gsl_vector_calloc(np);
		double estimates[randomstarts][np];
		double startingvalues[randomstarts][np];
		double likelihoods[randomstarts][1];
		if(detailedoutput==false) {
			ProgressBar(randomstarts);
		}
		for (int startnum=0;startnum<randomstarts;startnum++) {
			if (optimizationalgorithm==1) { //do nelder-mead simplex
				const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
				gsl_multimin_fminimizer *s = NULL;
				gsl_vector *ss, *x;
				size_t iter = 0, i;
				int status;
				double size;
				bool hitlimits=false;
				/* Initial vertex size vector */
				ss = gsl_vector_alloc (np);
				gsl_vector_set_all (ss, currentstepsize);
				
				/* Starting points */
				x = gsl_vector_calloc (np);
				if (startnum==0 || (gsl_ran_flat(r,0,1))>0.5 ) { //about 50% of the time, start from these values
					for (int i=0; i<numberoffreerates; i++) {
						gsl_vector_set (x,i,GSL_MIN(gsl_ran_exponential (r,0.5),gsl_ran_flat(r,0,1) )); //starting rate
						startingvalues[startnum][i]=gsl_vector_get(x,i);
					}
					for (int i=numberoffreerates; i<numberoffreerates+numberoffreefreqs; i++) {  
						gsl_vector_set (x,i,(1.0/localnumbercharstates)); //starting freqs are equal
						startingvalues[startnum][i]=gsl_vector_get(x,i);			
					}
				}
				else { //start again from near current point
					for (int i=0; i<numberoffreerates; i++) { 
						gsl_vector_set (x,i,gsl_ran_exponential(r,(estimates[startnum-1][i]))); //use a modified optimal value from the last run
						startingvalues[startnum][i]=gsl_vector_get(x,i);
					}
					for (int i=numberoffreerates; i<numberoffreerates+numberoffreefreqs; i++) {
						gsl_vector_set (x,i,(estimates[startnum-1][i])); //use optimal value from the last run
						startingvalues[startnum][i]=gsl_vector_get(x,i);			
					}			
				}
				if(nonnegvariables) { //N-M can get negative values for parameters. This is fine usually, but not with rates and frequencies, which must be nonnegative. Solution? NM variable x=log(true variable); true variable Y=exp(NM variable)
					for (int i=0; i<x->size; i++) {
						gsl_vector_set(x,i,log(gsl_vector_get(x,i)));
					}
				}
				BROWNIE *pt;
				pt=(this);
				double (*F)(const gsl_vector *, void *);
				F = &BROWNIE::GetDiscreteCharLnL_gsl;
				gsl_multimin_function minex_func;
				minex_func.f=*F;
				minex_func.params=pt;
				minex_func.n = np;
				s = gsl_multimin_fminimizer_alloc (T, np);
				gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
				do
				{
					iter++;
				/*	if (negbounceparam>-1) { //Means the previous iteration went too low at that position; we'll set the value to zero before iterating again
						if (detailedoutput) {
							cout<<"Had a negative value; now changing s->x from"<<endl;
							PrintVector(s->x);
							gsl_vector_set(s->x,negbounceparam,0.0);
							cout<<endl<<"to"<<endl;
							PrintVector(s->x);
						}
					} */
					status = gsl_multimin_fminimizer_iterate(s);
					if (status!=0) { //0 Means it's a success in c++, but not in C
						printf ("error: %s\n", gsl_strerror (status));
						break;
					}
					size = gsl_multimin_fminimizer_size (s);
					status = gsl_multimin_test_size (size, stoppingprecision); //since we want more precision
					if (status == GSL_SUCCESS)
					{
					//printf ("converged to minimum at\n");
					}
					
					if (detailedoutput) {
						if (iter<100 || (iter%25 ==0)) {
							printf ("%5d ", iter);
							for (i = 0; i < np; i++)
							{
								if (nonnegvariables) {
									printf ("%10.9e ", exp(gsl_vector_get (s->x, i)));
								}
								
								else {
									printf ("%10.9e ", gsl_vector_get (s->x, i));
								}
							}
							printf ("f() = %7.9f size = %.9f\n", s->fval, size);	
						}
					}
				}
				while (status == GSL_CONTINUE && iter < maxiterations);
				
			//cout<<"Current BEST matrix\n\n";
			//PrintMatrix(optimaldiscretecharQmatrix);
			//cout<<"\n\nLast matrix\n";
			//cout<<"Old best likelihood = ";
			//cout<<bestlikelihood;
			//cout<<"\n    Last likelihood = ";
			//cout<<s->fval;
			//cout<<endl<<endl;
			//PrintMatrix(currentdiscretecharQmatrix);
				if (s->fval<bestdiscretelikelihood) {
					globalbesthadfixedzerosorones=false;
					gsl_vector_memcpy(results,s->x);
					if(nonnegvariables) { //N-M can get negative values for parameters. This is fine usually, but not with rates and frequencies, which must be nonnegative. Solution? NM variable x=log(true variable); true variable Y=exp(NM variable)
						for (int i=0; i<x->size; i++) {
							gsl_vector_set(results,i,exp(gsl_vector_get(results,i)));
						}
					}					
					bestdiscretelikelihood=s->fval;
				//if (debugmode) {
				//	cout<<"likelihood score improved"<<endl;
				//}
					gsl_matrix_swap(currentdiscretecharQmatrix,optimaldiscretecharQmatrix);
					gsl_vector_swap(currentdiscretecharstatefreq,optimaldiscretecharstatefreq);
				}
				//The following section tries to round the parameter values: it's possible that 0 is a better value than some very small double
				gsl_vector *roundx=gsl_vector_calloc((s->x)->size);
				gsl_vector_memcpy(roundx,s->x);
				for (int i=0;i<roundx->size;i++) {
					if(nonnegvariables) {
						gsl_vector_set(roundx,i,exp(gsl_vector_get(roundx,i)));
					}
					if (gsl_vector_get(roundx,i)<8.0*BROWNIE_EPSILON) {
						gsl_vector_set(roundx,i,0.0);
					}
					else if (fabs(1.0-gsl_vector_get(roundx,i))<8.0*BROWNIE_EPSILON) {
						gsl_vector_set(roundx,i,1.0);
					}
				}
				bool roundedwasbetter=false;
				bool orignonnegvariables=nonnegvariables;
				nonnegvariables=false; //since roundx is untransformed
				double roundlikelihood=GetDiscreteCharLnL(roundx);
				nonnegvariables=orignonnegvariables;
				if (roundlikelihood<bestdiscretelikelihood) {
					roundedwasbetter=true;
					globalbesthadfixedzerosorones=true;
					bestdiscretelikelihood=roundlikelihood;
					gsl_matrix_swap(currentdiscretecharQmatrix,optimaldiscretecharQmatrix);
					gsl_vector_swap(currentdiscretecharstatefreq,optimaldiscretecharstatefreq);
					gsl_vector_memcpy(results,roundx);
				}
				if (detailedoutput) {
						printf ("fixed zeros %5d ", iter);
						for (i = 0; i < np; i++)
						{
				
								printf ("%10.9e ", gsl_vector_get (roundx, i));

						}
						printf ("f() = %7.9f \n", roundlikelihood);	
				}
				
				
				if (iter==maxiterations) {
					hitlimits=true;
					hitlimitscount++;
				}
				message="Replicate ";
				if (detailedoutput) {
					message+=startnum+1;
					if (hitlimits) {
						message+=" **WARNING**";
					}
					message+="\n   NM iterations needed = ";
					int iterationsrequired=iter;
					message+=iterationsrequired;
					if (hitlimits) {
						message+=" **Max iterations hit; see WARNING below**";
					}
					message+="\n   -LnL = ";
					char outputstring[60];
					sprintf(outputstring,"%60.45f",1.0*(s->fval));
					if (roundedwasbetter) {
						sprintf(outputstring,"%60.45f",1.0*roundlikelihood);					
					}
					message+=outputstring;
					message+="\n   Starts:    ";
					for (int parameternumber=0; parameternumber<np; parameternumber++) {
						message+=startingvalues[startnum][parameternumber];
						message+=" ";
					}				
					message+="\n   Estimates: ";
				}
				for (int parameternumber=0; parameternumber<np; parameternumber++) {
					if(roundedwasbetter) {
						estimates[startnum][parameternumber]=gsl_vector_get(roundx,parameternumber);
					}
					else {
						if(nonnegvariables) { //N-M can get negative values for parameters. This is fine usually, but not with rates and frequencies, which must be nonnegative. Solution? NM variable x=log(true variable); true variable Y=exp(NM variable)
							estimates[startnum][parameternumber]=exp(gsl_vector_get(s->x,parameternumber));
						}	
						else {
							estimates[startnum][parameternumber]=gsl_vector_get(s->x,parameternumber);
						}
					}
					if (detailedoutput) {
						message+=estimates[startnum][parameternumber];
						message+=" ";
					}
				}		
			// message+="\n   Rate = ";
			// message+=gsl_vector_get(s->x,0);
				if (detailedoutput) {
					PrintMessage();
				}
				gsl_vector_free(x);
				gsl_vector_free(ss);
				gsl_vector_free(roundx);
				gsl_multimin_fminimizer_free (s);
				if (hitlimits && redobad) {
					startnum--; //redo this rep (from a new starting point)
					if (hitlimitscount>giveupfactor*randomstarts) {
						message="\n----------------------------------------------------------------------------\n";
						message+= " ABORTING: You have chosen to keep restarting until you get ";
						message+=randomstarts;
						message+=" to\n";
						message+=" complete, but we've already tried ";
						message+=hitlimitscount;
						message+=" and only completed\n ";
						message+=startnum+1;
						message+=" starts. Maybe this is enough for you?\n You can change optimization settings with the NumOpt command.\n----------------------------------------------------------------------------";
						PrintMessage();
						startnum=randomstarts+1; // So we'll stop the run
					}
				}
				else {
					if (detailedoutput==false) {
						ProgressBar(0);
					}			
				}
			}
			else if (optimizationalgorithm==2) { //do simulated annealing
			//This uses the algorithm from the GSL siman.c, though completely rewritten
			//use info from example:
//			int ntries=200;             /* how many points do we try before stepping */
				int iters_fixed_t=200;      /* how many iterations for each T? */
				double K=1.0;                   /* Boltzmann constant */
				double mu_t=1.1;              /* damping factor for temperature */
				double t_min=stoppingprecision;
				gsl_vector *x = gsl_vector_calloc (np);
				double lastlikelihood=BROWNIE_MAXLIKELIHOOD;
				int numTEststarts=50;
				gsl_vector *initiallikelihoodscores=gsl_vector_calloc(numTEststarts);
				for (int Tstartrep=0; Tstartrep<numTEststarts; Tstartrep++) { //this gives us an estimate of what a good starting T would be
					lastlikelihood=BROWNIE_MAXLIKELIHOOD;
					while(isinf(lastlikelihood)) {
						for (int i=0; i<numberoffreerates; i++) {
							gsl_vector_set (x,i,GSL_MIN(gsl_ran_exponential (r,0.5),gsl_ran_flat(r,0,1) )); //starting rate
							startingvalues[startnum][i]=gsl_vector_get(x,i);
						}
						for (int i=numberoffreerates; i<numberoffreerates+numberoffreefreqs; i++) {  
							gsl_vector_set (x,i,(1.0/localnumbercharstates)); //starting freqs are equal
							startingvalues[startnum][i]=gsl_vector_get(x,i);			
						}
						if(nonnegvariables) { //N-M can get negative values for parameters. This is fine usually, but not with rates and frequencies, which must be nonnegative. Solution? NM variable x=log(true variable); true variable Y=exp(NM variable)
							for (int i=0; i<x->size; i++) {
								gsl_vector_set(x,i,log(gsl_vector_get(x,i)));
							}
						}
						
						lastlikelihood=GetDiscreteCharLnL(x);
					}
					gsl_vector_set(initiallikelihoodscores,Tstartrep,lastlikelihood);
				}
				double t_initial=2.0*(gsl_vector_max(initiallikelihoodscores)-gsl_vector_min(initiallikelihoodscores));       /* initial temperature */
				
				double T=t_initial;
				while (T>=t_min) {
					for (int j=0;j<iters_fixed_t;j++) {
						gsl_vector * newx=gsl_vector_calloc(numberoffreeparameters);
						gsl_vector_memcpy(newx,x);
						int i=int(gsl_ran_flat(r,0,numberoffreeparameters)); //changes only one parameter, we pick it randomly
						if (i<numberoffreerates) { //must be changing a rate; rates must be >=0 but no upper bound, so use an exponential distribution
							gsl_vector_set(newx,i,0.99*gsl_vector_get(newx,i)+0.01*gsl_ran_exponential (r,1.0/gsl_vector_get(newx,i)));
						}
						else {
							double p[numberoffreefreqs+1];
							unsigned int n[numberoffreefreqs+1];
							double totalfreq=0;
							for (int j=0;j<numberoffreefreqs;j++) {
								p[j]=gsl_vector_get(x,j);
								totalfreq+=p[j];
							}
							p[numberoffreefreqs]=1.0-totalfreq;
							int multinomialsamples=250;
							gsl_ran_multinomial (r, numberoffreefreqs+1, multinomialsamples, p, n); //modify the state freqs slightly
							for (int j=numberoffreerates;j<numberoffreeparameters;j++) {
								gsl_vector_set(newx,j,(1.0*n[j-numberoffreerates])/(1.0*multinomialsamples));
							}	
						}
						if (detailedoutput) {
							cout<<j<<endl<<"\t";
							for (int j=0;j<numberoffreeparameters;j++) {
								cout<<gsl_vector_get(x,j)<<"\t";
							}
							cout<<endl<<"\t";
							for (int j=0;j<numberoffreeparameters;j++) {
								cout<<gsl_vector_get(newx,j)<<"\t";
							}
							
						}
						double newlikelihood=GetDiscreteCharLnL(newx);
						if(isinf(T)) {
							T=0.5*newlikelihood;
						}
						if (newlikelihood<lastlikelihood || ( gsl_rng_uniform(r) > ((newlikelihood-lastlikelihood)/(K * T))) ) {
						//take the move
							gsl_vector_memcpy(x,newx);
							lastlikelihood=newlikelihood;
							if (newlikelihood<bestdiscretelikelihood) {
								gsl_vector_memcpy(results,newx);
								bestdiscretelikelihood=newlikelihood;
								gsl_matrix_swap(currentdiscretecharQmatrix,optimaldiscretecharQmatrix);
								gsl_vector_swap(currentdiscretecharstatefreq,optimaldiscretecharstatefreq);
							}
						}
						if(detailedoutput) {
							cout<<endl<<"\t"<<T<<" "<<newlikelihood<<" "<<lastlikelihood<<" ";
							for (int i=0;i<numberoffreeparameters;i++) {
								cout<<gsl_vector_get(newx,i)<<" ";
							}
							cout<<endl;
						}
					}
					T /= mu_t;
				}
				
				for (int parameternumber=0; parameternumber<numberoffreeparameters; parameternumber++) {
					if(nonnegvariables) { //N-M can get negative values for parameters. This is fine usually, but not with rates and frequencies, which must be nonnegative. Solution? NM variable x=log(true variable); true variable Y=exp(NM variable)
						estimates[startnum][parameternumber]=exp(gsl_vector_get(x,parameternumber));
					}	
					else {
						estimates[startnum][parameternumber]=gsl_vector_get(x,parameternumber);
					}
				}
				if (detailedoutput==false) {
					ProgressBar(0);
				}	
				gsl_vector_free(x);
			}
		}
		
		if (hitlimitscount>0) {
			if (redobad) {
				message="\n----------------------------------------------------------------------------\n WARNING: For ";
				message+=randomstarts;
				message+=" desired completed runs, I had to do ";
				message+=hitlimitscount;
				message+=" starts total.\n This could suggest that you should change some NumOpt options\n (\"numopt ?\" for help).\n----------------------------------------------------------------------------";
			}
			else {
				message="\n----------------------------------------------------------------------------\n WARNING: Out of ";
				message+=randomstarts;
				message+=" optimization starts, ";
				message+=hitlimitscount;
				if (hitlimitscount==1) {
					message+=" was ";
				}
				else {
					message+=" were ";
				}
				message+="stopped by hitting\n  the maximum # of iterations. This means that those replicates\n  may not even have hit the local maximum.\n\n  You can increase the maximum number of iterations or decrease the\n  precision with the NumOpt command. You could also consider\n  increasing the number of random starts using that same command.\n\n  If this happened on a small proportion of replicates, though,\n  or if the precision (below) is good enough, don't worry about it.\n----------------------------------------------------------------------------";
			}
			if (detailedoutput) {
				PrintMessage();
			}
			else if ((randomstarts-hitlimitscount)<10 && (hitlimitscount/randomstarts)>.1) {
				PrintMessage();
			}
		}
		for (int position=0; position<np; position++) {
			gsl_vector_set(finalvector,position,gsl_vector_get(results,position));
			double paramestimate[randomstarts];
			for (int startnumber=0;startnumber<randomstarts;startnumber++) {
				paramestimate[startnumber]=estimates[startnumber][position];
			}
			gsl_vector_set(finalvector,position+np,gsl_stats_sd(paramestimate,1,randomstarts));
		}
		gsl_vector_free(results);
	}
else {
		gsl_vector* emptyvector=gsl_vector_calloc(1);
		bestdiscretelikelihood=GetDiscreteCharLnL(emptyvector);	
		gsl_matrix_swap(currentdiscretecharQmatrix,optimaldiscretecharQmatrix);
		gsl_vector_swap(currentdiscretecharstatefreq,optimaldiscretecharstatefreq);
		gsl_vector_free(emptyvector);
}

/*if (globalbesthadfixedzerosorones) {
	message="\n--------------------------------------------------------------------------------------------\n";
	message+="Note that the best model estimates had some free parameter estimates fixed to zero or one;\n";
	message+="you may get a more precise estimate of other parameters by inputing such a model (but\n ";
	message+="use the current number of free parameters when calculating AIC/AICc).\n";
	message+="--------------------------------------------------------------------------------------------";
	PrintMessage();
}*/
	gsl_vector_set(finalvector,(2*np),bestdiscretelikelihood);
	//cout<<"Best rate matrix=\n";
	//PrintMatrix(optimaldiscretecharQmatrix);
	return finalvector;
}

/*double BROWNIE::E_discretegeneral(void *xp) {
	gsl_vector * inputvariables=gsl_vector_calloc(numberoffreeparameters);
	double *params = (double *) xp;
	for (int i=0;i++;i<numberoffreeparameters) {
		gsl_vector_set(inputvariables,i,params[i]);
	}
	double E=GetDiscreteCharLnL(inputvariables);
	gsl_vector_free(inputvariables);
	return E;
}

double BROWNIE::M_discretegeneral(void *xp, void *yp)
{
	double *params1 = (double *) xp, *params2 = (double *) yp;
	double distance = 0;
	int i;
	for (i = 0; i < numberoffreeparameters; ++i) {
		distance += GSL_MAX(params1[i],params2[i])-GSL_MAX(params1[i],params2[i]);
	}
	return distance;
}

void BROWNIE::S_discretegeneral(const gsl_rng * r, void *xp, double step_size)
{
	double old_x = *((double *) xp);
	double new_x=old_x;
	step_size=0;
	int i=int(gsl_ran_flat(r,0,numberoffreeparameters)); //changes only one parameter, we pick it randomly
	if (i<numberoffreerates) { //must be changing a rate; rates must be >=0 but no upper bound, so use an exponential distribution
		new_x[i]=0.9*old_x[i]+0.1*gsl_ran_exponential (r,1.0/old_x[i]);
	}
	else {
		double p[numberoffreefreqs+1];
		int n[numberoffreefreqs+1];
		double totalfreq=0;
		for (int j=0;j++;j<numberoffreefreqs) {
			p[j]=old_x[j];
			totalfreq+=old_x[j];
		}
		totalfreq[numberoffreefreqs]=1.0-totalfreq;
		int multinomialsamples=250;
		gsl_ran_multinomial (r, numberoffreefreqs+1, multinomialsamples, p, n); //modify the state freqs slightly
		for (int j=numberoffreerates;j++;j<numberoffreeparameters) {
			new_x[j]=(1.0*n[j-numberoffreerates])/(1.0*multinomialsamples);
		}	
	}
	
	memcpy(xp, &new_x, sizeof(new_x));
	
}
*/

//Calculates the joint ML mapping of ancestral states of a discrete on the tree, using the Pupko et al 2000 algorithm. Returns a pointer to the root of new tree with 
//ML  joint estimates of state orders and state times inferred; node labels give the reconstructed states. If breaksperbranch==0, each branch gets
//just one state assigned; if breaks per branch>0, extra nodes with degree 2 are introduced on the branch and states reconstructed at these nodes;
//this allows for a branch to be partly in one state and then change to another (if a branch starts in state 0, and ends in state 1, and the 0->1
//rate is lower than the 1->0 rate, for example, this will allow most of the branch to be reconstructed as having state 0).
NodePtr BROWNIE::EstimateMLDiscreteCharJointAncestralStates(gsl_matrix * RateMatrix, gsl_vector * ancestralstatevector, int breaksperbranch) {
	double lnL=0.0;
	double L=0;
	map<Node*, vector<vector<double> > > Lvector; //need vector of vectors to deal with breaks per branch
	map<Node*, vector<vector<int> > > Cvector; // gives C vector, as in Pupko et al algorithm
	//Tree T=intrees.GetIthTree(chosentree-1);
	Tree *Tptr=&(intrees.Trees[chosentree-1]);
	(*Tptr).Update();
	(*Tptr).GetNodeDepths();
	//(*Tptr).Draw(cout);
	//Tree OldFormat=(*Tptr);
	NodeIterator <Node> n ((*Tptr).GetRoot()); //Goes from tips down
	NodePtr currentnode = n.begin();
	while (currentnode)
	{
		double eachsegmentlength=(currentnode->GetEdgeLength())/(1.0*(breaksperbranch+1)); //adding breaksperbranch nodes of degree 2 divides the branch into breaksperbranch+1 segments
		gsl_matrix * Pmatrix=ComputeTransitionProb(RateMatrix,eachsegmentlength);
		if (currentnode->IsLeaf() ) {
			for (int breaknum=0;breaknum<(breaksperbranch+1);breaknum++) { //initialize the vectors
				vector<double> tempdouble;
				vector<int> tempint;
				for(int statepos=0;statepos<ancestralstatevector->size;statepos++) {
					tempdouble.push_back(0.0);
					tempint.push_back(0);
				}
				(Lvector[currentnode]).push_back(tempdouble);
				(Cvector[currentnode]).push_back(tempint);
			}
			for (int breaknum=0;breaknum<(breaksperbranch+1);breaknum++) { //Break nums are numbered from the tip down
				if (breaknum==0) { //means we're at the tip
					int statenumber=discretecharacters->GetInternalRepresentation(taxa->FindTaxon(currentnode->GetLabel()),discretechosenchar); //NOTE: for discrete chars, the number starts at 0
					for(int i=0;i<ancestralstatevector->size;i++) {
						(Cvector[currentnode])[breaknum][i]=statenumber; //we always must end up with the observed state
						(Lvector[currentnode])[breaknum][i]=(gsl_matrix_get(Pmatrix,i,statenumber));
					}
				}
				else { //at one of the degree two nodes on this "edge"
					for(int i=0;i<ancestralstatevector->size;i++) {
						double bestProb=0.0;
						int bestJ=-1;
						for(int j=0;j<ancestralstatevector->size;j++) {
							double  currentProb=(gsl_matrix_get(Pmatrix,i,j))*((Lvector[currentnode])[breaknum-1][j]); //Modify Pupko algorithm 2a: Lz(i)=maxj Pij(tz) x Lx(j)
							if (currentProb>bestProb) {
								bestProb=currentProb;
								bestJ=j;
							}
						}
						(Cvector[currentnode])[breaknum][i]=bestJ;
						(Lvector[currentnode])[breaknum][i]=bestProb;
					}
				}
			}
		}
		else if (currentnode!=(*Tptr).GetRoot()) { //must be an internal node, but not the root
			for (int breaknum=0;breaknum<(breaksperbranch+1);breaknum++) { //initialize the vectors
				vector<double> tempdouble;
				vector<int> tempint;
				for(int statepos=0;statepos<ancestralstatevector->size;statepos++) {
					tempdouble.push_back(0.0);
					tempint.push_back(0);
				}
				(Lvector[currentnode]).push_back(tempdouble);
				(Cvector[currentnode]).push_back(tempint);
			}
			for (int breaknum=0;breaknum<(breaksperbranch+1);breaknum++) { //Break nums are numbered from the tip down
				if (breaknum==0) { //means we're at the node of degree>2
					for(int i=0;i<ancestralstatevector->size;i++) {
						double bestProb=0.0;
						int bestJ=-1;
						for(int j=0;j<ancestralstatevector->size;j++) {
							double  currentProb=(gsl_matrix_get(Pmatrix,i,j));
							NodePtr descnode=currentnode->GetChild();
							while (descnode!=NULL) { 
								currentProb*=((Lvector[descnode])[breaksperbranch][j]); //Get the likelihood of state j at the earliest examined node on each of the descendant branches
								descnode=descnode->GetSibling(); //we're going to look at all descendant subtrees (even in case of polytomies)
							}
							if (currentProb>bestProb) {
								bestProb=currentProb;
								bestJ=j;
							}
						}
						(Cvector[currentnode])[breaknum][i]=bestJ;
						(Lvector[currentnode])[breaknum][i]=bestProb;
					}
				}
				else { //at one of the degree two nodes on this "edge"
					for(int i=0;i<ancestralstatevector->size;i++) {
						double bestProb=0.0;
						int bestJ=-1;
						for(int j=0;j<ancestralstatevector->size;j++) {
							double  currentProb=(gsl_matrix_get(Pmatrix,i,j))*((Lvector[currentnode])[breaknum-1][j]); //Modify Pupko algorithm 2a: Lz(i)=maxj Pij(tz) x Lx(j)
							if (currentProb>bestProb) {
								bestProb=currentProb;
								bestJ=j;
							}
						}
						(Cvector[currentnode])[breaknum][i]=bestJ;
						(Lvector[currentnode])[breaknum][i]=bestProb;
					}
				}
			}
		}
		else { //hooray! at root
			double bestProb=0.0;
			int bestK=-1;
			for(int k=0;k<ancestralstatevector->size;k++) {
				double currentProb=gsl_vector_get(ancestralstatevector,k);
				NodePtr descnode=currentnode->GetChild();
				while (descnode!=NULL) { 
					currentProb*=((Lvector[descnode])[breaksperbranch][k]); //Get the likelihood of state j at the earliest examined node on each of the descendant branches
					descnode=descnode->GetSibling(); //we're going to look at all descendant subtrees (even in case of polytomies)
				}
				if (currentProb>bestProb) {
					bestProb=currentProb;
					bestK=k;
				}
			}
			//message="LnLikelihood of this reconstruction is ";
			//message+=log(bestProb);
			//PrintMessage();
			nxsstring AncStateLabel="";
			AncStateLabel+=bestK;
			currentnode->SetLabel(AncStateLabel);
		}
		currentnode = n.next();
		gsl_matrix_free(Pmatrix);
	}
	
	//Now, back up the tree
	map<Node*, nxsstring> NewLabels;
	map<Node*, nxsstring> SimmapLabels;
	map<Node*, nxsstring> OriginalLabels;
	PreorderIterator <Node> m ((*Tptr).GetRoot()); //Goes from root up
	currentnode = m.begin();
	while (currentnode)
	{
		if (currentnode!=(*Tptr).GetRoot() ) {
			nxsstring newlabeltext="";
			nxsstring simmaplabeltext="";
			nxsstring oldlabeltext="";
			if (currentnode->IsLeaf()) {
				newlabeltext+=currentnode->GetLabel();
				simmaplabeltext+=currentnode->GetLabel();
				oldlabeltext+=currentnode->GetLabel();
			}
			newlabeltext+="[&S "; //Has list of states on subtending branch
			simmaplabeltext+=":{";
			double eachsegmentlength=(currentnode->GetEdgeLength())/(1.0*(breaksperbranch+1));
			double elapsedtime=0.0;
			int numchanges=0;
			vector<int> stateordervector; 
			vector<double> statetimesvector;
			vector<double> modelvector(maxModelCategoryStates,0.0);
			int i=-1;
			int j=-1;
			for (int breaknum=breaksperbranch;breaknum>=0;breaknum--) { //Break nums are numbered from the tip down, we want to go from the bottom up
				if (breaknum==breaksperbranch) { //at the rootmost interval
					i=atoi(((currentnode->GetAnc())->GetLabel()).c_str()); //The label has the last state
					if (debugmode) {
						(*Tptr).Update();
						(*Tptr).SetPathLengths();
						cout << endl << "Postorder traversal of tree" << endl;
						cout << "Node Length     Label                           Type" << endl;
						cout << "----------------------------------------------------" << endl;
						int count = 0;
						NodeIterator <Node> n ((*Tptr).GetRoot());
						Node *q = n.begin();
						Node *a=n.begin();
						Node *p=n.begin();
						Node *b=n.begin();
						Node *mrca=n.begin();
						nxsstring listofchanges="\n\nList of changes";
						bool notfirstleaf=false;
						while (q)
						{
							count++;
							cout << setiosflags(ios::right) << setw(4) << count					// arbitrary counter
								<< " " << setw(8) << q->GetEdgeLength() 							// edge length
								<< " " << setw(32) << setiosflags(ios::left) << q->GetLabel()	// node label
								;
							if (q->IsLeaf())
							{
								cout << " [LEAF]";
								if (q->IsMarked()) {
									cout <<" Node is marked\n\tModelCategorySize: ";
								}
								else {
									cout <<" Node not marked\n\tModelCategorySize: ";
								}
								vector<double> modelcatoutput(q->GetModelCategory());
								cout<<"( ";
								for (int i=0;i<modelcatoutput.size();i++) {
									cout<<modelcatoutput[i]<<" ";
								}
								cout<<") "<<endl;
								vector<int> stateordervector(q->GetStateOrder()); 
								vector<double> statetimesvector(q->GetStateTimes());
								cout<<"\tStateOrder: ( ";
								for (int i=0;i<stateordervector.size();i++) {
									cout<<stateordervector[i]<<" ";
									if (i>0) {
										listofchanges+="\nchange ";
										listofchanges+=stateordervector[i-1];
										listofchanges+=" -> ";
										listofchanges+=stateordervector[i];
										double timeofchange=(q->GetPathLength())-(q->GetEdgeLength());
										for (int j=0; j<i; j++) {
											timeofchange+=statetimesvector[j];
										}
										listofchanges+=" at time above root of ";
										listofchanges+=timeofchange;
									}
								}
								cout<<") "<<endl;			
								cout<<"\tStateTimes: ( ";
								for (int i=0;i<statetimesvector.size();i++) {
									cout<<statetimesvector[i]<<" ";
								}
								cout<<") "<<endl;			
								
								
								//gsl_vector modelcatoutput(q->GetModelCategory());
								//cout<<modelcatoutput->size;
								cout<<"\t length from root="<<q->GetPathLength();
								cout<<endl<<" length of edge="<<q->GetEdgeLength();
								float pathlength=0;
								a=q;
								if(notfirstleaf) {
									//  cout<<"Get MRCA of "<< q->GetLabel()<<" and "<<p->GetLabel()<<endl;
									bool mrcanotfound=true;
									while (a->GetAnc() && mrcanotfound) {
										a=a->GetAnc();
										b=p;
										while (b->GetAnc() && mrcanotfound) {
											b=b->GetAnc();
											if (a==b) {
												mrca=a;
												mrcanotfound=false;
												while(mrca->GetAnc()) {
													pathlength+=mrca->GetEdgeLength();
													//         cout<<"Edge length: "<<mrca->GetEdgeLength()<<" Total: "<<pathlength<<endl;
													mrca=mrca->GetAnc();
												}
											}
										}
									}
								}
								notfirstleaf=true;
								p=q;
								//  cout<<"Root to MRCA length="<<pathlength<<endl;
								//    while (aanc != Root && mrcanotfound) {
								//        banc=b;
								//        aanc=aanc->GetAnc();
								//       while (banc != Root && mrcanotfound) {
								//           banc=banc->GetAnc();
								//          if (aanc == banc) {
								//              mrcaptr=aanc;
								//               mrcanotfound=false;
								//           }
								//        }
								//    }
								//    return mrcaptr;
				
				
								//a=q;
								//cout<<endl<<"interim path length 1: "<<pathlength<<endl;
								//while (a->GetAnc()) {
								//   r=a->GetAnc();
								//  pathlength+=r->GetEdgeLength();
								// cout<<"interim: "<<pathlength<<" and edge "<<r->GetEdgeLength()<<endl;
								// a=r;
								//}
								//cout<<endl<<"Total length="<<pathlength<<endl;
							}
							else
							{
								cout << " [INTERNAL]";
								if (q->IsMarked()) {
									cout <<" Node is marked  ModelCategorySize: ";
								}
								else {
									cout <<" Node not marked  ModelCategorySize: ";
								}
								vector<double> modelcatoutput(q->GetModelCategory());
								cout<<"( ";
								for (int i=0;i<modelcatoutput.size();i++) {
									cout<<modelcatoutput[i]<<" ";
								}
								cout<<") "<<endl;
								vector<int> stateordervector(q->GetStateOrder()); 
								vector<double> statetimesvector(q->GetStateTimes());
								cout<<"\tStateOrder: ( ";
								for (int i=0;i<stateordervector.size();i++) {
									cout<<stateordervector[i]<<" ";
									if (i>0) {
										listofchanges+="\nchange ";
										listofchanges+=stateordervector[i-1];
										listofchanges+=" -> ";
										listofchanges+=stateordervector[i];
										double timeofchange=(q->GetPathLength())-(q->GetEdgeLength());
										for (int j=0; j<i; j++) {
											timeofchange+=statetimesvector[j];
										}
										listofchanges+=" at time above root of ";
										listofchanges+=timeofchange;
									}					
								}
								cout<<") "<<endl;			
								cout<<"\tStateTimes: ( ";
								for (int i=0;i<statetimesvector.size();i++) {
									cout<<statetimesvector[i]<<" ";
								}
								cout<<") "<<endl;	
								//gsl_vector *modelcatoutput;
								//modelcatoutput=gsl_vector_calloc(q->GetModelCategory());
								//gsl_vector modelcatoutput(q->GetModelCategory());
								//cout<<modelcatoutput->size;
								cout<<" length from root="<<q->GetPathLength();
								cout<<endl<<" length of edge="<<q->GetEdgeLength();
							}
							cout << endl;
							q = n.next();
						}
						cout << "----------------------------------------------------" << endl;
						message=listofchanges;
						PrintMessage();




						cout<<endl<<"currentnode "<<currentnode<<" currentnode->GetAnc() "<<currentnode->GetAnc()<<" ((currentnode->GetAnc())->GetLabel()) "<<((currentnode->GetAnc())->GetLabel())<<endl;
					}
				}
				j=(Cvector[currentnode])[breaknum][i];
				//stateordervector.push_back(j);
				//statetimesvector.push_back(eachsegmentlength);
				if (debugmode) {
					cout<<endl<<"breaksperbranch "<<breaksperbranch;
					cout<<endl<<"breaknum "<<breaknum;
					cout<<endl<<"i "<<i;
					cout<<endl<<"j "<<j;
					cout<<endl<<"eachsegmentlength "<<eachsegmentlength;
				}
				modelvector[j]+=eachsegmentlength;
				if (i!=j) { //There's been a change of state
					if (numchanges>0) {
						newlabeltext+=",";
						simmaplabeltext+=":";
					}
					newlabeltext+=i;
					newlabeltext+=":";
					newlabeltext+=elapsedtime+(0.5*eachsegmentlength);
					stateordervector.push_back(i);
					statetimesvector.push_back(elapsedtime+(0.5*eachsegmentlength));
					simmaplabeltext+=i;
					simmaplabeltext+=",";
					simmaplabeltext+=elapsedtime+(0.5*eachsegmentlength);
					elapsedtime=0.5*eachsegmentlength;
					numchanges++;
				}
				else {
					elapsedtime+=eachsegmentlength;
				}
				i=j;
				if (breaknum==0 && !(currentnode->IsLeaf())) { //change labels to include ancestral state
					nxsstring AncStateLabel="";
					AncStateLabel+=j;
					currentnode->SetLabel(AncStateLabel);
				}
				if (breaknum==0) {
					if (numchanges>0) {
						newlabeltext+=",";	
						simmaplabeltext+=":";
					}
					newlabeltext+=j;
					newlabeltext+=":";
					newlabeltext+=elapsedtime;
					newlabeltext+="]";
					stateordervector.push_back(j);
					statetimesvector.push_back(elapsedtime);
					simmaplabeltext+=j;
					simmaplabeltext+=",";
					simmaplabeltext+=elapsedtime;
					simmaplabeltext+="}";
				}
			}
			currentnode->SetStateOrder(stateordervector);
			currentnode->SetStateTimes(statetimesvector);
			currentnode->SetModelCategory(modelvector);
			if (debugmode) {
				cout<<"SetModelCategory(";
				for (int i=0;i<modelvector.size();i++) {
					cout<<" "<<modelvector[i];
				}
				cout<<" )"<<endl;
				vector<double> returnedmodelvector(currentnode->GetModelCategory());
				cout<<"GetModelCategory(";
				for (int i=0;i<returnedmodelvector.size();i++) {
					cout<<" "<<returnedmodelvector[i];
				}
				cout<<" )"<<endl;
				
				
			}
			NewLabels[currentnode]=newlabeltext;
			SimmapLabels[currentnode]=simmaplabeltext;
			OriginalLabels[currentnode]=oldlabeltext;
			if (debugmode) {
				cout<<"Node "<<currentnode<<" "<<currentnode->GetLabel()<<endl<<"\tState vectors: ";
				for (int vectorpos=0;vectorpos<stateordervector.size();vectorpos++) {
						cout<<stateordervector[vectorpos]<<":"<<statetimesvector[vectorpos]<<" ";
				}
				cout<<endl;
				for (int  vectorpos=0;vectorpos<modelvector.size();vectorpos++) {
					cout<<modelvector[vectorpos]<<" ";
				}
				cout<<endl;
				cout<<"NewLabel = "<<newlabeltext<<" NewLabels.size()="<<NewLabels.size()<<endl<<"SimmapLabel = "<<simmaplabeltext<<" SimmapLabels.size()="<<SimmapLabels.size()<<endl;
			}

		}
		currentnode = m.next();
	}
	(*Tptr).SetEdgeLengths(true);
//	OldFormat.SetEdgeLengths(true);
//	OldFormat.Update();
//	OldFormat.GetNodeDepths();
//	OldFormat.Draw(cout);
	(*Tptr).Update();
	(*Tptr).GetNodeDepths();
	//(*Tptr).Draw(cout);
	(*Tptr).SetInternalLabels(true);
	PreorderIterator <Node> q ((*Tptr).GetRoot()); //Goes from root up
	currentnode = q.begin();
	if(debugmode)  {
		cout<<" NewLabels.size()="<<NewLabels.size()<<endl<<" SimmapLabels.size()="<<SimmapLabels.size()<<endl;
	}
	while (currentnode)
	{
		if (currentnode!=(*Tptr).GetRoot()) {
			currentnode->SetLabel(NewLabels[currentnode]);
			if(debugmode) {
				cout<<"setting node "<<currentnode<<" to label "<<NewLabels[currentnode]<<endl;
				vector<double> returnedmodelvector(currentnode->GetModelCategory());
				cout<<"GetModelCategory(";
				for (int i=0;i<returnedmodelvector.size();i++) {
					cout<<" "<<returnedmodelvector[i];
				}
				cout<<" )"<<endl;
				
			}
		}
		currentnode = q.next();
	}
	((*Tptr).GetRoot())->SetLabel("");
	(*Tptr).Update();
	(*Tptr).GetNodeDepths();
	cout<<"tree input"<<chosentree<<" = ";
	(*Tptr).WriteNoQuote(cout);
	cout<<"\n";
	if( logf_open ) {
		logf<<"begin trees;\n[reconstruction of ancestral states]\ntree input"<<chosentree<<" = ";
		(*Tptr).WriteNoQuote(logf);
		logf<<"\nend;\n";
		
	}
	
	currentnode = q.begin();
	while (currentnode)
	{
		if (currentnode!=(*Tptr).GetRoot()) {
			currentnode->SetLabel(SimmapLabels[currentnode]);
			if(debugmode) {
				cout<<"setting node "<<currentnode<<" to label "<<SimmapLabels[currentnode]<<endl;
				vector<double> returnedmodelvector(currentnode->GetModelCategory());
				cout<<"GetModelCategory(";
				for (int i=0;i<returnedmodelvector.size();i++) {
					cout<<" "<<returnedmodelvector[i];
				}
				cout<<" )"<<endl;
				
			}
			
		}
		currentnode = q.next();
	}
	(*Tptr).Update();
	(*Tptr).GetNodeDepths();
	(*Tptr).SetEdgeLengths(false); //So it won't print out the edge lengths;
	cout<<"tree input"<<chosentree<<" = ";
	(*Tptr).WriteNoQuote(cout);
	cout<<endl;
	(*Tptr).Draw(cout);
	cout<<"\n";
	if( logf_open ) {
		logf<<"begin trees;\n[reconstruction of ancestral states]\ntree input"<<chosentree<<" = ";
		(*Tptr).WriteNoQuote(logf);
		logf<<"\nend;\n";
		
	}
	//Now go back to original form:
	currentnode = q.begin();
	while (currentnode)
	{
		if (currentnode!=(*Tptr).GetRoot()) {
			currentnode->SetLabel(OriginalLabels[currentnode]);
			if(debugmode) {
				cout<<"setting node "<<currentnode<<" to label "<<SimmapLabels[currentnode]<<endl;
				vector<double> returnedmodelvector(currentnode->GetModelCategory());
				cout<<"GetModelCategory(";
				for (int i=0;i<returnedmodelvector.size();i++) {
					cout<<" "<<returnedmodelvector[i];
				}
				cout<<" )"<<endl;				
			}
			
		}
		currentnode = q.next();
	}
	
	(*Tptr).Update();
	(*Tptr).GetNodeDepths();
	(*Tptr).SetEdgeLengths(true); //Since it has edge lengths;
	//(*Tptr).Draw(cout);
	
	
	/*ContainingTree A;
	A.SetRoot((*Tptr).CopyOfSubtree((*Tptr).GetRoot()));
	A.Update();
	A.GetNodeDepths();
	A.ReportTreeHealth();
*/

	
	return (*Tptr).GetRoot();	
}

/*  ORIGINAL FUNCTION FOR CALCULATING THE LIKELIHOOD: HAS ISSUES WITH REALLY TINY LIKELIHOODS, AS IT DOESN'T USE LOGS ON THE DOWNPASS
//Calculates the likelihood of discrete character discretechosenchar on tree chosentree
double BROWNIE::CalculateDiscreteCharLnL(gsl_matrix * RateMatrix, gsl_vector * ancestralstatevector)
{
	double neglnL=0;
	double Prob=0;
	if (variablecharonly) {
		Prob=CalculateDiscreteCharProbAllConstant(RateMatrix,ancestralstatevector);
	}			
	int startingdiscretechosenchar=discretechosenchar;
	int endingdiscretechosenchar=discretechosenchar+1;
	if (allchar) {
		startingdiscretechosenchar=0;
		endingdiscretechosenchar=discretecharacters->GetNChar();
	}
	int olddiscretechosenchar=discretechosenchar;
	for (discretechosenchar=startingdiscretechosenchar;discretechosenchar<endingdiscretechosenchar;discretechosenchar++) {
		if ((discretecharacters->GetObsNumStates(discretechosenchar))>1 || variablecharonly==false) { //so, ignore invariant characters if variablecharonly==true
			long double L=0;
			map<Node*, vector<long double> > stateprobatnodes;
			Tree T=intrees.GetIthTree(chosentree-1);
			NodeIterator <Node> n (T.GetRoot()); //Goes from tips down
			NodePtr currentnode = n.begin();
			while (currentnode)
			{
				if (currentnode->IsLeaf() ) {
					int statenumber=discretecharacters->GetInternalRepresentation(taxa->FindTaxon(currentnode->GetLabel()),discretechosenchar); //NOTE: for discrete chars, the number starts at 0
					for(int j=0;j<ancestralstatevector->size;j++) {
						long double probofstatej=0; //do all in straight prob, then convert to ln L
						if (j==statenumber) {
							probofstatej=1; 
						}
						(stateprobatnodes[currentnode]).push_back(probofstatej);
					}
				}
				else { //must be an internal node, including the root
					for(int i=0;i<ancestralstatevector->size;i++) { //do this for each possible state at the current node
						NodePtr descnode=currentnode->GetChild();
						long double probofstatei=1;
						while (descnode!=NULL) { 
							gsl_matrix * Pmatrix=ComputeTransitionProb(RateMatrix,descnode->GetEdgeLength());
							long double probofthissubtree=0;
							for(int j=0;j<ancestralstatevector->size;j++) {
								probofthissubtree+=(gsl_matrix_get(Pmatrix,i,j))*((stateprobatnodes[descnode])[j]); //Prob of going from i to j on desc branch times the prob of the subtree with root state j
							}
							probofstatei*=probofthissubtree;
							descnode=descnode->GetSibling(); //we're going to look at all descendant subtrees (even in case of polytomies)
							gsl_matrix_free(Pmatrix);
						}
						(stateprobatnodes[currentnode]).push_back(probofstatei);
						if (debugmode) {
							cout<<"CalculateDiscreteCharLnL: i = "<<i<<", probofstatei = "<<probofstatei<<", ln(probofstatei) = "<<log(probofstatei)<<endl;
						}
					}
				}
				currentnode = n.next();
				
			}
	//now, finish up by getting the weighted sum at the root
			for (int i=0;i<ancestralstatevector->size;i++) {
				L+=(gsl_vector_get(ancestralstatevector,i))*((stateprobatnodes[T.GetRoot()])[i]);
			}
			if (variablecharonly) {
				L=L/(1.0-Prob); //after equation 3 in Lewis 2001 and equation 8 in Felsenstein 1992
			}		
			if (debugmode) {
				cout<<"CalculateDiscreteCharLnL: original ("<<neglnL<<" * -1.0*log(L) [L="<<L<<"] = ";
			}
			neglnL+=-1.0*log(L);
			if (debugmode) {
				cout<<neglnL<<endl;
			}
			
		}
	}
	discretechosenchar=olddiscretechosenchar;
	return neglnL;
}
*/

/* A failed attempt at a solution
//Calculates the likelihood of discrete character discretechosenchar on tree chosentree
//Uses logs on the down pass and some crude approaches to prevent underflows
double BROWNIE::CalculateDiscreteCharLnL(gsl_matrix * RateMatrix, gsl_vector * ancestralstatevector)
{
	double neglnL=0;
	double Prob=0;
	if (variablecharonly) {
		Prob=CalculateDiscreteCharProbAllConstant(RateMatrix,ancestralstatevector);
	}			
	int startingdiscretechosenchar=discretechosenchar;
	int endingdiscretechosenchar=discretechosenchar+1;
	if (allchar) {
		startingdiscretechosenchar=0;
		endingdiscretechosenchar=discretecharacters->GetNChar();
	}
	int olddiscretechosenchar=discretechosenchar;
	for (discretechosenchar=startingdiscretechosenchar;discretechosenchar<endingdiscretechosenchar;discretechosenchar++) {
		if ((discretecharacters->GetObsNumStates(discretechosenchar))>1 || variablecharonly==false) { //so, ignore invariant characters if variablecharonly==true
			 double L=0;
			map<Node*, vector< double> > stateprobatnodes; //actually, -ln(probability)
			Tree T=intrees.GetIthTree(chosentree-1);
			NodeIterator <Node> n (T.GetRoot()); //Goes from tips down
			NodePtr currentnode = n.begin();
			while (currentnode)
			{
				if (currentnode->IsLeaf() ) {
					int statenumber=discretecharacters->GetInternalRepresentation(taxa->FindTaxon(currentnode->GetLabel()),discretechosenchar); //NOTE: for discrete chars, the number starts at 0
					for(int j=0;j<ancestralstatevector->size;j++) {
						 double probofstatej=-log(0); //do all in -ln(prob)
						if (j==statenumber) {
							probofstatej=-log(1.0); 
						}
						(stateprobatnodes[currentnode]).push_back(probofstatej);
						if (debugmode) {
							cout<<"CalculateDiscreteCharLnL: j = "<<j<<", probofstatej = "<<probofstatej<<", exp(-probofstatej) = "<<exp(-probofstatej)<<endl;
						}
						
					}
				}
				else { //must be an internal node, including the root
					for(int i=0;i<ancestralstatevector->size;i++) { //do this for each possible state at the current node
						NodePtr descnode=currentnode->GetChild();
						double probofstatei=0;
						while (descnode!=NULL) { 
							gsl_matrix * Pmatrix=ComputeTransitionProb(RateMatrix,descnode->GetEdgeLength()); //this is in terms of probabilities, not log probabilities
							if (debugmode) {
								PrintMatrix(Pmatrix);
							}
							double probofthissubtree=0;
							double minimumdescendantstateprob=GSL_POSINF; //This is used to rescale the probabilities
							for (int k=0; k<ancestralstatevector->size; k++) {
								minimumdescendantstateprob=GSL_MIN((stateprobatnodes[descnode])[k],minimumdescendantstateprob);
							}
							for(int j=0;j<ancestralstatevector->size;j++) {
								 //here, for ancestor with state i and descendant with state j, we want to calculate the probability of going from i to j multiplied by the probability of subtree j
								//that is, P(i,j) * P(subtree j)
								//and then add those for each possible j to get P(subtree i).
								//The problem is that P(subtree j) can get too small for even long doubles.
								//The solution is to calculate 
								//  ( P(i,j1)*P(subtree j1)*scaling_factor + P(i,j2)*P(subtree j2)*scaling_factor + P(i,j3)*P(subtree j3)*scaling_factor + ...) / scaling_factor
								// and then take the negative log of this for storing back as the -ln(probability) of subtree i
								// (stateprobatnodes[descnode])[j] = -ln(P(subtree j1))
								// P(subtree j1) = exp(-(stateprobatnodes[descnode])[j])
								// P(subtree j1) / scaling_factor = exp(-(stateprobatnodes[descnode])[j]) / scaling_factor
								// P(subtree j1) / scaling_factor = exp(-(stateprobatnodes[descnode])[j]) / exp(different_factor)
								// P(subtree j1) / scaling_factor = exp( -(stateprobatnodes[descnode])[j] - different_factor)
								//for different factor, use different_factor=min(stateprobatnodes[descnode]) 
								//[could have used max, but -ln(0) = +inf, and we'd see these for things like tips, where prababilities of some states is 0]
								
								
								probofthissubtree+=(gsl_matrix_get(Pmatrix,i,j))*exp(-1.0*((stateprobatnodes[descnode])[j]-minimumdescendantstateprob)); //Prob of going from i to j on desc branch times the prob of the subtree with root state j
								if (debugmode) {
									cout<<"going from "<<i<<" to "<<j<<" has move prob "<<gsl_matrix_get(Pmatrix,i,j)<<" descendant prob "<<exp(-1.0*((stateprobatnodes[descnode])[j]))<<" total prob "<<(gsl_matrix_get(Pmatrix,i,j))*exp(-1.0*((stateprobatnodes[descnode])[j]))<<" and rescaled prob "<<(gsl_matrix_get(Pmatrix,i,j))*exp(-1.0*((stateprobatnodes[descnode])[j]-minimumdescendantstateprob))<<" with scaling factor "<<minimumdescendantstateprob<<endl;
								}
							}
							//Now take -ln( P(i,j1)*P(subtree j1)*scaling_factor + P(i,j2)*P(subtree j2)*scaling_factor + P(i,j3)*P(subtree j3)*scaling_factor + ...) / scaling_factor)
							// = -ln (probofthissubtree / scaling_factor)
							// = -ln (probofthisubtree) - (-ln(scaling_factor))
							// = -ln (probofthisubtree) - different_factor
							probofstatei+=-(log(probofthissubtree))-minimumdescendantstateprob;
							descnode=descnode->GetSibling(); //we're going to look at all descendant subtrees (even in case of polytomies)
							gsl_matrix_free(Pmatrix);
						}
						(stateprobatnodes[currentnode]).push_back(probofstatei);
						if (debugmode) {
							cout<<"CalculateDiscreteCharLnL: i = "<<i<<", probofstatei = "<<probofstatei<<", exp(-probofstatei) = "<<exp(-probofstatei)<<endl;
						}
					}
				}
				currentnode = n.next();
				
			}
	//now, finish up by getting the weighted sum at the root
			double minimumdescendantstateprob=GSL_POSINF; //This is used to rescale the probabilities
			for (int k=0; k<ancestralstatevector->size; k++) {
				minimumdescendantstateprob=GSL_MIN((stateprobatnodes[T.GetRoot()])[k],minimumdescendantstateprob);
			}
			
			for (int i=0;i<ancestralstatevector->size;i++) {
				L+=(gsl_vector_get(ancestralstatevector,i))*exp(-1.0*((stateprobatnodes[T.GetRoot()])[i]-minimumdescendantstateprob));
			}
			if (variablecharonly) {
				L=L/(1.0-Prob); //after equation 3 in Lewis 2001 and equation 8 in Felsenstein 1992
			}		

			if (debugmode) {
				cout<<"CalculateDiscreteCharLnL: original ("<<neglnL<<" * -1.0*log(L) [L="<<L<<"] = ";
			}
			neglnL+=-(log(L))-minimumdescendantstateprob;
			if (debugmode) {
				cout<<neglnL<<endl;
			}
			
		}
	}
	discretechosenchar=olddiscretechosenchar;
	return neglnL;
}
*/


//Calculates the likelihood of discrete character discretechosenchar on tree chosentree
//Deals with underflow issues by using superdouble
double BROWNIE::CalculateDiscreteCharLnL(gsl_matrix * RateMatrix, gsl_vector * ancestralstatevector)
{
	double neglnL=0;
	double Prob=0;
	if (variablecharonly) {
		Prob=CalculateDiscreteCharProbAllConstant(RateMatrix,ancestralstatevector);
	}			
	int startingdiscretechosenchar=discretechosenchar;
	int endingdiscretechosenchar=discretechosenchar+1;
	if (allchar) {
		startingdiscretechosenchar=0;
		endingdiscretechosenchar=discretecharacters->GetNChar();
	}
	int olddiscretechosenchar=discretechosenchar;
	for (discretechosenchar=startingdiscretechosenchar;discretechosenchar<endingdiscretechosenchar;discretechosenchar++) {
		if ((discretecharacters->GetObsNumStates(discretechosenchar))>1 || variablecharonly==false) { //so, ignore invariant characters if variablecharonly==true
			Superdouble L=0;
			map<Node*, vector<Superdouble> > stateprobatnodes;
			Tree T=intrees.GetIthTree(chosentree-1);
			NodeIterator <Node> n (T.GetRoot()); //Goes from tips down
			NodePtr currentnode = n.begin();
			while (currentnode)
			{
				if (currentnode->IsLeaf() ) {
					int statenumber=discretecharacters->GetInternalRepresentation(taxa->FindTaxon(currentnode->GetLabel()),discretechosenchar); //NOTE: for discrete chars, the number starts at 0
					for(int j=0;j<ancestralstatevector->size;j++) {
						Superdouble probofstatej=0; //do all in straight prob, then convert to ln L
						if (j==statenumber) {
							probofstatej=1; 
						}
						(stateprobatnodes[currentnode]).push_back(probofstatej);
/*						if (debugmode) {
							cout<<"CalculateDiscreteCharLnL: j = "<<j<<", probofstatej = "<<probofstatej.getMantissa()<<" x 10^"<<probofstatej.getExponent()<<", -ln(probofstatej) = "<<-1.0*probofstatej.getLn()<<endl<<endl;
						}*/
						
					}
					
				}
				else { //must be an internal node, including the root
					for(int i=0;i<ancestralstatevector->size;i++) { //do this for each possible state at the current node
						NodePtr descnode=currentnode->GetChild();
						Superdouble probofstatei=1;
						while (descnode!=NULL) { 
							gsl_matrix * Pmatrix=ComputeTransitionProb(RateMatrix,descnode->GetEdgeLength());
							Superdouble probofthissubtree=0;
/*							if (debugmode) {
								Tree PrunedTree;
								PrunedTree.SetRoot((intrees.GetIthTree(chosentree-1)).CopyOfSubtree(descnode));
								cout<<endl;
								PrunedTree.Write(cout);
							}*/
							for(int j=0;j<ancestralstatevector->size;j++) {
								Superdouble transitionprob=gsl_matrix_get(Pmatrix,i,j);
								probofthissubtree+=transitionprob*((stateprobatnodes[descnode])[j]); //Prob of going from i to j on desc branch times the prob of the subtree with root state j
							/*	if (debugmode) {
									cout<<"From "<<i<<" to "<<j<<" has move prob "<<gsl_matrix_get(Pmatrix,i,j)<<" and input subtree prob "<<((stateprobatnodes[descnode])[j])<<" product "<<transitionprob*((stateprobatnodes[descnode])[j])<<" and cumulative prob "<<probofthissubtree<<endl;
								} */
							}
							probofstatei*=probofthissubtree;
							descnode=descnode->GetSibling(); //we're going to look at all descendant subtrees (even in case of polytomies)
							gsl_matrix_free(Pmatrix);
						}
						(stateprobatnodes[currentnode]).push_back(probofstatei);
				/*		if (debugmode) {
							cout<<"CalculateDiscreteCharLnL: i = "<<i<<", probofstatei = "<<probofstatei<<", -ln(probofstatei) = "<<-1.0*probofstatei.getLn()<<endl<<endl;
						} */
					}
				}
		/*		if (debugmode) {
					Superdouble totalprob=0;
					cout<<"(stateprobatnodes[currentnode])[i] ";
					for(int i=0;i<ancestralstatevector->size;i++) {
						cout<<"state "<<i<<" "<<(stateprobatnodes[currentnode])[i]<<" [totalprob is "<<totalprob<<"], ";
						totalprob+=(stateprobatnodes[currentnode])[i];
					}
					cout<<" mantissa="<<totalprob.getMantissa()<<"total = "<<totalprob<<endl;
					//assert(totalprob.getMantissa()>0);
				} */
				currentnode = n.next();
				
			}
	//now, finish up by getting the weighted sum at the root
			for (int i=0;i<ancestralstatevector->size;i++) {
				Superdouble ancestralprob=gsl_vector_get(ancestralstatevector,i);
				L+=ancestralprob*((stateprobatnodes[T.GetRoot()])[i]);
			}
			if (variablecharonly) {
				L=L/Superdouble(1.0-Prob); //after equation 3 in Lewis 2001 and equation 8 in Felsenstein 1992
			}		
/*			if (debugmode) {
				cout<<"CalculateDiscreteCharLnL: original ("<<neglnL<<" * -1.0*log(L) [L="<<L<<"] = ";
			}  */
			neglnL+=-1.0*L.getLn();
/*			if (debugmode) {
				cout<<neglnL<<endl;
			} */
			
		}
	}
	discretechosenchar=olddiscretechosenchar;
	if (1==isnan(neglnL)) { //this is not a number, which makes optimization difficult
		if (debugmode) {
			message="\nWarning: The negative ln likelihood in CalculateDiscreteCharLnL was NaN, so a very large value (BROWNIE_MAXLIKELIHOOD) was returned instead.\n";
			PrintMessage();
		}
		neglnL=BROWNIE_MAXLIKELIHOOD ;
	}
	else if (neglnL<0) {
		if (debugmode) {
			message="\nWarning: The negative ln likelihood in CalculateDiscreteCharLnL was less than zero (";
			message+=neglnL;
			message+="), so a very large value (BROWNIE_MAXLIKELIHOOD) was returned instead.\n";
			PrintMessage();
		}
		neglnL=BROWNIE_MAXLIKELIHOOD ;

	}
	return neglnL;
}

//Calculates the likelihood of getting only constant characters (all 0, or all 1, or all...)
double BROWNIE::CalculateDiscreteCharProbAllConstant(gsl_matrix * RateMatrix, gsl_vector * ancestralstatevector)
{
	double Prob=0;
	for (int tipstate=0;  tipstate<localnumbercharstates; tipstate++) { //we loop over all possible tip states
		double L=0;
		map<Node*, vector<double> > stateprobatnodes;
		Tree T=intrees.GetIthTree(chosentree-1);
		NodeIterator <Node> n (T.GetRoot()); //Goes from tips down
		NodePtr currentnode = n.begin();
		while (currentnode)
		{
			if (currentnode->IsLeaf() ) {
				int statenumber=tipstate; //we force all tips to have the same state
				for(int j=0;j<ancestralstatevector->size;j++) {
					double probofstatej=0; //do all in straight prob, then convert to ln L
					if (j==statenumber) {
						probofstatej=1; 
					}
					(stateprobatnodes[currentnode]).push_back(probofstatej);
				}
			}
			else { //must be an internal node, including the root
				for(int i=0;i<ancestralstatevector->size;i++) { //do this for each possible state at the current node
					NodePtr descnode=currentnode->GetChild();
					double probofstatei=1;
					while (descnode!=NULL) { 
						gsl_matrix * Pmatrix=ComputeTransitionProb(RateMatrix,descnode->GetEdgeLength());
						double probofthissubtree=0;
						for(int j=0;j<ancestralstatevector->size;j++) {
							probofthissubtree+=(gsl_matrix_get(Pmatrix,i,j))*((stateprobatnodes[descnode])[j]); //Prob of going from i to j on desc branch times the prob of the subtree with root state j
						}
						probofstatei*=probofthissubtree;
						descnode=descnode->GetSibling(); //we're going to look at all descendant subtrees (even in case of polytomies)
						gsl_matrix_free(Pmatrix);
					}
					(stateprobatnodes[currentnode]).push_back(probofstatei);
				}
			}
			currentnode = n.next();
		}
	//now, finish up by getting the weighted sum at the root
		for (int i=0;i<ancestralstatevector->size;i++) {
			L+=(gsl_vector_get(ancestralstatevector,i))*((stateprobatnodes[T.GetRoot()])[i]);
		}
		Prob+=L;
	}
	return Prob;
}

double BROWNIE::CalculateDiscreteLindy1(double rateA)
{
	double neglnL=CalculateDiscreteLindy2(rateA, rateA);
	return neglnL;
}

//Calculates the likelihood of going from state 1 to state 0 on one tree across all chars
double BROWNIE::CalculateDiscreteLindy2(double rateA, double rateB)
{
//	gsl_matrix * testmat = gsl_matrix_calloc(2,2);
//	gsl_matrix_set(testmat,0,0,-0.5);
//	gsl_matrix_set(testmat,0,1,0.5);
//	gsl_matrix_set(testmat,1,0,0.7);
//	gsl_matrix_set(testmat,1,1,-0.7);
//	gsl_matrix * testmatout=ComputeTransitionProb(testmat,1.0);
//	cout<<"testmatnow with brlen 1\n"<<gsl_matrix_get(testmatout,0,0)<<"\t"<<gsl_matrix_get(testmatout,0,1)<<"\n"<<gsl_matrix_get(testmatout,1,0)<<"\t"<<gsl_matrix_get(testmatout,1,1)<<endl;
//	testmatout=ComputeTransitionProb(testmat,0.0);
//	cout<<"testmatnow with brlen 0\n"<<gsl_matrix_get(testmatout,0,0)<<"\t"<<gsl_matrix_get(testmatout,0,1)<<"\n"<<gsl_matrix_get(testmatout,1,0)<<"\t"<<gsl_matrix_get(testmatout,1,1)<<endl;

	
	double neglnL=0.0;
	if (rateA<0.0 || rateB<0.0) {
		neglnL=GSL_POSINF;
	}
	else {
//		cout<<"rateA is "<<rateA<<endl;
		gsl_vector* ancestralstatevector=gsl_vector_calloc(2);
		gsl_vector_set(ancestralstatevector,0,0.0); //state 1 at root
		gsl_vector_set(ancestralstatevector,1,1.0);
		gsl_matrix* ratematrixA=gsl_matrix_calloc(2,2);
		//cout<<"Qmatrix starts as\n"<<gsl_matrix_get(ratematrixA,0,0)<<"\t"<<gsl_matrix_get(ratematrixA,0,1)<<"\n"<<gsl_matrix_get(ratematrixA,1,0)<<"\t"<<gsl_matrix_get(ratematrixA,1,1)<<endl;
		gsl_matrix_set(ratematrixA,1,0,rateA);
		//cout<<"Qmatrix now is is\n"<<gsl_matrix_get(ratematrixA,0,0)<<"\t"<<gsl_matrix_get(ratematrixA,0,1)<<"\n"<<gsl_matrix_get(ratematrixA,1,0)<<"\t"<<gsl_matrix_get(ratematrixA,1,1)<<endl;
		gsl_matrix_set(ratematrixA,1,1,0.0-rateA);
//		cout<<"ratematrixA now is is\n"<<gsl_matrix_get(ratematrixA,0,0)<<"\t"<<gsl_matrix_get(ratematrixA,0,1)<<"\n"<<gsl_matrix_get(ratematrixA,1,0)<<"\t"<<gsl_matrix_get(ratematrixA,1,1)<<endl;
		gsl_matrix* ratematrixB=gsl_matrix_calloc(2,2);
		gsl_matrix_set(ratematrixB,1,0,rateB);
		gsl_matrix_set(ratematrixB,1,1,0.0-rateB);
//		cout<<"ratematrixB now is is\n"<<gsl_matrix_get(ratematrixB,0,0)<<"\t"<<gsl_matrix_get(ratematrixB,0,1)<<"\n"<<gsl_matrix_get(ratematrixB,1,0)<<"\t"<<gsl_matrix_get(ratematrixB,1,1)<<endl;
		int nchartotal=discretecharacters->GetNChar();
		for (int currentchar=1;currentchar<=nchartotal;currentchar++) {
			//cout<<"rateA is "<<rateA<<endl;
			map<Node*, vector<double> > stateprobatnodes;
			Tree T=intrees.GetIthTree(chosentree-1);
			PreorderIterator <Node> n (T.GetRoot());
			NodePtr currentnode = n.begin();
			vector<double> ancstatevec;
			for (int i=0;i<ancestralstatevector->size;i++) {
				ancstatevec.push_back(gsl_vector_get(ancestralstatevector,i));
//				cout<<"anc state "<<gsl_vector_get(ancestralstatevector,i)<<endl;
			}
			stateprobatnodes[T.GetRoot()]=ancstatevec;
			while (currentnode)
			{
				if (currentnode!=T.GetRoot() && !currentnode->IsLeaf()) {
					gsl_matrix * Pmatrix=gsl_matrix_calloc(2,2);
					if ((currentnode->GetStateOrder())[0]==0) {
//						cout<<"\n\nQmatrix is\n"<<gsl_matrix_get(ratematrixA,0,0)<<"\t"<<gsl_matrix_get(ratematrixA,0,1)<<"\n"<<gsl_matrix_get(ratematrixA,1,0)<<"\t"<<gsl_matrix_get(ratematrixA,1,1)<<endl;
						Pmatrix=ComputeTransitionProb(ratematrixA,currentnode->GetEdgeLength());
//						cout<<"Qmatrix is\n"<<gsl_matrix_get(ratematrixA,0,0)<<"\t"<<gsl_matrix_get(ratematrixA,0,1)<<"\n"<<gsl_matrix_get(ratematrixA,1,0)<<"\t"<<gsl_matrix_get(ratematrixA,1,1)<<endl;

					}
					else {
						Pmatrix=ComputeTransitionProb(ratematrixB,currentnode->GetEdgeLength());
					}
//					cout<<"edgelength was "<<currentnode->GetEdgeLength()<<endl;
//					cout<<"Pmatrix is\n"<<gsl_matrix_get(Pmatrix,0,0)<<"\t"<<gsl_matrix_get(Pmatrix,0,1)<<"\n"<<gsl_matrix_get(Pmatrix,1,0)<<"\t"<<gsl_matrix_get(Pmatrix,1,1)<<endl;
					vector<double> statesatthisnode;
					for(int j=0;j<ancestralstatevector->size;j++) {
						double probofstatej=0;
						vector<double> statesatpreviousnode=stateprobatnodes[currentnode->GetAnc()];
						for (int i=0;i<ancestralstatevector->size;i++) {
//							cout<<"for j="<<j<<" prob of state "<<i<<" at prev node = "<<statesatpreviousnode[i]<<" and Pij="<<gsl_matrix_get(Pmatrix,i,j)<<endl;
							probofstatej+=(statesatpreviousnode[i])*gsl_matrix_get(Pmatrix,i,j);
						}
						statesatthisnode.push_back(probofstatej);
//						cout<<"j="<<j<<" prob of state j="<<probofstatej<<endl;
					}
					stateprobatnodes[currentnode]=statesatthisnode;
					gsl_matrix_free(Pmatrix);
				}
				else if (currentnode!=T.GetRoot() && currentnode->IsLeaf()) {
					int statenumber=discretecharacters->GetInternalRepresentation(taxa->FindTaxon(currentnode->GetLabel()),currentchar-1);
					gsl_matrix * Pmatrix=gsl_matrix_calloc(2,2);
					if ((currentnode->GetStateOrder())[0]==0) {
						Pmatrix=ComputeTransitionProb(ratematrixA,currentnode->GetEdgeLength());
					}
					else {
						Pmatrix=ComputeTransitionProb(ratematrixB,currentnode->GetEdgeLength());
					}
					for(int j=0;j<ancestralstatevector->size;j++) {
						double probofstatej=0;
						if (j==statenumber) {
							vector<double> statesatpreviousnode=stateprobatnodes[currentnode->GetAnc()];
							for (int i=0;i<ancestralstatevector->size;i++) {
								probofstatej+=(statesatpreviousnode[i])*gsl_matrix_get(Pmatrix,i,j);
							}
							neglnL+=(-1.0)*log(probofstatej);
							//cout<<"For tip with state "<<statenumber<<" probofstate="<<probofstatej<<" prob at previousnode="<<statesatpreviousnode[statenumber]<<" -lnL total so far = "<<neglnL<<endl;
						} 
					}
					gsl_matrix_free(Pmatrix);
				}
				currentnode = n.next();
			}
			//PreorderIterator <Node> m (T.GetRoot());
			//NodePtr currentnode2 = m.begin();
			//message="\nAncestralStates\n0\t1\n";
			//PrintMessage();
			//while (currentnode2) {
			//	message="";
			//	message+=(stateprobatnodes[currentnode])[0];
			//	message+="\t";
			//	message+=(stateprobatnodes[currentnode])[1];
			//	PrintMessage();
			//	currentnode2=m.next();
			//}
			

		}
	}
	return neglnL;
}

//uses pagel 1994's formula, P(t)=exp(Qt)+c=C*exp(Dt)*C^-1, where C is eigenvectors of Q and D is eigenvalues (in diagonal matrix)
gsl_matrix * BROWNIE::ComputeTransitionProb(gsl_matrix *RateMatrix, double brlen) {
	int dimension=RateMatrix->size1;
	gsl_matrix *ratematrixTMP=gsl_matrix_calloc(dimension,dimension);
	gsl_matrix_memcpy (ratematrixTMP, RateMatrix);
//	cout<<"ratematrixTMP is\n"<<gsl_matrix_get(ratematrixTMP,0,0)<<"\t"<<gsl_matrix_get(ratematrixTMP,0,1)<<"\n"<<gsl_matrix_get(ratematrixTMP,1,0)<<"\t"<<gsl_matrix_get(ratematrixTMP,1,1)<<endl;
	gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (dimension);
	gsl_vector_complex *eigenvalues=gsl_vector_complex_calloc(dimension);
	gsl_matrix_complex *eigenvectors=gsl_matrix_complex_calloc(dimension,dimension); //Matrix C as in Pagel 1994
	gsl_matrix *transitionmatrix=gsl_matrix_calloc(dimension,dimension);
	gsl_matrix_complex *transitionmatrixcomplex=gsl_matrix_complex_calloc(dimension,dimension);

	gsl_eigen_nonsymmv (ratematrixTMP, eigenvalues, eigenvectors, w);
//	cout<<"eigenvalues are "<<GSL_REAL(gsl_vector_complex_get(eigenvalues,0))<<" "<<GSL_REAL(gsl_vector_complex_get(eigenvalues,1))<<endl;
	gsl_matrix_complex *eigenvaluematrix=gsl_matrix_complex_calloc(dimension,dimension);
	//gsl_matrix *eigenvectorsreal=gsl_matrix_calloc(dimension,dimension);
	for (int i=0;i<dimension;i++) {
		gsl_matrix_complex_set(eigenvaluematrix,i,i,gsl_complex_exp(gsl_complex_mul_real(gsl_vector_complex_get(eigenvalues,i),brlen)));
	//	for (int j=0;j<dimension;j++) {
		//	gsl_matrix_set(eigenvectorsreal,i,j,GSL_REAL(gsl_matrix_complex_get(eigenvectors,i,j)));
	//	}
	}
	
	
	gsl_matrix_complex *inverseeigenvectorsstart=gsl_matrix_complex_calloc(dimension,dimension);
	gsl_matrix_complex *inverseeigenvectors=gsl_matrix_complex_calloc(dimension,dimension);
//	cout<<"eigenvaluematrix (after taking exp(t*lambda): \n"<<GSL_REAL(gsl_matrix_complex_get(eigenvaluematrix,0,0))<<" "<<GSL_REAL(gsl_matrix_complex_get(eigenvaluematrix,0,1))<<endl<<GSL_REAL(gsl_matrix_complex_get(eigenvaluematrix,1,0))<<"\t"<<GSL_REAL(gsl_matrix_complex_get(eigenvaluematrix,1,1))<<endl;

//	cout<<"eigenvectors: \n"<<GSL_REAL(gsl_matrix_complex_get(eigenvectors,0,0))<<" "<<GSL_REAL(gsl_matrix_complex_get(eigenvectors,0,1))<<endl<<GSL_REAL(gsl_matrix_complex_get(eigenvectors,1,0))<<"\t"<<GSL_REAL(gsl_matrix_complex_get(eigenvectors,1,1))<<endl;
	gsl_permutation * p = gsl_permutation_alloc (dimension);
	int signum;
	gsl_matrix_complex_memcpy (inverseeigenvectorsstart, eigenvectors);
	gsl_linalg_complex_LU_decomp (inverseeigenvectorsstart,p, &signum);
//	cout<<"inverseeigenvectorsstart: \n"<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectorsstart,0,0))<<" "<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectorsstart,0,1))<<endl<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectorsstart,1,0))<<"\t"<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectorsstart,1,1))<<endl;
	gsl_linalg_complex_LU_invert (inverseeigenvectorsstart,p, inverseeigenvectors);
//	cout<<"inverseeigenvectors: \n"<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectors,0,0))<<" "<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectors,0,1))<<endl<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectors,1,0))<<"\t"<<GSL_REAL(gsl_matrix_complex_get(inverseeigenvectors,1,1))<<endl;

	//cout<<"inverseeigenvectors: "<<gsl_matrix_get(inverseeigenvectors,0,0)<<" "<<gsl_matrix_get(inverseeigenvectors,0,1)<<endl<<gsl_matrix_get(inverseeigenvectors,1,0)<<"\t"<<gsl_matrix_get(inverseeigenvectors,1,1)<<endl;

	gsl_matrix_complex *stepA = gsl_matrix_complex_calloc(dimension,dimension);
	//gsl_complex one=gsl_complex_rect (1.0, 1.0);
	gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, eigenvaluematrix, inverseeigenvectors,GSL_COMPLEX_ZERO, stepA);
	//gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, one, eigenvaluematrix, inverseeigenvectors,one, stepA);
//	cout<<"stepA: \n"<<GSL_REAL(gsl_matrix_complex_get(stepA,0,0))<<" "<<GSL_REAL(gsl_matrix_complex_get(stepA,0,1))<<endl<<GSL_REAL(gsl_matrix_complex_get(stepA,1,0))<<"\t"<<GSL_REAL(gsl_matrix_complex_get(stepA,1,1))<<endl;
	gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, eigenvectors, stepA,GSL_COMPLEX_ZERO, transitionmatrixcomplex);
	gsl_eigen_nonsymmv_free(w);
	gsl_vector_complex_free(eigenvalues);
	gsl_matrix_complex_free(eigenvectors);
	gsl_matrix_complex_free(eigenvaluematrix);
	gsl_matrix_complex_free(inverseeigenvectorsstart);
	gsl_matrix_complex_free(inverseeigenvectors);
	gsl_matrix_complex_free(stepA);
	gsl_matrix_free(ratematrixTMP);
	gsl_permutation_free(p);
//	gsl_matrix *transitionmatrix=gsl_matrix_calloc(dimension,dimension);
	for (int i=0;i<dimension;i++) {
		for (int j=0;j<dimension;j++) {
			gsl_matrix_set(transitionmatrix,i,j,GSL_REAL(gsl_matrix_complex_get(transitionmatrixcomplex,i,j)));
			/*if (1==isnan(gsl_matrix_get(transitionmatrix,i,j))) {
				message="\nWarning: Got NaN in ComputeTransitionProb\n";
				PrintMessage();
			} */
		}
	}
	gsl_matrix_complex_free(transitionmatrixcomplex);
//	cout<<"transitionmatrixcomplex: \n"<<GSL_REAL(gsl_matrix_complex_get(transitionmatrixcomplex,0,0))<<" "<<GSL_REAL(gsl_matrix_complex_get(transitionmatrixcomplex,0,1))<<endl<<GSL_REAL(gsl_matrix_complex_get(transitionmatrixcomplex,1,0))<<"\t"<<GSL_REAL(gsl_matrix_complex_get(transitionmatrixcomplex,1,1))<<endl;
//	cout<<"transitionmatrix: \n"<<gsl_matrix_get(transitionmatrix,0,0)<<" "<<gsl_matrix_get(transitionmatrix,0,1)<<endl<<gsl_matrix_get(transitionmatrix,1,0)<<"\t"<<gsl_matrix_get(transitionmatrix,1,1)<<endl;
	return transitionmatrix;
	
}


/** @method HandleTimeSlice
*
* Allows user to specify how to slice trees so that different time periods have different rates.
*
*/
void BROWNIE::HandleTimeSlice( NexusToken& token )
{
	staterestrictionvector.clear();
	timeslicetimes.clear();
	timeslicemodels.clear();
	for (int state=0;state<=(maxModelCategoryStates-1);state++) {
        staterestrictionvector.push_back(state);
        timeslicetimes.push_back(GSL_POSINF);
        timeslicemodels.push_back(0);
    }
	
    bool noerror=true;
    bool adequateinput=false;
	bool enteredsplit=false;
	bool enteredmodel=false;
    for(;;)
    {
        token.GetNextToken();
        if( token.Equals(";") ) {
            if (adequateinput==false) {
                message="Insufficient input: type \"timeslice ?\" for help";
                PrintMessage();
            }
            break;
        }
        else if( token.Equals("?") ) {
            adequateinput=true;
            message="Usage: Timeslice splits=(<times>) models=(<numbers>) [options...]\n\n";
            message+="Assigns models to intervals.\n\n";
            message+="Available options:\n\n";
            message+="Keyword ---- Option type ------------------------ Current setting --";
            message+="\nSplits       (<numbers>)                          *None";
            message+="\nModels       (<numbers>)                          *None";
            // message+="\nTaxset       <taxset name>                        *All";
            // message+="\nFromRoot     No|Yes                               *No";
            message+="\n                                                 *Option is nonpersistent\n\n";
            message+="This will allow assignment of different models on the specified intervals across all trees.\nFor example, to specify one rate parameter for the interval from the present to 50 MYA,\nanother rate parameter from 50 MYA to 65 MYA, and the first rate parameter again from 65 MYA\nto the root of the phylogeny:\n\n  Time Model\n   0     1\n   5     1\n   .     .\n  45     1\n__50_____1_\n| 50     2 |\n| 55     2 |\n| 60     2 |\n|_65_____2_|\n  65     1\n  70     1\n   .     .\n\none would type\n\nTimeslice splits=(50 65) models=(1 2 1);\n\nSplits specifies where one model changes to the other;\nthere should be one fewer split than there are intervals.\nEdges spanning a split are properly divided\n(so a given path may have more than one model).\n\nBy default, times are measured from the most recent taxon. \nTo specify that they should be measured from the root instead, specify fromroot=yes [NOT IMPLEMENTED YET];\nThis becomes most important when the root node may be of different depths in different trees.\nNote that you are currently limited to only 9 splits (if this becomes a problem, email omeara.brian@gmail.com).\nDoing this command will overwrite any previous splits or discrete character mappings.";
            PrintMessage();
        }
        else if( token.Abbreviation("Splits") ) {
			enteredsplit=true;
            token.GetNextToken();
			token.GetNextToken(); //eat the equals sign
            if (!token.Equals("(")) {
                errormsg="Expecting next token to be a \'(\' but instead got ";
                errormsg+=token.GetToken();
                throw XNexus( errormsg);
            }
            token.GetNextToken(); //Should be a number
            int splitposcount=0;
			while (!token.Equals(")")) {
				//cout<<"token is "<<token.GetToken();
                nxsstring numbernexus;
                numbernexus=token.GetToken();
				token.GetNextToken();
                timeslicetimes[splitposcount]=atof( numbernexus.c_str() );
                splitposcount++;
            }
        }
        else if( token.Abbreviation("Models") ) {
            token.GetNextToken();
			token.GetNextToken(); //eat the equals sign
			enteredmodel=true;
            if (!token.Equals("(")) {
                errormsg="Expecting next token to be a \'(\' but instead got ";
                errormsg+=token.GetToken();
                throw XNexus( errormsg);
            }
            token.GetNextToken(); //Should be a number
			int modelposcount=0;
            while (!token.Equals(")")) {
				//cout<<"token is "<<token.GetToken();
                nxsstring numbernexus;
                numbernexus=token.GetToken();
				token.GetNextToken();
				//cout<<" modelposcount="<<modelposcount<<" number="<<numbernexus.c_str()<<" numbernexus="<<numbernexus<<endl;
                timeslicemodels[modelposcount]=atoi( numbernexus.c_str() )-1;
                modelposcount++;
                //token.GetNextToken();
            }
        }
        else if (token.Abbreviation("Taxset") ) {
        }
        else if (token.Abbreviation("FromRoot")  ) {
        }
        else {
            errormsg = "Unexpected keyword (";
            errormsg += token.GetToken();
            errormsg += ") encountered reading RateTest command. Type \"RateTest ?\" for help.";
            throw XNexus( errormsg, token.GetFilePosition(), token.GetFileLine(), token.GetFileColumn() );
        }
        //timeslicetimes[9]=GSL_POSINF;
        //timeslicemodels;
		
    }
    //if (__________EVERYTHING IS OKAY__________) {
	if (enteredsplit && enteredmodel) { //Everything is okay?
		adequateinput=true;
		int originalchosentree=chosentree;
		for (chosentree = 1; chosentree <= trees->GetNumTrees(); chosentree++) {
			Tree *Tptr=&(intrees.Trees[chosentree-1]);
			(*Tptr).SetPathLengths();
			double MaxLength=(*Tptr).GetMaxPathLength();
			//cout<<"MaxPathLength = "<<MaxLength<<endl;
			NodeIterator <Node> n ((*Tptr).GetRoot());
			Node *q = n.begin();
			while (q)
			{
				if (q!=(*Tptr).GetRoot()) {
					double BeginTime=MaxLength-(q->GetPathLength());
					double EndTime=BeginTime+(q->GetEdgeLength());
				//cout<<endl<<endl<<"---------------------------------"<<endl<<"Begin time="<<BeginTime<<" EndTime="<<EndTime<<endl<<"PathLength="<<q->GetPathLength()<<endl<<"EdgeLength="<<q->GetEdgeLength()<<endl;
					vector<double> modelvector(maxModelCategoryStates,0.0); //maxModelCategoryStates is in TreeLib.h
					vector<int> stateordervector;
					vector<double> statetimesvector; 
					double RegimeStartTime=0;
					double RegimeEndTime=MaxLength;
					for (int i=0; i<timeslicemodels.size(); i++) {
						if (i>0) {
							RegimeStartTime=timeslicetimes[i-1];
							if (gsl_isinf(timeslicetimes[i-1])!=0) {
								RegimeStartTime=MaxLength;
							}
							
						}
						if (gsl_isinf(timeslicetimes[i])==0) {
							RegimeEndTime=timeslicetimes[i];
						}
						else {
							RegimeEndTime=MaxLength;
						}
						//cout<<endl<<"Regime "<<i<<" state "<<timeslicemodels[i]<<" "<<RegimeStartTime<<"-"<<RegimeEndTime;
						if (BeginTime<RegimeEndTime && BeginTime<EndTime) {
					//	cout<<"  IN REGIME";
							double timeinstate=GSL_MIN(RegimeEndTime,EndTime)-BeginTime;
							modelvector[(timeslicemodels[i])]+=timeinstate;
							stateordervector.push_back(timeslicemodels[i]);
							statetimesvector.push_back(timeinstate);
							BeginTime=GSL_MIN(RegimeEndTime,EndTime);
						}
					}
					q->SetModelCategory(modelvector); 
					q->SetStateOrder(stateordervector); 
					q->SetStateTimes(statetimesvector);
					/*	cout<<endl<<"modelvector"<<endl;
					for (int i=0; i<modelvector.size(); i++) {
						cout<<modelvector[i]<<"\t";
					}
					cout<<endl<<"stateordervector"<<endl;
					for (int i=0; i<stateordervector.size(); i++) {
						cout<<stateordervector[i]<<"\t";
					}
					
					cout<<endl<<"statetimesvector"<<endl;
					for (int i=0; i<statetimesvector.size(); i++) {
						cout<<statetimesvector[i]<<"\t";
					} */
					
				}
				q = n.next();
			}
		}
		chosentree=originalchosentree;
	}
    //Do the set model thing
    //}
	
}



int main(int argc, char *argv[])
{
    //Code on random number generation ripped from GSL documentation



    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_mt19937;
    r = gsl_rng_alloc (T);
    BROWNIE brownie;
	//gsl_set_error_handler_off();
    bool inputfilegiven=false;
    if (argc>1) {
        for (int i = 1; i < argc; i++) {
            inputfilegiven=true;
            nxsstring fn=argv[i];
            brownie.RunCmdLine(inputfilegiven, fn );
        }
    }
    else {
        nxsstring fn="a";
        brownie.RunCmdLine(inputfilegiven, fn);
    }
    return 0;
}


