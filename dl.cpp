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



// globals and structures

// function declarations

using namespace std;


/*
* File for testing how to get tree information from the BROWNIE object 
* (for eventual return to R)
*/
int main()
{
	
	// tyring using some important stuff	
	printf("Initializing BROWNIE object...\n");
	BROWNIE brownie;
	brownie.Init();
	printf("done\n");
	
	std::string pre = "execute ";
	std::string fin;
	cout<<"please enter a file to execute: ";
	cin>>fin;
	
	// load in text file
	printf("Executing text file...");
	//cout << "preload status: "<< brownie.intrees.GetNumTrees()<<endl;
	strcpy(brownie.next_command,pre.append(fin).c_str());
	brownie.PreprocessNextCommand();
	printf("\n .. conditioned command is: %s\n",brownie.next_command);
   	brownie.HandleNextCommand();
	printf(" ...DONE\n");
	
	// How to get information from a BROWNIE object
	LabelList nn;
	int nchar = (*brownie.characters).GetNCharTotal();
	int nchard = 0;
	if(brownie.discretecharloaded)
		nchard = (*brownie.discretecharacters).GetNCharTotal();
	int ntaxa = (*brownie.taxa).GetNumTaxonLabels();
	int ntrees = (*brownie.trees).GetNumTrees();
	int nass = (*brownie.assumptions).GetNumTaxSets();
	cout << "Number of trees: "<< brownie.intrees.GetNumTrees()<< " ("<< ntrees <<")"<<endl;
	cout << "Number of taxa: "<< ntaxa <<endl;
	cout << "Number of assumptions: "<< nass <<endl;
	cout << "Number of characters total: "<<nchar<<endl;
	cout << "Number of taxa total: "<< (*brownie.characters).GetNTaxTotal()<<endl;
	cout << "Discrete chars loaded? "<< brownie.discretecharloaded<<endl;
	cout << "--Number of disc. chars: " << nchard <<endl;
	
	//(*brownie.assumptions).GetTaxSetNames(nn);
	//(*brownie.assumptions).Report(cout);
	//(*brownie.trees).Report(cout);
	
	
	if(ntrees > 0)
	{
		// Print translated trees
		//
		for (int j=0; j<ntrees; j++)
		{
			//cout<<(*brownie.trees).GetTranslatedTreeDescription(j)<<endl;
			//cout << brownie.intrees.GetIthTree(j).GetWeight() <<endl;
		}
		cout<<"Writing trees back....";
		ofstream output;
		cout<<"is_open=="<<output.is_open()<<"....";
		output.open("asdfasdf.txt");
		if(brownie.intrees.WriteTrees(output))
			cout<<"done";
		else
			cout<<"FAILED";
		
		output.close();
		cout<<endl;
	} else {
		cout << "no trees" << endl;
	}
	

	if(nchar > 0)
	{	
		// Print characters vectors, including their names
		cout<<"Character labels: "<<endl;
		for(int i=0;i<nchar;i++)
		{
			cout<<(*brownie.characters).GetCharLabel(i)<< " ";
		}
		cout<<endl;	
		
		// Characters Block
		for(int j =0; j < ntaxa; j++)
		{
			
			for(int i=0;i<nchar;i++)
			{
				if(!brownie.discretecharloaded)
				{
					cout << (*brownie.characters).GetValue(j,i,false) << " ";
				} else {
					cout << "State: "<<(*brownie.characters).GetState(j,i);
					cout << " Value: "<< (*brownie.characters).GetValue(j,i,false);
					cout << " -- ";
				}
			}
			cout << endl;
		}
	} else {
		cout << "no characters" << endl;
	}
	
	// Assumptions Block
	// --- set nexusdefs.h for LabelList and IntSet definitions
	//
	if(nass > 0)
	{
		LabelList assnames(nass);  // std::vector<nxsstring>
		(*brownie.assumptions).GetTaxSetNames(assnames);
		IntSet* isets = new IntSet[nass];  // std::set<int,<less>>
		
		for(int j=0; j < nass; j++)
		{
			cout << assnames[j] << " - " ;
			if(assnames[j].substr(0,3).compare("NOT")!=0 && assnames[j].compare("ALL")!=0)
			{	
				isets[j] = (*brownie.assumptions).GetTaxSet(assnames[j]);
				cout << isets[j].size()<<endl;
				if(!isets[j].empty())
				{
					cout << "Taxa:";
					for(IntSet::iterator kk=isets[j].begin(); kk != isets[j].end(); kk++)
						cout << " " << (*brownie.taxa).GetTaxonLabel(*kk); // *kk is an int, gettaxonlabel returns nxsstring
					cout << endl;
				}
			} else {
				cout << "do not return"<<endl;
			}
			
		}		
		delete [] isets;
	} else {
		cout <<"no assumptions"<<endl;
	}
	
	return 0;
}



