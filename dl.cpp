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
	
	cout << "\n\n------------\nSUMMARY INFORMATION:\n-------------\n";;	
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
	
	cout << "\n\n------------\nTREE STUFF\n-------------\n";;
	if(ntrees > 0)
	{
		// Print translated trees
		//
		cout<<"WEIGHTS:";
		for (int j=0; j<ntrees; j++)
		{
			//cout<<(*brownie.trees).GetTranslatedTreeDescription(j)<<endl;
			//cout<<(*brownie.trees).GetTreeDescription(j)<<endl;
			//cout << brownie.intrees.GetIthTree(j).Weight() <<endl;
		}
		cout<<endl;
		
		cout<<"TREE NEXUS STRING:";
		cout<<" (the chosen tree is at index: "<<brownie.chosentree<<")"<<endl;
		//stringstream ss(stringstream::in | stringstream::out); 
		ostringstream ss;
		//if(brownie.intrees.WriteTrees(ss))
			//cout<<ss.str();
		//cout<<endl;
		if(ss.str().length() == 0)
			cout<<"Nothing in there yet!"<<endl;
		
			// Try to get simmap formatted tree..
		brownie.intrees.GetIthTree(brownie.chosentree-1).WriteNoQuote(ss);
		cout<<ss.str()<<endl;
		ss.str(std::string());
		cout<<"Num chars now: "<<ss.str().length()<<endl;
		cout<<endl;
		
		cout<<"TREE LABELS:";
		//brownie.intrees.ShowLabelList(cout);
		cout<<endl;
		
		cout<<"Writing trees to file....";
		ofstream output;
		cout<<"is_open=="<<output.is_open()<<"....";
		output.open("asdfasdf.txt");
		if(brownie.intrees.WriteTrees(output))
			cout<<"done";
		else
			cout<<"FAILED";
		output.close();
		cout<<endl;
		
		// Check rettree:
		if(brownie.rettree.str().length() != 0)
			cout << brownie.rettree.str() << endl;
		
	} else {
		cout << "no trees" << endl;
	}
	
	cout << "\n\n------------\nCHARACTER STUFF\n-------------\n";;
	if(nchar > 0)
	{	
		// Print characters vectors, including their names
		cout<<"Character labels: "<<endl;
		for(int i=0;i<nchar;i++)
		{
			cout << (*brownie.characters).GetCharLabel(i) << " ";
		}
		cout<<endl;	
		
		// Characters Block
		for(int j =0; j < ntaxa; j++)
		{
			
			for(int i=0;i<nchar;i++)
			{
				
				// NOTE: GetState returns char and GetValues return float
				if(!brownie.discretecharloaded)
				{
					cout << (*brownie.characters).GetValue(j,i,false) << " "; // GetValue returns a float
				} else {
					cout << "State: "<<(*brownie.characters).GetState(j,i); // this one is the one to use
					cout << " Value: "<< (*brownie.characters).GetValue(j,i,false);
					cout << " -- ";
				}
			}
			cout << endl;
		}
	} else {
		cout << "no characters" << endl;
	}
	
	vector<string> taxanames(nass);
	vector< vector<string> > taxasets(nass);
	int currindex = 0;
	
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
			if(assnames[j].substr(0,3).compare("NOT")!=0 && assnames[j].compare("ALL")!=0)
			{	
				taxanames[currindex] = (std::string)assnames[j];
				cout << taxanames[currindex] << " - " ;
				// IntSets:
				isets[j] = (*brownie.assumptions).GetTaxSet(assnames[j]);
				
				if(!isets[j].empty())
				{
					cout << "Taxa:";
					for(IntSet::iterator kk=isets[j].begin(); kk != isets[j].end(); kk++)
					{
						cout << " " << (*brownie.taxa).GetTaxonLabel(*kk); // *kk is an int, gettaxonlabel returns nxsstring
						taxasets[currindex].push_back((*brownie.taxa).GetTaxonLabel(*kk));
					}
					
					//cout << isets[j].size() << endl;
					cout << "Total: " << taxasets[currindex].size() << endl;
					cout << endl;
					currindex++;
				}
			} else {
				cout << "do not return"<<endl;
			}
			
		}		
		delete [] isets;
	} else {
		cout <<"no assumptions"<<endl;
	}
	
	taxanames.resize(currindex);
	taxasets.resize(currindex);
	
	cout<<currindex<<" taxa sets used:"<<endl;
	for(int iii = 0; iii < taxanames.size(); iii++)
		cout<<taxanames[iii]<<endl;
	
	return 0;
}



