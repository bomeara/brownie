/*
 * TreeLib
 * A library for manipulating phylogenetic trees.
 * Copyright (C) 2001 Roderic D. M. Page <r.page@bio.gla.ac.uk>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307, USA.
 */

// $Id: TreeLib.cpp,v 1.6 2006/06/18 22:01:09 bcomeara Exp $

/**
* This file has been modified by
 * Brian O'Meara. August 2005.
 * My modifications generally have
 * BCO preceding them or //added by BCO following them.
 */


#include "TreeLib.h"
#include "Parse.h"

#include <vector>
//BCO added iostream just to give more feedback on errors
#include <iostream>
#include <cstdlib>
#include <climits>
#include <cstring>
#include <memory>
#include <cassert>

// Convert a string to a NEXUS format string
std::string NEXUSString (const std::string s)
{
	std::string outputString ="";
	bool enclose = false;
	int i = 0;
	
	if (isalpha (s[0]))
	{
		// First character is a letter, check the rest
		i = 1;
		while ((i < s.length()) && !enclose)
		{
			if (!isalnum (s[i]) && (s[i] != ' ') && (s[i] != '_') && (s[i] != '.'))
				enclose = true;
			i++;
		}
	}
	else
		enclose = true;
	
	if (enclose)
	{
		outputString = "'";
	}
	
	i = 0;
	while (i < s.length())
	{
		if (s[i] == '\'')
			outputString += "''";
		else if ((s[i] == ' ') && !enclose)
			outputString += '_';
		else
			outputString += s[i];
		i++;
	}
	
	
	if (enclose)
	{
		outputString += "'";
	}
	
	return outputString;
}




//------------------------------------------------------------------------------
Node::Node ()
{
	Child = Anc = Sib = NULL;
	Label 	= "";
	Weight 	= 0;
	vector<double> ModelCategory(1,1.0); //added by BCO
	vector<int> StateOrder; //Added by BCO
	vector<double> StateTimes; //Added by BCO
	Length 	= 0.0;
	Leaf 	= false;
	Height 	= 0;
	Marked 	= false;
	Depth   = 0;
	Degree	= 0;
	PathLength = 0.0;
	LeafNumber = 0;
	LabelNumber = 0;
}

//------------------------------------------------------------------------------
void Node::Copy (Node *theCopy)
{
	theCopy->SetLeaf (IsLeaf ());
	theCopy->SetLabel (Label);
	theCopy->SetIndex (Index);
	theCopy->SetLeafNumber (LeafNumber);
	theCopy->SetLabelNumber (LabelNumber);
	theCopy->SetEdgeLength (Length);
	theCopy->SetModelCategory (ModelCategory); //Added by BCO
	theCopy->SetStateOrder (StateOrder); //Added by BCO
	theCopy->SetStateTimes (StateTimes); //Added by BCO
}

///------------------------------------------------------------------------------
// True if node is a left descendant of q
bool Node::IsALeftDescendantOf (Node *q)
{
	NodePtr r = this;
	while ((r) && (r != q->Child) && (r != q)) r = r->Anc;
	return (r == q->Child);
}

//------------------------------------------------------------------------------
// Return sibling node that is immediately to the LEFT of this node
Node *Node::LeftSiblingOf ()
{
	NodePtr q = Anc->Child;
	while (q->Sib != this)
		q = q->Sib;
	return q;
}

//------------------------------------------------------------------------------
// Return the rightmost sibling of the node
Node *Node::GetRightMostSibling ()
{
	Node *p = this;
	
	while (p->Sib) p = p->Sib;
	return (p);
}

//------------------------------------------------------------------------------
Tree::Tree ()
{
	Root 			= NULL;
	CurNode 		= NULL;
	Leaves 		= 0;
	Internals 		= 0;
	Error 		= 0;
	InternalLabels 	= false;
	EdgeLengths 	= false;
	Nodes 		= NULL;
	Name 		= "";
	Rooted 		= false;
	Weight		= 1.0;
}

//------------------------------------------------------------------------------
// Copy constructor
Tree::Tree (const Tree &t)
{
    if (t.GetRoot() == NULL)
    {
		Root = NULL;
        CurNode 		= NULL;
        Leaves 		= 0;
        Internals 		= 0;
        Error 		= 0;
        InternalLabels 	= false;
        EdgeLengths 	= false;
        Nodes 		= NULL;
        Name 		= "";
        Rooted 		= false;
        Weight		= 1.0;
    }
    else
    {
		CurNode 		= t.GetRoot();   		
		NodePtr placeHolder;  				
		t.copyTraverse (CurNode, placeHolder );
		Root 			= placeHolder;  
		Leaves    		= t.GetNumLeaves ();
		Internals 		= t.GetNumInternals ();
		Name	  		= t.GetName ();
		CurNode 		= NULL;
		Error 		= 0;
		InternalLabels 	= t.GetHasInternalLabels ();;
		EdgeLengths 	= t.GetHasEdgeLengths ();
		Nodes 		= NULL;
		Rooted 		= t.IsRooted();
		Weight		= t.GetWeight();
		
	}
}


//------------------------------------------------------------------------------
Tree::~Tree ()
{
	deletetraverse (Root);
	delete [] Nodes;
	Nodes = NULL; //added by BCO; recommended to prevent invalid memory reference (by web C++ tutorials)
}

//------------------------------------------------------------------------------
void Tree::deletetraverse (NodePtr p)
{
	if (p)
	{
		deletetraverse (p->GetChild());
		deletetraverse (p->GetSibling());
		delete p;
	}
}

//------------------------------------------------------------------------------
/*
	Return a copy of the subtree in rooted at RootedAt.
 The function result is the root of the subtree.

 Tree

            x---y    b---c               x'--y'
             \        \                   \
              \        \                   \
    RootedAt ->w--------a       ==>        CopyOfSubTree -> NULL
                \
                 \

 which corresponds to copying

 x'    y'             x     y b     c
  \   /                \   /   \   /
   \ /                  \ /     \ /
    w'       from        w       a
                          \     /
                           \   /
                            \ /
                             v

 */

//------------------------------------------------------------------------------
// This code needs Tree to be a friend of Node
void Tree::copyTraverse (NodePtr p1, NodePtr &p2) const
{
	if (p1)
	{
		p2 = NewNode ();
		p1->Copy (p2);
		
		// Note the call to p2->child, not p2->GetChild(). Calling the field
		// is essential because calling GetChild merely passes a temporary
		// copy of the child, hence we are not actually creating a child of p2.
		// We can access child directly by making Tree a friend of Node (see
		// TreeLib.h).
		copyTraverse (p1->GetChild(), p2->Child);
		if (p2->GetChild())
			p2->GetChild()->SetAnc (p2);
		
		// Ensure we don't copy RootedAt sibling. If the sibling is NULL then
		// we won't anyway, but this line ensures this for all cases.
		if (p1 != CurNode)
			copyTraverse (p1->GetSibling(), p2->Sib);  // note sib
		if (p2->GetChild ())
		{
			NodePtr q = p2->GetChild()->GetSibling();
			while (q)
			{
				q->SetAnc (p2);
				q = q->GetSibling();
			}
		}
	}
}

//------------------------------------------------------------------------------
NodePtr Tree::CopyOfSubtree (NodePtr RootedAt) 
{
	CurNode = RootedAt;   // Store this to avoid copying too much of the tree
	NodePtr placeHolder;  // This becomes the root of the subtree
	copyTraverse (CurNode, placeHolder);
	return placeHolder;
}


//------------------------------------------------------------------------------
string Tree::GetErrorMsg ()
{
	string s = "No error";
	
	switch (Error)
	{
		case errSYNTAX:
			s = "Syntax error";
			break;
		case errENDOFSTRING:
			s = "Syntax error";
			break;
		case errMISSINGLPAR:
			s = "Missing '('";
			break;
		case errUNBALANCED:
			s = "Unbalance parentheses";
			break;
		case errSTACKNOTEMPTY:
			s = "Stack not empty";
			break;
		case errSEMICOLON:
			s = "Expecting a semicolon";
			break;
		case errTOOMANYMODELCATS: //added by BCO
			s = "Too many states in Simmap code"; //added by BCO
			break; //added by BCO
		default:
			break;
	}
	return s;
}

//------------------------------------------------------------------------------
// Make CurNode a leaf, called by tree reading code
void Tree::MakeCurNodeALeaf (int i)
{
	Leaves++;
	CurNode->SetLeaf(true);
	CurNode->SetWeight(1);
	//	CurNode->SetLabelNumber(i);
	CurNode->SetLeafNumber(Leaves);
	//	Nodes[i - 1] = CurNode;
}

//m

//------------------------------------------------------------------------------
// Make a child of CurNode and make it CurNode
void Tree::MakeChild ()
{
	NodePtr	q = NewNode();
	CurNode->SetChild(q);
	q->SetAnc(CurNode);
	CurNode->IncrementDegree();
	CurNode = q;
	Internals++;
}

//------------------------------------------------------------------------------
// Create a new node and make it the root of the tree, and set CurNode=Root.
void Tree::MakeRoot ()
{
	CurNode   = NewNode();
	Root 	  = CurNode;
}

//------------------------------------------------------------------------------
// Make a sibling of CurNode and make CurNode the new node.
void Tree::MakeSibling ()
{
	NodePtr	q = NewNode ();
	NodePtr	ancestor = CurNode->GetAnc();
	CurNode->SetSibling(q);
	q->SetAnc(ancestor);
	ancestor->AddWeight(CurNode->GetWeight());
	ancestor->IncrementDegree();
	CurNode = q;
}



//------------------------------------------------------------------------------
int Tree::Parse (const char *TreeDescr)
{
	enum {stGETNAME, stGETINTERNODE, stNEXTMOVE,
		stFINISHCHILDREN, stQUIT, stACCEPTED} state;
	
	stack< NodePtr, vector<NodePtr> > stk;
	NodePtr		q;
	Parser		p ((char *)TreeDescr);
	tokentype	token;
	float 		f;
	
	// Initialise tree variables
	Root 		= NULL;
	Leaves 		= 0;
	Internals 	= 0;
	Error 		= 0;
	
	// Create first node
	CurNode		= NewNode();
	Root		= CurNode;
	
	// Initialise FSA that reads tree
	state = stGETNAME;
	//std::cout<<endl<<"400 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
	token = p.NextToken ();
	while ((state != stQUIT) && (state != stACCEPTED))
	{
		switch (state)
		{
			
			case stGETNAME:
				switch (token)
				{
					case SPACE:
					case TAB:
					case NEWLINE:
						//std::cout<<endl<<"413 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
						token = p.NextToken ();
						break;
					case STRING:
						// to do: handle translation
						Leaves++;
						CurNode->SetLeaf (true);
						CurNode->SetLeafNumber (Leaves);
						CurNode->SetWeight (1);
						//std::cout<<endl<<"425 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
						CurNode->SetLabel (p.GetToken());
						CurNode->SetDegree (0);
						//CurNode->SetModelCategory(vector<double>(1,1.0)); //Added by BCO
						//std::cout<<endl<<"425 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
						token = p.NextToken ();
						state = stGETINTERNODE;
						break;
					case NUMBER:
						// to do: handle translation
						Leaves++;
						CurNode->SetLeaf (true);
						CurNode->SetLeafNumber (Leaves);
						CurNode->SetWeight (1);
						CurNode->SetLabel (p.GetToken());
						CurNode->SetDegree (0);
						//  CurNode->SetModelCategory(vector<double>(1,1.0)); //Added by BCO
						
						//std::cout<<endl<<"439 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
						token = p.NextToken ();
						state = stGETINTERNODE;
						break;
					case LPAR:
						state = stNEXTMOVE;
						break;
					case ENDOFSTRING:
						Error = errENDOFSTRING;
						state = stQUIT;
						break;
					default:
						Error = errSYNTAX;
						std::cerr<<"Had error in Parse() in TreeLib.cpp line 447"; //added by BCO
						state = stQUIT;
						break;
				}
				break;
				
			case stGETINTERNODE:
				switch (token)
				{
					case SPACE:
					case TAB:
					case NEWLINE:
						//std::cout<<endl<<"464 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
						token = p.NextToken ();
						break;
					case COLON:
					case COMMA:
					case RPAR:
						state = stNEXTMOVE;
						break;
					case ENDOFSTRING:
						Error = errENDOFSTRING;
						state = stQUIT;
						break;
					default:
						Error = errSYNTAX;
						std::cerr<<"Had error in Parse() in TreeLib.cpp line 472"; //added by BCO
						state = stQUIT;
						break;
				}
				break;
				
			case stNEXTMOVE:
				switch (token)
				{
					case COLON: //BCO, have it here deal separately with the next token being a number (use default behavior of treating it as a brlen) and having it a left curly bracket, indicating that simmap output is following.
								//std::cout<<endl<<"488 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
						token = p.NextToken ();
						char matchchar;
						matchchar=p.GetTokenAsCstr()[0];
						if (matchchar=='{') { //added by BCO
							
							//  wrote everything in this if statement
							//Note that the simmap parser only deals with "morphological" data -- states should be integers, in other words. //BCO
							//State order goes from rootmost end of the branch to the tipmost end //BCO
							//StateTimes[i] is amount of time in StateOrder[i] //BCO
							vector<double> modelvector(maxModelCategoryStates,0.0); //maxModelCategoryStates is in TreeLib.h
							vector<int> stateordervector; //Added by BCO
							vector<double> statetimesvector; //Added by BCO
							double totaledgelength=0;
							char internalmatchchar;
							internalmatchchar=p.GetTokenAsCstr()[0];
							while (internalmatchchar!='}') {
								//std::cout<<endl<<"502 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
								token=p.NextToken();
								char colonmatchchar;
								colonmatchchar=p.GetTokenAsCstr()[0];
								if (colonmatchchar==':') {
									//std::cout<<endl<<"506 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
									token=p.NextToken(); //Get rid of colon
									
								}
								int state=atoi(p.GetTokenAsCstr());
								if (state>=maxModelCategoryStates) {
									Error=errTOOMANYMODELCATS;
									break;
								}
								//std::cout<<endl<<"515 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
								token=p.NextToken(); //get rid of comma
													 //std::cout<<endl<<"517 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
								token=p.NextToken(); //This should have the branch length spent in the above state
								double brleninstate = atof (p.GetTokenAsCstr());
								modelvector[state]+=brleninstate;
								stateordervector.push_back(state);
								statetimesvector.push_back(brleninstate);
								totaledgelength+=brleninstate;
								//std::cout<<endl<<"523 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO								
								token=p.NextToken();  //take in next token
								internalmatchchar=p.GetTokenAsCstr()[0];
							}
							CurNode->SetEdgeLength (totaledgelength);
							CurNode->SetModelCategory(modelvector); //Added by BCO
							CurNode->SetStateOrder(stateordervector); //Added by BCO
							CurNode->SetStateTimes(statetimesvector); //Added by BCO
							EdgeLengths = true;
							//std::cout<<endl<<"527 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
							token = p.NextToken ();
							break;
						} //BCO
							else { //BCO. Will work when not dealing with a simmap tree.
								f = atof (p.GetTokenAsCstr());
								CurNode->SetEdgeLength (f);
								CurNode->SetModelCategory(vector<double>(1,f)); //Added by BCO
								CurNode->SetStateOrder(vector<int>(1,0)); //Added by BCO
								CurNode->SetStateTimes(vector<double>(1,f)); //Added by BCO
								EdgeLengths = true;
								//std::cout<<endl<<"536 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
								token = p.NextToken ();
							} //BCO
							break;
                    case SPACE:
					case TAB:
					case NEWLINE:
						//std::cout<<endl<<"543 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
						token = p.NextToken ();
						break;
						// The next node encountered will be a sibling
						// of Curnode and a descendant of the node on
						// the top of the node stack.
					case COMMA:
						q = NewNode();
						CurNode->SetSibling (q);
						if (stk.empty())
						{
							Error = errMISSINGLPAR;
							state = stQUIT;
						}
							else
							{
								q->SetAnc (stk.top());
								stk.top()->AddWeight (CurNode->GetWeight());
								stk.top()->IncrementDegree ();
								CurNode = q;
								state = stGETNAME;
								//std::cout<<endl<<"564 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
								token = p.NextToken ();
							}
							break;
						// The next node will be a child of CurNode, hence
						// we create the node and push CurNode onto the
						// node stack.
					case LPAR:
						Internals++;
						stk.push (CurNode);
						q = NewNode();
						CurNode->SetChild (q);
						q->SetAnc (CurNode);
						CurNode->IncrementDegree ();
						CurNode = q;
						//std::cout<<endl<<"579 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
						token = p.NextToken ();
						state = stGETNAME;
						break;
						// We've finished ready the descendants of the node
						// at the top of the node stack so pop it off.
					case RPAR:
						if (stk.empty ())
						{
							Error = errUNBALANCED;
							state = stQUIT;
						}
						else
						{
							q = stk.top();
							q->AddWeight (CurNode->GetWeight());
							CurNode = q;
							stk.pop();
							state = stFINISHCHILDREN;
							//std::cout<<endl<<"598 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
							token = p.NextToken ();
						}
						break;
						// We should have finished the tree
					case SEMICOLON:
						if (stk.empty())
						{
							state = stACCEPTED;
						}
						else
						{
							Error = errSTACKNOTEMPTY;
							state = stQUIT;
						}
						break;
					case ENDOFSTRING:
						Error = errENDOFSTRING;
						state = stQUIT;
						break;
					default:
						Error = errSYNTAX;
						std::cerr<<"Had error in Parse() in TreeLib.cpp line 597"; //added by BCO
						state = stQUIT;
						break;
				}
				break;
				
			case stFINISHCHILDREN:
				switch (token)
				{
					case STRING:
					case NUMBER:
						// internal label
						InternalLabels = true;
						CurNode->SetLabel (p.GetToken());
						//std::cout<<endl<<"634 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
						token = p.NextToken ();
						break;
					case COLON: //BCO, have it here deal separately with the next token being a number (use default behavior of treating it as a brlen) and having it a left curly bracket, indicating that simmap output is following.
								//std::cout<<endl<<"640 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
						token = p.NextToken ();
						char matchchar;
						matchchar=p.GetTokenAsCstr()[0];
						if (matchchar=='{') { //added by BCO
							
							//BCO wrote everything in this if statement
							//Note that the simmap parser only deals with discrete data -- states should be integers, in other words. //BCO
							
							vector<double> modelvector(maxModelCategoryStates,0.0); //maxModelCategoryStates is in TreeLib.h
							vector<int> stateordervector; //added by BCO
							vector<double> statetimesvector; //added by BCO
							double totaledgelength=0;
							char internalmatchchar;
							internalmatchchar=p.GetTokenAsCstr()[0];
							while (internalmatchchar!='}') {
								//std::cout<<endl<<"654 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
								token=p.NextToken();
								char colonmatchchar;
								colonmatchchar=p.GetTokenAsCstr()[0];
								if (colonmatchchar==':') {
									//std::cout<<endl<<"659 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
									token=p.NextToken(); //Get rid of colon
									
								}
								int state=atoi(p.GetTokenAsCstr());
								if (state>=maxModelCategoryStates) {
									Error=errTOOMANYMODELCATS;
									break;
								}
								//std::cout<<endl<<"668 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
								token=p.NextToken(); //get rid of comma
													 //std::cout<<endl<<"670 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
								token=p.NextToken(); //This should have the branch length spent in the above state
								double brleninstate = atof (p.GetTokenAsCstr());
								modelvector[state]+=brleninstate;
								stateordervector.push_back(state);
								statetimesvector.push_back(brleninstate);
								totaledgelength+=brleninstate;
								//std::cout<<endl<<"675 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO								
								token=p.NextToken();  //take in next token
								internalmatchchar=p.GetTokenAsCstr()[0];
							}
							CurNode->SetEdgeLength (totaledgelength);
							CurNode->SetModelCategory(modelvector); //Added by BCO
							CurNode->SetStateOrder(stateordervector); //Added by BCO
								CurNode->SetStateTimes(statetimesvector); //Added by BCO
							EdgeLengths = true;
							//std::cout<<endl<<"682 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
							token = p.NextToken ();
							break;
						} //BCO
							else { //BCO. Will work when not dealing with a simmap tree.
								f = atof (p.GetTokenAsCstr());
								CurNode->SetEdgeLength (f);
								CurNode->SetModelCategory(vector<double>(1,f)); //Added by BCO
								CurNode->SetStateOrder(vector<int>(1,0)); //Added by BCO
								CurNode->SetStateTimes(vector<double>(1,f)); //Added by BCO
								EdgeLengths = true;
								//std::cout<<endl<<"691 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
								token = p.NextToken ();
							} //BCO
							break;
					case SPACE:
					case TAB:
					case NEWLINE:
						//std::cout<<endl<<"653 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
						token = p.NextToken ();
						break;
						
						// We've completed traversing the descendants of the
						// node at the top of the stack, so pop it off.
					case RPAR:
						if (stk.empty())
						{
							Error = errUNBALANCED;
							state = stQUIT;
						}
						else
						{
							q = stk.top();
							q->AddWeight (CurNode->GetWeight());
							CurNode = q;
							stk.pop();
							//std::cout<<endl<<"671 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
							token = p.NextToken ();
						}
						break;
						
						// The node at the top of the stack still has some
						// descendants.
					case COMMA:
						q = NewNode();
						CurNode->SetSibling (q);
						if (stk.empty())
						{
							Error = errMISSINGLPAR;
							state = stQUIT;
						}
							else
							{
								q->SetAnc (stk.top());
								stk.top()->AddWeight (CurNode->GetWeight());
								stk.top()->IncrementDegree ();
								CurNode = q;
								state = stGETNAME;
								//std::cout<<endl<<"693 Current token is "<<p.GetTokenAsCstr()<<endl; //added by BCO
								token = p.NextToken ();
							}
							break;
					case SEMICOLON:
						state = stNEXTMOVE;
						break;
					default:
						if (stk.empty())
						{
							Error = errSEMICOLON;
						}
						else
						{
							Error = errSYNTAX;
							std::cerr<<"Had error in Parse() in TreeLib.cpp line 679"; //added by BCO
						}
						state = stQUIT;
						break;
				}
				break;
		}
	}
	// Handle errors
	if (state == stQUIT)
	{
		// Clean up to go here
		return (p.GetPos());
	}
	else
	{
		Root->SetWeight(Leaves);
		getPathLengths(Root); //added by BCO
		return (0);
	}
}



//------------------------------------------------------------------------------
int Tree::Read (istream &f)
{
	char ch = '\0';
	string str = "";
	while (f.good() && !f.eof() && (ch != ';'))
	{
		f.get (ch);
		str += ch;
	}
	return Parse (str.c_str());
}


//------------------------------------------------------------------------------
void Tree::Write (ostream &f)
{
	treeStream = &f;
	traverse (Root);
	f << ";";
}

//------------------------------------------------------------------------------
void Tree::traverse (NodePtr p)
{
	
	if (p)
	{
		if (p->IsLeaf())
		{
			*treeStream << NEXUSString (p->GetLabel());
			
			if (EdgeLengths)
			{
				*treeStream << ':' << p->GetEdgeLength ();
			}
		}
		else
		{
			*treeStream << "(";
		}
		
		traverse (p->GetChild());
		if (p->GetSibling())
		{
			*treeStream << ",";
		}
		else
		{
			if (p != Root)
			{
				*treeStream << ")";
				// 29/3/96
				if ((p->GetAnc()->GetLabel() != "") && InternalLabels)
				{
					*treeStream << '\'' << NEXUSString (p->GetAnc()->GetLabel ()) << '\'';
				}
				if (EdgeLengths && (p->GetAnc () != Root))
				{
					*treeStream << ':' << p->GetAnc()->GetEdgeLength ();
				}
			}
		}
		traverse (p->GetSibling());
	}
	
}

//Added by BCO
//------------------------------------------------------------------------------
void Tree::WriteNoQuote (ostream &f)
{
	treeStream = &f;
	traversenoquote (Root);
	f << ";";
}

//Added by BCO
//------------------------------------------------------------------------------
void Tree::traversenoquote (NodePtr p)
{
	
	if (p)
	{
		if (p->IsLeaf())
		{
			for (int i=0; i<(NEXUSString (p->GetLabel())).length(); i++) {
				if ((NEXUSString (p->GetLabel()))[i]!= '\'') {
					*treeStream <<(NEXUSString (p->GetLabel()))[i];
				}
			}
			//*treeStream<<NEXUSString (p->GetLabel());
			if (EdgeLengths)
			{
				*treeStream << ':' << p->GetEdgeLength ();
			}
		}
		else
		{
			*treeStream << "(";
		}
		
		traversenoquote (p->GetChild());
		if (p->GetSibling())
		{
			*treeStream << ",";
		}
		else
		{
			if (p != Root)
			{
				*treeStream << ")";
				// 29/3/96
				if ((p->GetAnc()->GetLabel() != "") && InternalLabels)
				{
					*treeStream <<  NEXUSString (p->GetAnc()->GetLabel ()) ; //here's the change from traverse
				}
				if (EdgeLengths && (p->GetAnc () != Root))
				{
					*treeStream << ':' << p->GetAnc()->GetEdgeLength ();
				}
			}
		}
		traversenoquote (p->GetSibling());
	}
	
}



//------------------------------------------------------------------------------
void Tree::drawInteriorEdge (NodePtr p)
{
	NodePtr r = p->GetAnc ();
	int stop = r->GetHeight();
	if (p->IsTheChild ())
	{
		// Visiting ancestor for the first time, so draw the
		// end symbol
		if (r == Root)
		{
			//			if (IsRooted ())
			Line[stop] = TEE;   // «
								//			else
								//				Line[stop] = VBAR;  // Ò
		}
		else
		{
			Line[stop] = TEE;      // «
		}
		
		
		// Draw branch itself
		if (r != Root)
		{
			// Line
			int start = r->GetAnc()->GetHeight();
			for (int i = start + 1; i < stop; i++)
			{
				Line[i] = HBAR; // 
			}
			// Start symbol
			if (start == stop)
				Line[start] = VBAR;     // Ò
			else if (r->IsTheChild ())
				Line[start] = LEFT;     // ò
			else if (r->GetSibling ())
				Line[start] = SIB;      // Ì
			else Line[start] = RIGHT;  // Ë
			
			//
			fillInAncestors (r);
		}
	}
	else
	{
		// Just draw nodes below
		Line[stop] = VBAR;
		fillInAncestors (p->GetSibling());
	}
	
	// Output the line
	Line[stop + 1] = '\0';
	
    drawLine (r, p->IsTheChild());
	
	/*
	 
	 *treeStream << Line.c_str();
	 
	 //	*treeStream		<< "h=" << r->GetHeight() << " w= "
	 //			<< r->GetWeight() << "-- ";
	 
	 
	 // Draw internal label, if present
	 string s = r->GetLabel();
	 if (s != "" && p->IsTheChild ())
	 *treeStream  << r->GetLabel();
	 *treeStream  << endl;
	 */
	// Clear the line for the next pass
	for (int i = 0; i < (Leaves + 2); i++) // to do: get a better limit
		Line[i] = ' ';
}

void Tree::drawLine (NodePtr p, bool isChild)
{
	*treeStream << Line.c_str();
	
	if (p->IsLeaf())
    {
		* treeStream<< " " << p->GetLabel () <<endl;
    }
    else
    {
        // Draw internal label, if present
        string s = p->GetLabel();
        if (s != "" && isChild)
            *treeStream  << p->GetLabel();
		
	 	*treeStream  << endl;
    }
}



//------------------------------------------------------------------------------
void Tree::drawPendantEdge (NodePtr p)
{
	NodePtr q = p->GetAnc();
	if (q == NULL)
	{
		// Handle the degenerate case of a tree with a single leaf
		Line = (char) HBAR;
		drawLine (p);	
	}
	else
	{
		int start = q->GetHeight();
		int stop = p->GetHeight();
		char	symbol;
		
		// Draw line between p and its ancestor
		int i;
		for (i = start + 1; i <= stop; i++)
			Line[i] = (char)HBAR;      // 
		
		// Find appropriate symbol for link to ancestor
		if (p == q->GetChild())
		{
			symbol = (char)LEFT;       // ò
		}
		else
		{
			// p is a sibling
			if (p->GetSibling())
				symbol = (char)SIB;     // Ì
			else symbol = (char)RIGHT; // Ë
		}
		Line[start] = symbol;
		
		// Fill in ancestors
		fillInAncestors (p);
		
		// Terminate line
		Line[stop + 1] = '\0';
		
		
		drawLine (p);
		/*		// Output line and taxon name
			*treeStream << Line.c_str() << " " << p->GetLabel ()
			//			<< "h=" << p->GetHeight() << " w= "
			//			<< p->GetWeight()
			<<endl;
		*/
		
		// Clear the line for the next pass
		for (i = 0; i < (Leaves + 2); i++) // to do: get a better limit for this
			Line[i] = ' ';
	}
}



//------------------------------------------------------------------------------
/*
	{ Put a VBAR symbol ("Ò") into buffer where needed to
		represent branches of tree under p. Two cases:
		
		.....  Ancestor
		>  Sibling
		<----  Child
		
		1. r is a sibling of the child of q, and q has a sibling, or
		
		xsib>r
		ò    .    .
		\   .   .
		\  .  .
		child\ . .|
						  \.. |
		q|sib>
							 \ |
		\| <- this branch passes under r
								\
		
		
		2. r is q's child and q is a sibling of another node.
		
						  r
		ò
							 \
							 |\
							 | \
		xsib|>q
							 | .
							 |.<- this branch passes under r
							 .
		.
		
	} */

//------------------------------------------------------------------------------
void Tree::fillInAncestors (NodePtr p)
{
	
	NodePtr q = p->GetAnc ();
	NodePtr r = p;
	while (q != Root)
	{
		if (
			(q->GetSibling () && !(r->IsTheChild ()))
			|| (!(q->IsTheChild ()) && r->IsTheChild())
			)
		{
			if (r ==p && q->GetHeight() == q->GetAnc()->GetHeight())
				Line[q->GetAnc()->GetHeight()] = SIB;
			else
            	Line[q->GetAnc()->GetHeight()] = VBAR;
		}
		r = q;
		q = q->GetAnc ();
	}
}





//------------------------------------------------------------------------------
void Tree::drawAsTextTraverse (NodePtr p)
{
	if (p)
	{
		drawAsTextTraverse (p->GetChild ());
		if (p->IsLeaf ())
			drawPendantEdge (p);
		if (p->GetSibling ())
			drawInteriorEdge (p);
		drawAsTextTraverse (p->GetSibling ());
	}
}

//------------------------------------------------------------------------------
void Tree::Draw (ostream &f)
{
	treeStream = &f;
	
	if (Root)
	{
		string s (Leaves + 2, ' ');
        Line = s;
		MaxHeight = 0;
		getNodeHeights (Root);
		drawAsTextTraverse (Root);
	}
	else f << "(No tree)" << endl;
}

//------------------------------------------------------------------------------
void Tree::getNodeHeights(NodePtr p)
{
	if (p)
	{
		p->SetHeight (Leaves - p->GetWeight ());
		if (p->GetHeight() > MaxHeight)
			MaxHeight = p->GetHeight();
		getNodeHeights (p->GetChild());
		getNodeHeights (p->GetSibling());
	}
}


//------------------------------------------------------------------------------
// Compute node depth (i.e, height above root). Based on COMPONENT 2.0 code, 
// assumes count is set to 0 prior to calling code
void Tree::getNodeDepth(NodePtr p)
{
	if (p)
    {
		
    	p->SetDepth (count);
    	count++;
        getNodeDepth (p->GetChild());
		count--;
        getNodeDepth (p->GetSibling());
    }
}

//------------------------------------------------------------------------------
// Compute node depth (i.e, height above root).
void Tree::GetNodeDepths ()
{
	count = 0;
	getNodeDepth (Root);
}


//------------------------------------------------------------------------------
void Tree::markNodes(NodePtr p, bool on)
{
	if (p)
	{
		p->SetMarked (on);
		markNodes (p->GetChild(), on);
		markNodes (p->GetSibling(), on);
	}
}

//------------------------------------------------------------------------------
void Tree::MarkNodes (bool on)
{
	markNodes (Root, on);
}

//------------------------------------------------------------------------------
void Tree::makeNodeList (NodePtr p)
{
	if (p)
	{
		makeNodeList (p->GetChild ());
		makeNodeList (p->GetSibling ());
		if (p->IsLeaf())
		{
			int leafnumberposition=p->GetLeafNumber()-1;//modified by BCO
			string plabel=p->GetLabel(); //modified by BCO
			LeafList[plabel] = leafnumberposition; //modified by BCO
			assert((Leaves+Internals)>leafnumberposition);
			Nodes[leafnumberposition] = p; //modified by BCO
			p->SetIndex (leafnumberposition); //modified by BCO
		}
		else
		{
			Nodes[count] = p;
			p->SetIndex (count);
			count++;
		}
		if (p != Root)
		{
		}
	}
}

//------------------------------------------------------------------------------
void Tree::MakeNodeList ()
{
	Update();//to ensure number of leaves and internals is correct //Added by BCO
	if (Nodes!=NULL) { //Added by BCO 
		delete [] Nodes; //Added by BCO
		Nodes=NULL; //Added by BCO
	} //Added by BCO
//	if (Nodes == NULL)// deleted by BCO: we will always delete the old node arrary and make a new one
//	{
		Nodes = new NodePtr [Leaves + Internals];
//	}
	count = Leaves;
	makeNodeList (Root);
}

//------------------------------------------------------------------------------
NodePtr Tree::GetLeafWithLabel (string s)
{
	NodePtr result = NULL;
	
	map<string, int, less<string> >::iterator index  = LeafList.find (s);
	if (index != LeafList.end())
		result = Nodes[LeafList[s]];
	
	return result;
}

//------------------------------------------------------------------------------
// Add Node below Below. Doesn't update any clusters, weights, etc.
void Tree::AddNodeBelow (NodePtr Node, NodePtr Below)
{
	NodePtr Ancestor = NewNode ();
	Ancestor->SetChild (Node);
	Node->SetAnc (Ancestor);
	NodePtr q = Below->GetAnc ();
	Internals++;
	if (Node->IsLeaf())
		Leaves++;
	if (q == NULL || Below == q->GetChild())
	{
		Node->SetSibling (Below);
		Ancestor->SetAnc (q);
		Ancestor->SetSibling (Below->GetSibling());
		Below->SetSibling (NULL);
		Below->SetAnc (Ancestor);
		if (q == NULL)
			Root = Ancestor;
		else
			q->SetChild (Ancestor);
	}
	else
	{
		// Get left sibling of Below
		NodePtr r = Below->LeftSiblingOf();
		while (Below != r->GetSibling())
			r = r->GetSibling();
		Node->SetSibling (Below);
		Ancestor->SetAnc (q);
		Ancestor->SetSibling (Below->GetSibling());
		Below->SetSibling (NULL);
		Below->SetAnc (Ancestor);
		r->SetSibling (Ancestor);
	}
}

//------------------------------------------------------------------------------
// Compute nodal heights based on path length from root, and store maximum
// value in plot.maxheight. Used by drawing routines.
void Tree::getPathLengths (NodePtr p)
{
	if (p)
	{
		if (p != Root)
		{
			float l = p->GetEdgeLength();
			if (l < 0.000001) // suppress negative branch lengths
				l = 0.0;
			p->SetPathLength (p->GetAnc()->GetPathLength() + l);
		}
		if (p->GetPathLength() > MaxPathLength)
			MaxPathLength = p->GetPathLength();
		getPathLengths (p->GetChild());
		getPathLengths (p->GetSibling());
	}
}

//------------------------------------------------------------------------------
// Fill in weight, degree, etc.
void Tree::buildtraverse (NodePtr p)
{
	if (p)
	{
		p->SetWeight (0);
		/*p->SetModelCategory(vector<double>(1,1.0)); //added by BCO
		p->SetStateOrder(vector<int>(1,0)); //Added by BCO
		p->SetStateTimes(vector<double>(1,0.0)); //Added by BCO*/
		p->SetDegree (0);
		buildtraverse (p->GetChild ());
		buildtraverse (p->GetSibling ());
		if (p->IsLeaf())
		{
			Leaves++;
			p->SetWeight (1);
			/*p->SetModelCategory(vector<double>(1,1.0)); //Added by BCO
			p->SetStateOrder(vector<int>(1,0)); //Added by BCO
			p->SetStateTimes(vector<double>(1,0.0)); //Added by BCO*/
			
		}
		else
		{
			Internals++;
		}
		if (p != Root)
		{
			p->GetAnc()->AddWeight (p->GetWeight());
			p->GetAnc()->IncrementDegree();
		}
	}
}

//------------------------------------------------------------------------------
// Ensure fields like weight, degree, etc are correct
void Tree::Update ()
{
	count = 0;
	Leaves = Internals = 0;
	buildtraverse (Root);
}


NodePtr Tree::RemoveNode (NodePtr Node)
{
	NodePtr result = NULL;
	
	if (Node == Root)
	{
		if (Leaves == 1)
		{
			Root = NULL;
			Node->SetAnc (NULL);
			Leaves = Internals = 0;
		}					
		return result;
	}
	
	NodePtr p;
	NodePtr Ancestor = Node->GetAnc();
	
	if (Ancestor->GetDegree() == 2)
	{
		// ancestor is binary, so remove node and its ancestor
		if (Node->IsTheChild ())
			p = Node->GetSibling();
		else
			p = Ancestor->GetChild();
		NodePtr q = Ancestor->GetAnc();
		p->SetAnc (q);
		if (q != NULL)
		{
			if (Ancestor->IsTheChild())
				q->SetChild (p);
			else
			{
				NodePtr r = Ancestor->LeftSiblingOf ();
				r->SetSibling (p);
			}
			p->SetSibling (Ancestor->GetSibling());
			result = p;
		}
		else
		{
			// Ancestor is the root
			Root = p;
			p->SetSibling (NULL);
			result = p;
		}
		delete Ancestor;
		Internals--;
		if (Node->IsLeaf())
			Leaves--;
		Node->SetAnc (NULL);
		Node->SetSibling (NULL);
	}
	else
	{
		// polytomy, just remove node
		NodePtr q;
		if (Node->IsTheChild())
		{
			Ancestor->SetChild (Node->GetSibling());
			q = Node->GetSibling ();
		}
		else
		{
			q = Node->LeftSiblingOf ();
			q->SetSibling (Node->GetSibling ());
		}
		Node->SetSibling (NULL);
		Node->SetAnc (NULL);
		if (Node->IsLeaf())
			Leaves--;
		Ancestor->SetDegree (Ancestor->GetDegree() - 1);
		result = q;
	}
}


//BCO function for getting the pointer to the node of the mrca
// for a pair of taxa.
//Node *GetMrcaPtr(Node *a, Node *b)
//{
//    NodePtr aanc;
//    aanc=a;
//    NodePtr banc;
//    NodePtr mrcaptr;
//    NodePtr Root;
//    bool mrcanotfound=true;
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
//}
	



//BCO added this function. It gets path distance
//from the given node to the root.
//double ReturnNodeDepth (NodePtr p)
//{
//    double totallength;
//    totallength=0;
//    if (p)
//    {
//        while (p != Root)
//        {
//            double l = p->GetEdgeLength();
//            totallength+=l;
//            p=p->GetAnc();
//        }
//    }
//    return totallength;
//}
	




