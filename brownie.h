
#ifndef __BROWNIE_H
#define __BROWNIE_H

#define COMMAND_MAXLEN  1000000// used to be 255
#define PI 3.141592653589793238462643383279502884197169399375
#define maxModelCategoryStates         10
#define BROWNIE_EPSILON 0.00001
#define BROWNIE_MAXLIKELIHOOD 1000000000 //Big but not big enough to blow up numerical optimization (I think).
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "containingtree.h"
#include "charactersblock2.h"



class BROWNIE : public NexusBlock, public Nexus
{
    friend class OptimizationFn;
    friend class OptimizationFnMultiModel;
    friend class ContainingTree;
	friend class CDFvectorholder;
public:
        bool inf_open;
    bool logf_open;
    bool echof_open;
    bool quit_now;
	bool quit_onerr;
    ofstream logf;
    ofstream echof;
    ofstream exportf;
    nxsstring echofname;
    nxsstring logfname;
    nxsstring message;
	nxsstring jackknifetreestooutput;
    unsigned long int gslseed;
    int gslseedtoprint;
    char next_command[COMMAND_MAXLEN+1];
    int maxnumspecies;
    int minnumspecies;
	int minsamplesperspecies;
	int maxstartstops;
	int rearrlimit;
    bool steepest;
	bool exhaustive;
    bool status;
	bool jackknifesearch;
    TreesBlock* trees;
    TaxaBlock* taxa;
    AssumptionsBlock* assumptions;
    CharactersBlock* characters;
	CharactersBlock2* characters2;
	CharactersBlock* discretecharacters;
	bool discretecharloaded;
	CharactersBlock* continuouscharacters;
    int chosentree;
    int chosenchar;
	int discretechosenchar;
	int numberoffreeparameters;
	int numberoffreerates;
	int numberoffreefreqs;
	vector<double> ratematfixedvector; 
	vector<int> ratematassignvector; 
	vector<double> userstatefreqvector;
	string freerateletterstring;
	nxsstring usermatrix;
	int negbounceparam; //equals -1, or the index of the parameter to have a negative (and therefore out of bounds) value
	bool nonnegvariables;
	gsl_matrix *optimaldiscretecharQmatrix;
	gsl_vector *optimaldiscretecharstatefreq;
	gsl_matrix *currentdiscretecharQmatrix;
	gsl_vector *currentdiscretecharstatefreq;	
	int discretechosenmodel;
/*	double E_discretegeneral(void *xp);
	double M_discretegeneral(void *xp, void *yp);
	void S_discretegeneral(const gsl_rng * r, void *xp, double step_size);*/
	int discretechosenstatefreqmodel;
	int numbercharstates;
	int localnumbercharstates;
	bool allchar;
	bool globalstates;
	bool variablecharonly;
	double bestdiscretelikelihood;
	int optimizationalgorithm;
    bool matrixsingular;
    bool debugmode;
    double structwt;
	bool triplettoohigh;
	bool gtptoohigh;
	bool infinitescore;
	double tripletdistthreshold;
    double pthreshold;
	double chosensubsampling;
    vector<double> movefreqvector;
    int unrooted;
    int maxiterations;
    double stoppingprecision;
    int randomstarts;
	nxsstring treefilename;
	bool useCOAL;
	bool useMS;
	int msbasereps;
	int contourBrlenToExport;
	int contourMaxRecursions;
	double contourstartingnumbersteps;
	double contourstartingwidth;
	bool exportalltrees;
	int COALaicmode;
	double markedmultiplier;
	double brlensigma;
	int numbrlenadjustments;
	map< vector<int>, int> quartetcounts;
	map< vector<int>, int> qualifiedquartets;
	vector<int> quartetspertaxon;
    double stepsize;
    bool detailedoutput;
	bool redobad;
	int giveupfactor;
    bool citationarray[20];
    vector<int> convertsamplestospecies; //Indexed by sample number, each entry is the species to which it's assigned
	vector<int> jackknifevector;
	vector<int> geneidvector;
	vector<vector<int> > BestConversions; //Vector of the best convertsamplestospeciesvectors
	vector<double> TotalScores;
	vector<double> StructScores;
	vector<double> GTPScores;
    bool sppnumfixed;
    bool outputallatonce;
    vector<int> staterestrictionvector;
    vector<double> timeslicetimes;
    vector<int> timeslicemodels;
    int chosenmodel;
	gsl_vector *optimalvaluescontinuouschar; //first values are parameters; last is lnL
	gsl_matrix *optimalVCV;
	//vector<nxsstring> optimalvalueslabels;
	gsl_vector *optimalTraitMeans;
	nxsstring globalchosentaxset;
    map<nxsstring, IntSet> MrcaMap;
    Profile<Tree> intrees;
    int tipvariancetype;
    int progressbartotal;
    int progressbarcount;
    int progressbarprinted;
	gsl_matrix *TaxonDistance; //used for taxon-taxon distances, measured in terms of triplets
	gsl_matrix *TaxonProportDistance;
	map<nxsstring,int> TripletCounts; //If triplet is (5,3),6, the triplet string is 3_5_6 (first two digits are the closest pair, ordered with the smaller first). If (6,5),3, the string is 5_6_3
	map<nxsstring,int> TripleCounts; //List of number of times each triple (set of three taxa, regardless of topology) appears: this may different for different sets of taxa, due to some taxa being missing from some treees 
	ContainingTree ComputeTripletNJTree();
	virtual void InitializeQuartetCounts();
	vector<nxsstring> ReturnClade(Tree *T, nxsstring a, nxsstring b, nxsstring c);
	virtual void DoNast();
	virtual void HandleNast(NexusToken&);
	double npercent;
    int nreps;
	int jreps;
	int jackrep;
	double pctdelete;
    bool showtries;
    vector<ContainingTree> RawBestTrees; //Contains all the best trees with just species labels; cleared every time a better tree is found
    vector<ContainingTree> FormattedBestTrees; //Uses samples as labels. Is useful as this way you don't have to remember the convertsamplestospecies for each best tree
	vector<nxsstring> ContourSearchDescription; //Stores a tree description and tab-delimited file with branch lengths and likelihood scores
	vector<vector <double> > ContourSearchVector;
	vector<ContainingTree> BestBranchlengthTreeForThisNextTree; //For optimizing branch lengths
    double bestscore;
	double bestscorelocal;
        //We use a struct so we can pass one set of params (a VCV matrix and two vectors) to a multimin function
    struct MatrixVectorVector {
        gsl_matrix *Matrix1;
        gsl_vector *Vector1;
        gsl_vector *Vector2;
    };
		//This rather clunky struct is for passing a set of VCV matrices (one per model type)
		//If we increase the number of models, we should also increase the number of matrices passed (or pass one thing more intelligently)
    struct MultiMatrixVectorVector {
        gsl_matrix *Matrix0;
        gsl_matrix *Matrix1;
        gsl_matrix *Matrix2;
        gsl_matrix *Matrix3;
        gsl_matrix *Matrix4;
        gsl_matrix *Matrix5;
        gsl_matrix *Matrix6;
        gsl_matrix *Matrix7;
        gsl_matrix *Matrix8;
        gsl_matrix *Matrix9;
        gsl_vector *Vector1;
        gsl_vector *Vector2;
    };

public:
        map<string, double> SimulateBrownian(double trend,double rate,double rootstate);
	double browniesafe_gsl_sf_exp(double x);
    void PreOrderTraversal( NexusToken& token);
    double GetTripletScore(ContainingTree *SpeciesTreePtr);
	void GetTaxonTaxonTripletDistances();
    vector<int> GetTripletOverlap(ContainingTree *t1, ContainingTree *t2, int taxaincommon);
	void DelDupes();
	void HandleAccuracy( NexusToken& token );
	void ComputeAccuracy();
	void HandlePartitionedEdgeSupport ( NexusToken& token);
	void BatchPartitionedEdgeSupport (int numberofpartitions);
        double ComputeTripletCost(int numberdisagree,int maxnumber,int ntaxincommon,double Tree1Wt,int Tree1Ntax,double Tree2Wt,int Tree2Ntax, int numberofgenes);
    void HandleExport ( NexusToken& token );
    void ProgressBar(int total);
    int CharLabelToNumber( nxsstring s );
    bool FileExists( const char* fn );
    nxsstring GetFileName( NexusToken& token );
    nxsstring GetNumber( NexusToken& token );
    nxsstring GetNumberOnly( NexusToken& token );
    nxsstring& blanks_to_underscores( nxsstring& s );
    nxsstring& underscores_to_blanks( nxsstring& s );
    void FactoryDefaults();
    void HandleEndblock( NexusToken& token );
    void HandleBlocks( NexusToken& token );
    void HandleDebug( NexusToken& token );
	void HandleNoQuitOnErr( NexusToken& token );
	void HandleQuitOnErr( NexusToken& token );
	double AIC(double neglnL, int K);
	double AICc(double neglnL, int K, int N);
    void HandleDebugOptimization (NexusToken & token);
    void HandleHelp( NexusToken& token );
    void HandleLog( NexusToken& token );
    void HandleEcho( NexusToken& token );
    void HandleRateTest( NexusToken& token);
	void HandleOrderByTree( NexusToken& token);
    void HandleExecute( NexusToken& token );
    void HandleChoose( NexusToken& token );
    void HandleSet( NexusToken& token );
    void HandleModel( NexusToken& token );
	void HandleDiscrete( NexusToken& token );
	void HandleSimulateCharacters( NexusToken& token);
	void SimulateCharacters(int n, int chartype, nxsstring outputfilename, bool treeloop);
	void HandlePagelDiscrete( NexusToken& token);
	void HandleLoss( NexusToken& token );
	void FindFixedDiscreteModel();
    void NumOpt( NexusToken& token);
	void DoExhaustiveSearch();
    void DoHeuristicSearch();
	void HandleDettmanCollapse( NexusToken& token );
	void DoDettmanCollapse();
    void HandleHeuristicSearch( NexusToken& token );
    void HandleShowtree( NexusToken& token );
    void HandleTipValues(NexusToken& token);
    void HandlePrintEdgeLengths(NexusToken& token);
    void HandleVCV(NexusToken& token);
	virtual double CalculateDiscreteCharLnL(gsl_matrix * RateMatrix, gsl_vector * ancestralstatevector);
	virtual double CalculateDiscreteCharLnLHetero(gsl_matrix * RateMatrixHetero, gsl_vector * ancestralstatevector);
	virtual double CalculateDiscreteCharProbAllConstant(gsl_matrix * RateMatrix, gsl_vector * ancestralstatevector);
	virtual NodePtr EstimateMLDiscreteCharJointAncestralStates(gsl_matrix * RateMatrix, gsl_vector * ancestralstatevector, int breaksperbranch);
	virtual double CalculateDiscreteLindy2(double rateA, double rateB);
	virtual double CalculateDiscreteLindy1(double rateA);
	static double GetLikelihoodUnderLindy2_gsl( const gsl_vector * variables, void *obj) ;
	static double GetDiscreteCharLnL_gsl( const gsl_vector * variables, void *obj) ;
	double GetDiscreteCharLnL(const gsl_vector * variables);

	double GetLikelihoodUnderLindy2(const gsl_vector * variables);
	static double GetLikelihoodUnderLindy1_gsl( const gsl_vector * variables, void *obj) ;
	double GetLikelihoodUnderLindy1(const gsl_vector * variables);
	gsl_vector * LindyGeneralOptimization(int ChosenModel);
    gsl_vector* DiscreteGeneralOptimization();	
	gsl_matrix* ComputeTransitionProb(gsl_matrix *RateMatrix, double brlen);
    void HandleTimeSlice( NexusToken& token );
    void HandleTipVariance( NexusToken& token );
    void HandleTestLikelihoodWithTipVariance(NexusToken& token);
	void Assign(NexusToken& token);
    void HandleCitation( NexusToken& token );
	void PrintCitations();
    void OutputComment( nxsstring& token );
    void HandleGarland( NexusToken& token );
    void HandleMRCA ( NexusToken& token );
    void HandleGrepCount(NexusToken& token);
	virtual nxsstring PipeLeafGTP ();
	virtual nxsstring PipeLeafSpeciesTree ();
	virtual nxsstring PipeGTP (Tree intree);
	virtual nxsstring PipeSpeciesTree (ContainingTree *SpeciesTreePtr);
	virtual nxsstring PipeEndOfTreeGTP ();
	virtual nxsstring PipeEndOfTreeSpeciesTree ();	
	virtual nxsstring PipeLeftParenthesis();
	virtual nxsstring PipeRightParenthesis();
	virtual nxsstring PipeSiblingSymbol();
	virtual nxsstring PipeInternal();
        virtual nxsstring ReturnFinalSpeciesTree (Tree t);
        virtual int PruneToOverlappingLeaves(ContainingTree *t1, ContainingTree *t2);
        virtual int PrepareTreesForTriplet(ContainingTree *t1, ContainingTree *t2);
        virtual void HandleCompareRandomTrees(NexusToken& token);
        virtual nxsstring PipeLeafFinalSpeciesTree();
        virtual vector<ContainingTree> MakePrettyForOutput();
        virtual void FormatAndStoreBestTree(ContainingTree *NewBestTree,vector<double> scorevector);
        int badgtpcount;
        vector<vector<double> > CDFvector;
        vector<vector<int> > CladeVector;
		vector<double> CladeVectorNJBrlen; //Stores initial NJ brlen
		vector<double> CladeVectorTripletSupport; //Stores proportion (0-1) of triplets relevant to the clade which support that clade
		double meaninternalbrlen; //stores the mean internal brlen of the starting nj tree, used for deciding on initial splits
        vector<int> SamplesToMove;
        vector<int> SampleDestinations;
        virtual bool MoveSamples(vector<int> Originalconvertsamplestospecies);
        virtual bool TestMoveSamples(int Sample, int Destination);
        virtual bool CheckConvertSamplesToSpeciesVector(bool fix);
		virtual bool CombineSpeciesWithTooFewSamples(bool fix);
		virtual bool CheckConvertSamplesToSpeciesVectorSpNum(int actualmaxspeciesnum);
		virtual bool CheckConvertSamplesToSpeciesTooManySpecies();
		virtual bool CheckConvertSamplesToSpeciesTooFewSpecies();
		virtual double DoAllAssignments(double bestscore, int maxspecies, ContainingTree *SpeciesTree );
     //   virtual double GetGTPScore(ContainingTree *SpeciesTreePtr); //took out as no longer use external GTP
        virtual double GetGTPScoreNew(ContainingTree *SpeciesTreePtr);
        virtual vector<double> GetCombinedScore(ContainingTree *SpeciesTreePtr);
	Node *cur;
	std::stack < Node *, std::vector<Node *> > stk;
    void PurgeBlocks();
    virtual void Read( NexusToken& token );
    virtual void Reset();
    int TaxonLabelToNumber( nxsstring s );
    gsl_vector* GetTipValues(nxsstring chosentaxset, int charnumber);
    gsl_vector* SimulateTips(gsl_matrix * VCV, double rate, gsl_vector *MeanValues);
	virtual void GetOptimalVCVAndTraitsContinuous();
    gsl_matrix* GetVCV(nxsstring chosentaxset);
	void PrintMatrix(gsl_matrix *VCV);
	void PrintVector(gsl_vector *somevector);
	gsl_matrix* GetVCVwithKappa(nxsstring chosentaxset,double kappa); //don't forget to use DeleteStem
	gsl_matrix* ConvertVCVwithDelta(gsl_matrix * VCVorig,double delta); //takes VCV as input; could use DeleteStem(GetVCV(chosentaxset)) as input
	gsl_matrix* ConvertVCVwithLambda(gsl_matrix * VCVorig,double lambda);//takes VCV as input; could use DeleteStem(GetVCV(chosentaxset)) as input
	gsl_matrix* GetVCVwithTree(nxsstring chosentaxset, Tree t);
    gsl_matrix* GetVCVforOneModel(nxsstring chosentaxset, int selectedmodel);
    gsl_matrix* GetStartStopTimesforOneState(nxsstring chosentaxset, int selectedstate);
    gsl_matrix* GetVCVforChangeNoChange(nxsstring chosentaxset, bool wantchangeedges);
	//gsl_matrix* GetVCVforIntervalCrossing(nxsstring chosentaxset, bool wantchangeedges);
    gsl_matrix* DeleteStem(gsl_matrix * VCVorig);
    gsl_matrix* AddTipVarianceVectorToRateTimesVCV(gsl_matrix * VCVorig,gsl_vector * TipVariance);
    double GetLikelihoodWithGivenTipVariance(const gsl_vector * variables, void *params);
	char* OutputForGTP(ContainingTree *SpeciesTreePtr);
    gsl_rng * ReturnR();
       // double GetLikelihoodWithGivenTipVariance(const gsl_vector * variables, void *params);
    
      //double Wrapper_For_GetLikelihoodWithGivenTipVariance(const gsl_vector * variables, void *params);
    void MakeCombinedVCV(gsl_matrix *VCVcombined, gsl_matrix *VCVtoadd, int ntaxprocessed);
    void MakeCombinedTips(gsl_vector *tipscombined, gsl_vector *tipstoadd, int ntaxprocessed);
    double GetLScore(gsl_matrix * VCV,gsl_vector *tipresid,double rate);
       // gsl_vector* ConvertVCVMatrixToVector(gsl_matrix * VCV, gsl_vector *OrigVector);
       // gsl_matrix* ConvertVCVVectorToMatrix(gsl_vector *VCVvector);
       // gsl_matrix* ExtractMatrixFromVector(gsl_vector *VCVvector, int ntaxprocessed, int currentntax);
    
      // gsl_matrix amoeba(gsl_matrix startingsimplex,int chosenFn,int repnumber,gsl_matrix likelihoodAndOutputParametersMatrix,double toler,int numberofparameters,int maxsteps);
      //  double amotry(gsl_matrix **startingsimplex,gsl_vector y, gsl_vector psum, int numberofparameters, int chosenFn,int ihi, double fac);
      //   double GetFunctionValue(int chosenFn,gsl_vector  inputvector);


    
//        double GetAncestralState(gsl_matrix VCV, gsl_vector tips);

    
       // Node *GetMrcaPtr (Node *a, Node *b);

public:
        BROWNIE();
    ~BROWNIE();
    void HandleExecuteCmdLine(nxsstring fn );
    double GetAncestralState(gsl_matrix *VCV, gsl_vector *tips); //should be protected?
    gsl_vector* GetTipResiduals(gsl_vector *tips, double ancestralstate);
    double EstimateRate(gsl_matrix *VCV, gsl_vector *tipresiduals);
    void HandleGettrees( NexusToken& token );
    void EnteringBlock( nxsstring blockName );
    void ExitingBlock( nxsstring blockName );
    void ExecuteStarting() {}
    void ExecuteStopping() {}
    void OutputComment( nxsstring ) {}
    void HandleNextCommand();
    void NexusError( nxsstring& msg, streampos pos, long line, long col );
    void PreprocessNextCommand();
    void PrintMessage( bool linefeed = true );
    virtual void Report( ostream& out );
	//void Run();
    void RunCmdLine(bool inputfilegiven, nxsstring fn);
    void SkippingBlock( nxsstring blockName );
    void SkippingCommand( nxsstring commandName );
    void SkippingDisabledBlock( nxsstring blockName );
    virtual bool UserSaysOk( nxsstring mb_message, nxsstring mb_title );
    nxsstring GetTaxonLabel( int i );
           //gsl_vector MatrixTimesVector(gsl_matrix A, gsl_vector b);
           // gsl_vector NMsimplex(int chosenFn,gsl_matrix inputparameters,int numberofstarts,int maxsteps, bool verbose);
    
};

#endif

