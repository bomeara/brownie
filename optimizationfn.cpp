/*
 *  optimizationfns.cpp
 *
 *
 *  Created by Brian O'Meara on 4/7/06.
 *  GPL2
 *
 */
#include <strstream>
#include <fstream.h>
#include <iomanip.h>
#include <unistd.h>
#include <stdio.h>
#include <set>
#include <math.h>

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
#include "gport.h"
#include "profile.h"
#include "ntree.h"
#include "nodeiterator.h"
#include "treeorder.h"
#include "treedrawer.h"
#include "TreeLib.h"
#include "gtree.h"
#include "treereader.h"
#include "treewriter.h"
#include <time.h>
#include <map>
#include "brownie.h"
#include "optimizationfn.h"
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
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>

//constructor
OptimizationFnMultiModel::OptimizationFnMultiModel( gsl_matrix *InMatrix0,  gsl_matrix *InMatrix1,  gsl_matrix *InMatrix2,  gsl_matrix *InMatrix3,  gsl_matrix *InMatrix4,  gsl_matrix *InMatrix5,  gsl_matrix *InMatrix6,  gsl_matrix *InMatrix7,  gsl_matrix *InMatrix8,  gsl_matrix *InMatrix9, gsl_vector *InVector1,gsl_vector *InVector2, int Inmaxiterations, double Instoppingprecision, int Inrandomstarts, double Instepsize, bool Indetailedoutput) :brownie()
{
	//gsl_set_error_handler_off();
    //id = "OptimizationFn";
    int ntax=InMatrix0->size1;
	Matrix0=gsl_matrix_calloc(InMatrix0->size1,InMatrix0->size2);
    Matrix1=gsl_matrix_calloc(InMatrix1->size1,InMatrix1->size2);
	Matrix2=gsl_matrix_calloc(InMatrix2->size1,InMatrix2->size2);
	Matrix3=gsl_matrix_calloc(InMatrix3->size1,InMatrix3->size2);
    Matrix4=gsl_matrix_calloc(InMatrix4->size1,InMatrix4->size2);
	Matrix5=gsl_matrix_calloc(InMatrix5->size1,InMatrix5->size2);
	Matrix6=gsl_matrix_calloc(InMatrix6->size1,InMatrix6->size2);
    Matrix7=gsl_matrix_calloc(InMatrix7->size1,InMatrix7->size2);
	Matrix8=gsl_matrix_calloc(InMatrix8->size1,InMatrix8->size2);
	Matrix9=gsl_matrix_calloc(InMatrix9->size1,InMatrix9->size2);
    Vector1=gsl_vector_calloc(ntax);
    Vector2=gsl_vector_calloc(ntax);
    int CopyResult= gsl_matrix_memcpy(Matrix0, InMatrix0);
	CopyResult= gsl_matrix_memcpy(Matrix1, InMatrix1);
	CopyResult= gsl_matrix_memcpy(Matrix2, InMatrix2);
	CopyResult= gsl_matrix_memcpy(Matrix3, InMatrix3);
	CopyResult= gsl_matrix_memcpy(Matrix4, InMatrix4);
    CopyResult= gsl_matrix_memcpy(Matrix5, InMatrix5);
	CopyResult= gsl_matrix_memcpy(Matrix6, InMatrix6);
	CopyResult= gsl_matrix_memcpy(Matrix7, InMatrix7);
	CopyResult= gsl_matrix_memcpy(Matrix8, InMatrix8);
	CopyResult= gsl_matrix_memcpy(Matrix9, InMatrix9);	
	CopyResult=gsl_vector_memcpy(Vector1, InVector1);
    CopyResult=gsl_vector_memcpy(Vector2, InVector2);
    maxiterations=Inmaxiterations;
    stoppingprecision=Instoppingprecision;
    randomstarts=Inrandomstarts;
    stepsize=Instepsize;
    detailedoutput=Indetailedoutput;
	fixedparams=gsl_vector_calloc(1);
	// cout<<"First entry in Matrix1 is "<<gsl_matrix_get(Matrix1,0,0)<<endl;
}

OptimizationFnMultiModel::~OptimizationFnMultiModel()
{
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
	gsl_vector_free(fixedparams);
}


//constructor
OptimizationFn::OptimizationFn( gsl_matrix *InMatrix1, gsl_vector *InVector1,gsl_vector *InVector2, int Inmaxiterations, double Instoppingprecision, int Inrandomstarts, double Instepsize, bool Indetailedoutput) :brownie()
{
	//gsl_set_error_handler_off();
    //id = "OptimizationFn";
    int ntax=InMatrix1->size1;
    Matrix1=gsl_matrix_calloc(ntax,ntax);
    Vector1=gsl_vector_calloc(ntax);
    Vector2=gsl_vector_calloc(ntax);
    int CopyResult= gsl_matrix_memcpy(Matrix1, InMatrix1);
    CopyResult=gsl_vector_memcpy(Vector1, InVector1);
    CopyResult=gsl_vector_memcpy(Vector2, InVector2);
    maxiterations=Inmaxiterations;
    stoppingprecision=Instoppingprecision;
    randomstarts=Inrandomstarts;
    stepsize=Instepsize;
    detailedoutput=Indetailedoutput;
	// cout<<"First entry in Matrix1 is "<<gsl_matrix_get(Matrix1,0,0)<<endl;
}

/**
* @destructor
 *
 */
OptimizationFn::~OptimizationFn()
{
	gsl_matrix_free(Matrix1);
	gsl_vector_free(Vector1);
	gsl_vector_free(Vector2);

}
	

//constructor
//OptimizationFnDiscreteRate::OptimizationFnDiscreteRate(Tree *InBrownieTree, gsl_matrix *InMatrix0,  gsl_matrix *InMatrix1,  gsl_matrix *InMatrix2,  gsl_matrix *InMatrix3,  gsl_matrix *InMatrix4,  gsl_matrix *InMatrix5,  gsl_matrix *InMatrix6,  gsl_matrix *InMatrix7,  gsl_matrix *InMatrix8,  gsl_matrix *InMatrix9, gsl_vector *InVector1, int Inmaxiterations, double Instoppingprecision, int Inrandomstarts, double Instepsize, bool Indetailedoutput, int Instatefrequencychoice) :brownie()
//{
 //   int nstates=InMatrix0->size1;
//	Tree *BrownieTree=InBrownieTree;
//	Matrix0=gsl_matrix_calloc(InMatrix0->size1,InMatrix0->size2);
//    Matrix1=gsl_matrix_calloc(InMatrix1->size1,InMatrix1->size2);
//	Matrix2=gsl_matrix_calloc(InMatrix2->size1,InMatrix2->size2);
//	Matrix3=gsl_matrix_calloc(InMatrix3->size1,InMatrix3->size2);
  //  Matrix4=gsl_matrix_calloc(InMatrix4->size1,InMatrix4->size2);
//	Matrix5=gsl_matrix_calloc(InMatrix5->size1,InMatrix5->size2);
//	Matrix6=gsl_matrix_calloc(InMatrix6->size1,InMatrix6->size2);
//    Matrix7=gsl_matrix_calloc(InMatrix7->size1,InMatrix7->size2);
//	Matrix8=gsl_matrix_calloc(InMatrix8->size1,InMatrix8->size2);
//	Matrix9=gsl_matrix_calloc(InMatrix9->size1,InMatrix9->size2);
//    Vector1=gsl_vector_calloc(ntax);
//    int CopyResult= gsl_matrix_memcpy(Matrix0, InMatrix0);
//	CopyResult= gsl_matrix_memcpy(Matrix1, InMatrix1);
//	CopyResult= gsl_matrix_memcpy(Matrix2, InMatrix2);
//	CopyResult= gsl_matrix_memcpy(Matrix3, InMatrix3);
//	CopyResult= gsl_matrix_memcpy(Matrix4, InMatrix4);
  //  CopyResult= gsl_matrix_memcpy(Matrix5, InMatrix5);
	//CopyResult= gsl_matrix_memcpy(Matrix6, InMatrix6);
//	CopyResult= gsl_matrix_memcpy(Matrix7, InMatrix7);
//	CopyResult= gsl_matrix_memcpy(Matrix8, InMatrix8);
//	CopyResult= gsl_matrix_memcpy(Matrix9, InMatrix9);	
//	CopyResult=gsl_vector_memcpy(Vector1, InVector1);
//    CopyResult=gsl_vector_memcpy(Vector2, InVector2);
//    maxiterations=Inmaxiterations;
//    stoppingprecision=Instoppingprecision;
//    randomstarts=Inrandomstarts;
//    stepsize=Instepsize;
//    detailedoutput=Indetailedoutput;
//	statefrequencychoice=Instatefrequencychoice;
	// cout<<"First entry in Matrix1 is "<<gsl_matrix_get(Matrix1,0,0)<<endl;
//}

//OptimizationFnDiscreteRate::~OptimizationFnDiscreteRate()
//{
//	gsl_matrix_free(Matrix0);
//	gsl_matrix_free(Matrix1);
//	gsl_matrix_free(Matrix2);
//	gsl_matrix_free(Matrix3);
//	gsl_matrix_free(Matrix4);
//	gsl_matrix_free(Matrix5);
//	gsl_matrix_free(Matrix6);
//	gsl_matrix_free(Matrix7);
//	gsl_matrix_free(Matrix8);
//	gsl_matrix_free(Matrix9);
//	gsl_vector_free(Vector1);
//}


double OptimizationFn::GetLikelihoodWithGivenTipVariance(const gsl_vector * variables) // where params points to the MatrixVectorVector struct //gsl_matrix * VCV, gsl_vector *tipvariance, gsl_vector * observedtips)
{
    //cout<<"Now in OptimizationFn::GetLikelihoodWithGivenTipVariance"<<endl;
    //cout<<"Size of variables vector = "<<variables->size<<endl;
    //gsl_matrix *VCV=((struct MatrixVectorVector *)params)->Matrix1;
    //gsl_vector *tipvariance=((struct MatrixVectorVector *)params)->Vector2;
    //gsl_vector *observedtips=((struct MatrixVectorVector *)params)->Vector1;
	double rate=gsl_vector_get(variables,0);
    //cout<<"Rate is "<<rate<<endl;
	int ntax=Matrix1->size1;
	//cout<<"Ntax = "<<ntax<<endl;
	gsl_matrix *VCV=gsl_matrix_calloc(ntax,ntax);
	gsl_vector *tipvariance=gsl_vector_calloc(ntax);
	gsl_vector *observedtips=gsl_vector_calloc(ntax);
	gsl_matrix_memcpy(VCV,Matrix1);
	gsl_vector_memcpy(tipvariance,Vector2);
	gsl_vector_memcpy(observedtips,Vector1);
	// cout<<"Observed state in first taxon = "<<gsl_vector_get(observedtips,0)<<endl;
	// cout<<"Tip variance in first taxon = "<<gsl_vector_get(tipvariance,0)<<endl;
	gsl_matrix *RateTimesVCV=gsl_matrix_calloc(ntax,ntax);
	gsl_matrix_memcpy(RateTimesVCV, VCV);
	gsl_matrix_scale(RateTimesVCV,rate);
	gsl_matrix *VCVfinal=gsl_matrix_calloc(ntax,ntax);
	
	gsl_matrix_add_constant (RateTimesVCV, -1.0*gsl_matrix_min (RateTimesVCV)); //replaces DeleteStem
	gsl_matrix_memcpy(VCVfinal,RateTimesVCV);
	for (int r=0;r<ntax;r++) {
        gsl_matrix_set(VCVfinal,r,r,(gsl_matrix_get(RateTimesVCV,r,r)+(gsl_vector_get(tipvariance,r)))); //replaces AddTipVarianceVectorToRateTImesVCV
    }
	
	//gsl_matrix_memcpy(VCVfinal,brownie.AddTipVarianceVectorToRateTimesVCV(brownie.DeleteStem(RateTimesVCV),tipvariance));
	gsl_vector *tipresiduals=gsl_vector_calloc(ntax);
    //cout<<"RateTimesVCV cell 0,0 ="<<gsl_matrix_get(RateTimesVCV,0,0)<<endl;
    //cout<<"VCVfinal cell 0,0 = "<<gsl_matrix_get(VCVfinal,0,0)<<endl;
    //cout<<"Size VCVfinal ="<<VCVfinal->size1<<" size of observedtips = "<<observedtips->size<<endl;
	double ancestralstate=brownie.GetAncestralState(VCVfinal,observedtips);
    //cout<<"Ancestral state = "<<ancestralstate<<endl;
	gsl_vector_memcpy(tipresiduals,brownie.GetTipResiduals(observedtips,ancestralstate));
	double likelihood=(brownie.GetLScore(VCVfinal,tipresiduals,1)); //returns -lnL
	if (gsl_vector_min(tipvariance)<0 || rate<0) {
		likelihood=GSL_POSINF;
	}
	
	//cout<<"rate = "<<rate<<" ancstate ="<<ancestralstate<<" likelihood = "<<likelihood<<endl;
	gsl_matrix_free(VCV);
	gsl_vector_free(tipvariance);
	gsl_vector_free(tipresiduals);
	gsl_vector_free(observedtips);
	gsl_matrix_free(RateTimesVCV);
	gsl_matrix_free(VCVfinal);
	return likelihood; //-lnL actually
}


double OptimizationFn::GetLikelihoodWithGivenTipVariance_gsl(const gsl_vector * variables, void *obj)
{
    //cout<<"Now in OptimizationFn::GetLikelihoodWithGivenTipVariance_gsl"<<endl;
    //cout<<"Size of variables vector = "<<variables->size<<endl;
    //This one is for passing to the GSL integration functions
	double temp;
	temp= ((OptimizationFn*)obj)->GetLikelihoodWithGivenTipVariance(variables);
	return temp;
}

gsl_vector * OptimizationFn::OptimizeRateWithGivenTipVariance()
{
	//cout<<"Float" <<sizeof(float)<<endl<<"Double "<<sizeof(double)<<endl<<"Long Double "<<sizeof(long double)<<endl;
	//cout<<"maxiterations = "<<maxiterations<<endl<<"step size = "<<stepsize<<endl;
	//gsl_vector * inputrate=gsl_vector_calloc(1);
	//gsl_vector_set(inputrate,0,1);
	double bestlikelihood=GSL_POSINF;
	int hitlimitscount=0;
	size_t np = 1;
	gsl_vector * results=gsl_vector_calloc(np);
	const gsl_rng_type * T_rng_type;
	gsl_rng * r_rng;
	gsl_rng_env_setup();
	T_rng_type = gsl_rng_default;
	r_rng = gsl_rng_alloc (T_rng_type);
	double rate=brownie.EstimateRate(Matrix1,brownie.GetTipResiduals(Vector1,brownie.GetAncestralState(Matrix1,Vector1)));
	double rateestimates[randomstarts];
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
        //cout<<"Now in OPtimizaRateWithGivenTipVariance in OptimizationFn"<<endl;
		x = gsl_vector_alloc (np);
		//gsl_vector_set (x,0,gsl_ran_exponential (r_rng,rate));
		//gsl_vector_set (x,0,gsl_ran_exponential (brownie.ReturnR(),rate));
		gsl_vector_set (x,0,gsl_ran_exponential (r,rate));
		double startingvalue=gsl_vector_get(x,0);
		// cout<<"Starting value = "<<gsl_vector_get(x,0)<<endl;
		OptimizationFn *pt;
		pt=(this);
		double (*F)(const gsl_vector *, void *);
		F = &OptimizationFn::GetLikelihoodWithGivenTipVariance_gsl;
		/* Initialize method and iterate */
		gsl_multimin_function minex_func;
		minex_func.f=*F;
		minex_func.params=pt;
		minex_func.n = np;
		s = gsl_multimin_fminimizer_alloc (T, np);
		gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
		do
		{
			// cout<<"Now on iteration "<<iter<<endl;
			iter++;
			status = gsl_multimin_fminimizer_iterate(s);
			if (status!=0) { //0 Means it's a success
				printf ("error: %s\n", gsl_strerror (status));
				break;
			}
			size = gsl_multimin_fminimizer_size (s);
			//status = gsl_multimin_test_size (size, 1e-2);
			status = gsl_multimin_test_size (size, stoppingprecision); //since we want more precision
			if (status == GSL_SUCCESS)
			{
				//printf ("converged to minimum at\n");
			}
			//printf ("%5d ", iter);
			//	for (i = 0; i < np; i++)
			//	{
			//printf ("%10.3e ", gsl_vector_get (s->x, i));
			//	}
			//printf ("f() = %7.3f size = %.3f\n", s->fval, size);
		}
		while (status == GSL_CONTINUE && iter < maxiterations);
		if (s->fval<bestlikelihood) {
			gsl_vector_memcpy(results,s->x);
			bestlikelihood=s->fval;
		}
		rateestimates[startnum]=gsl_vector_get(s->x,0); //Store each rate estimate
		if (iter==maxiterations) {
			hitlimits=true;
			hitlimitscount++;
		}
		brownie.message="Replicate ";
		brownie.message+=startnum+1;
		if (hitlimits) {
			brownie.message+=" **WARNING**";
		}
		brownie.message+="\n   Starting value = ";
		brownie.message+=startingvalue;
		brownie.message+="\n   NM iterations needed = ";
		int iterationsrequired=iter;
		brownie.message+=iterationsrequired;
		if (hitlimits) {
			brownie.message+=" **Max iterations hit; see WARNING below**";
		}
		brownie.message+="\n   LnL = ";
		//brownie.message+=s->fval;
		char outputstring[60];
		sprintf(outputstring,"%60.45f",-1*(s->fval));
		brownie.message+=outputstring;
		brownie.message+="\n   Rate = ";
		brownie.message+=gsl_vector_get(s->x,0);
		if (detailedoutput) {
			brownie.PrintMessage();
		}
		//cout<<"Rep "<<startnum+1<<" Iter "<<iter<<" LnL "<<s->fval<<" Rate "<<gsl_vector_get(s->x,0)<<endl;
		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free (s);
	}
	if (hitlimitscount>0) {
		brownie.message="\n----------------------------------------------------------------------------\n WARNING: Out of ";
		brownie.message+=randomstarts;
		brownie.message+=" optimization starts, ";
		brownie.message+=hitlimitscount;
		if (hitlimitscount==1) {
			brownie.message+=" was ";
		}
		else {
			brownie.message+=" were ";
		}
		brownie.message+="stopped by hitting\n  the maximum # of iterations. This means that those replicates\n  may not even have hit the local maximum.\n\n  You can increase the maximum number of iterations or decrease the\n  precision with the NumOpt command. You could also consider\n  increasing the number of random starts using that same command.\n\n  If this happened on a small proportion of replicates, though,\n  or if the precision (below) is good enough, don't worry about it.\n----------------------------------------------------------------------------";
		if (detailedoutput) {
			brownie.PrintMessage();
		}
		else if ((randomstarts-hitlimitscount)<10 && (hitlimitscount/randomstarts)>.1) {
			brownie.PrintMessage();
		}
	}
	brownie.message="\n\nRate ML estimate = ";
	brownie.message+=gsl_vector_get(results,0);
	brownie.message+="\nMean estimate across starts = ";
	brownie.message+=gsl_stats_mean(rateestimates,1,randomstarts);
	
	if (detailedoutput) {
		brownie.PrintMessage();
	}
	brownie.message="\nStandard deviation of estimates across ";
	brownie.message+=randomstarts;
	brownie.message+=" starts = ";
	brownie.message+=gsl_stats_sd(rateestimates,1,randomstarts);
	brownie.message+="\n[This is the precision of the rate estimate: numerical optimization does not give an exact value]";
	brownie.PrintMessage();
	return results;
};

double OptimizationFn::GetLikelihoodUnderACDC(const gsl_vector * variables)
{
	double rate=gsl_vector_get(variables,0);
	double ancestralstate=gsl_vector_get(variables,1);
	double g=gsl_vector_get(variables,2);
	int ntax=Matrix1->size1;
	gsl_matrix *VCV=gsl_matrix_calloc(ntax,ntax);
	gsl_vector *tipvariance=gsl_vector_calloc(ntax);
	gsl_vector *observedtips=gsl_vector_calloc(ntax);
	gsl_matrix_memcpy(VCV,Matrix1);
	gsl_vector_memcpy(tipvariance,Vector2);
	gsl_vector_memcpy(observedtips,Vector1);
	gsl_matrix *TransformedMatrix=gsl_matrix_calloc(ntax,ntax);
	for (int rowtaxon=0;rowtaxon<ntax;rowtaxon++) {
		for (int coltaxon=rowtaxon;coltaxon<ntax;coltaxon++) {
			double newVCVvalue=rate*(1-(pow(g,(-1*gsl_matrix_get(VCV,rowtaxon,coltaxon)))))/(1-(1/g));
			gsl_matrix_set(TransformedMatrix,rowtaxon,coltaxon,newVCVvalue);
			gsl_matrix_set(TransformedMatrix,coltaxon,rowtaxon,newVCVvalue);
		}
	}
	gsl_matrix *VCVfinal=gsl_matrix_calloc(ntax,ntax);
	gsl_matrix_add_constant (TransformedMatrix, -1.0*gsl_matrix_min (TransformedMatrix)); //replaces DeleteStem
	gsl_matrix_memcpy(VCVfinal,TransformedMatrix);
	for (int r=0;r<ntax;r++) {
        gsl_matrix_set(VCVfinal,r,r,(gsl_matrix_get(TransformedMatrix,r,r)+(gsl_vector_get(tipvariance,r)))); //replaces AddTipVarianceVectorToRateTImesVCV
    }
	
	//gsl_matrix_memcpy(VCVfinal,brownie.AddTipVarianceVectorToRateTimesVCV(brownie.DeleteStem(TransformedMatrix),tipvariance));
	gsl_vector *tipresiduals=gsl_vector_calloc(ntax);
	//cout<<"RateTimesVCV cell 0,0 ="<<gsl_matrix_get(RateTimesVCV,0,0)<<endl;
	//cout<<"VCVfinal cell 0,0 = "<<gsl_matrix_get(VCVfinal,0,0)<<endl;
	//cout<<"Size VCVfinal ="<<VCVfinal->size1<<" size of observedtips = "<<observedtips->size<<endl;
	// double ancestralstate=brownie.GetAncestralState(VCVfinal,observedtips);
	//cout<<"Ancestral state = "<<ancestralstate<<endl;
	gsl_vector_memcpy(tipresiduals,brownie.GetTipResiduals(observedtips,ancestralstate));
	double likelihood=(brownie.GetLScore(VCVfinal,tipresiduals,1)); //-lnL actually
	if (rate<0) {
		likelihood=GSL_POSINF;
	}
	//cout<<"likelihood "<<likelihood<<endl;
	gsl_matrix_free(VCV);
	gsl_vector_free(tipvariance);
	gsl_vector_free(observedtips);
	gsl_matrix_free(TransformedMatrix);
	gsl_matrix_free(VCVfinal);
	gsl_vector_free(tipresiduals);
	return likelihood; //-lnL actually
}


double OptimizationFn::GetLikelihoodUnderOU1(const gsl_vector * variables)
{
	double rate=gsl_vector_get(variables,0);
	double ancestralstate=gsl_vector_get(variables,1);
	double d=gsl_vector_get(variables,2);
	int ntax=Matrix1->size1;
	gsl_matrix *VCV=gsl_matrix_calloc(ntax,ntax);
	gsl_vector *tipvariance=gsl_vector_calloc(ntax);
	gsl_vector *observedtips=gsl_vector_calloc(ntax);
	gsl_matrix_memcpy(VCV,Matrix1);
	gsl_vector_memcpy(tipvariance,Vector2);
	gsl_vector_memcpy(observedtips,Vector1);
	gsl_matrix *TransformedMatrix=gsl_matrix_calloc(ntax,ntax);
	for (int rowtaxon=0;rowtaxon<ntax;rowtaxon++) {
		for (int coltaxon=rowtaxon;coltaxon<ntax;coltaxon++) {
			double newVCVvalue;
			if (coltaxon==rowtaxon) {
				newVCVvalue=(1-pow(d,2*gsl_matrix_get(VCV,rowtaxon,coltaxon)))*rate/(1-pow(d,2));
			}
			else {
				newVCVvalue=(pow(d,(gsl_matrix_get(VCV,coltaxon,coltaxon)+gsl_matrix_get(VCV,rowtaxon,rowtaxon)-(2*(gsl_matrix_get(VCV,rowtaxon,coltaxon))))))*(1-pow(d,2*gsl_matrix_get(VCV,rowtaxon,coltaxon)))*rate/(1-pow(d,2));
			}
			gsl_matrix_set(TransformedMatrix,rowtaxon,coltaxon,newVCVvalue);
			gsl_matrix_set(TransformedMatrix,coltaxon,rowtaxon,newVCVvalue);
		}
	}
	gsl_matrix *VCVfinal=gsl_matrix_calloc(ntax,ntax);
	gsl_matrix_add_constant (TransformedMatrix, -1.0*gsl_matrix_min (TransformedMatrix)); //replaces DeleteStem
	gsl_matrix_memcpy(VCVfinal,TransformedMatrix);
	for (int r=0;r<ntax;r++) {
        gsl_matrix_set(VCVfinal,r,r,(gsl_matrix_get(TransformedMatrix,r,r)+(gsl_vector_get(tipvariance,r)))); //replaces AddTipVarianceVectorToRateTImesVCV
    }
	
	//gsl_matrix_memcpy(VCVfinal,brownie.AddTipVarianceVectorToRateTimesVCV(brownie.DeleteStem(TransformedMatrix),tipvariance));
	gsl_vector *tipresiduals=gsl_vector_calloc(ntax);
	//cout<<"RateTimesVCV cell 0,0 ="<<gsl_matrix_get(RateTimesVCV,0,0)<<endl;
	//cout<<"VCVfinal cell 0,0 = "<<gsl_matrix_get(VCVfinal,0,0)<<endl;
	//cout<<"Size VCVfinal ="<<VCVfinal->size1<<" size of observedtips = "<<observedtips->size<<endl;
	// double ancestralstate=brownie.GetAncestralState(VCVfinal,observedtips);
	//cout<<"Ancestral state = "<<ancestralstate<<endl;
	gsl_vector_memcpy(tipresiduals,brownie.GetTipResiduals(observedtips,ancestralstate));
	double likelihood=(brownie.GetLScore(VCVfinal,tipresiduals,1)); //-lnL actually
	if (rate<0) {
		likelihood=GSL_POSINF;
	}
	//cout<<"likelihood "<<likelihood<<endl;
	gsl_matrix_free(VCV);
	gsl_vector_free(tipvariance);
	gsl_vector_free(observedtips);
	gsl_matrix_free(TransformedMatrix);
	gsl_matrix_free(VCVfinal);
	gsl_vector_free(tipresiduals);
	return likelihood; //-lnL actually
}

double OptimizationFn::GetLikelihoodUnderDelta(const gsl_vector * variables)
{
	double rate=gsl_vector_get(variables,0);
	double ancestralstate=gsl_vector_get(variables,1);
	double delta=gsl_vector_get(variables,2);
	int ntax=Matrix1->size1;
	gsl_matrix *VCV=gsl_matrix_calloc(ntax,ntax);
	gsl_vector *tipvariance=gsl_vector_calloc(ntax);
	gsl_vector *observedtips=gsl_vector_calloc(ntax);
	gsl_matrix_memcpy(VCV,Matrix1);
	gsl_vector_memcpy(tipvariance,Vector2);
	gsl_vector_memcpy(observedtips,Vector1);
	gsl_matrix *TransformedMatrix=brownie.ConvertVCVwithDelta(VCV,delta);
	gsl_matrix *VCVfinal=gsl_matrix_calloc(ntax,ntax);
	gsl_matrix_add_constant (TransformedMatrix, -1.0*gsl_matrix_min (TransformedMatrix)); //replaces DeleteStem
	gsl_matrix_memcpy(VCVfinal,TransformedMatrix);
	for (int r=0;r<ntax;r++) {
        gsl_matrix_set(VCVfinal,r,r,(gsl_matrix_get(TransformedMatrix,r,r)+(gsl_vector_get(tipvariance,r)))); //replaces AddTipVarianceVectorToRateTImesVCV
    }
	
	//gsl_matrix_memcpy(VCVfinal,brownie.AddTipVarianceVectorToRateTimesVCV(brownie.DeleteStem(TransformedMatrix),tipvariance));
	gsl_vector *tipresiduals=gsl_vector_calloc(ntax);
	//cout<<"RateTimesVCV cell 0,0 ="<<gsl_matrix_get(RateTimesVCV,0,0)<<endl;
	//cout<<"VCVfinal cell 0,0 = "<<gsl_matrix_get(VCVfinal,0,0)<<endl;
	//cout<<"Size VCVfinal ="<<VCVfinal->size1<<" size of observedtips = "<<observedtips->size<<endl;
	// double ancestralstate=brownie.GetAncestralState(VCVfinal,observedtips);
	//cout<<"Ancestral state = "<<ancestralstate<<endl;
	gsl_vector_memcpy(tipresiduals,brownie.GetTipResiduals(observedtips,ancestralstate));
	double likelihood=(brownie.GetLScore(VCVfinal,tipresiduals,1)); //-lnL actually
	if (rate<0 || delta<0) {
		likelihood=GSL_POSINF;
	}
	//cout<<"likelihood "<<likelihood<<endl;
	//cout<<likelihood<<" "<<delta<<" "<<rate<<" "<<ancestralstate<<endl;
	gsl_matrix_free(VCV);
	gsl_vector_free(tipvariance);
	gsl_vector_free(observedtips);
	gsl_matrix_free(TransformedMatrix);
	gsl_matrix_free(VCVfinal);
	gsl_vector_free(tipresiduals);
	return likelihood; //-lnL actually
}

double OptimizationFn::GetLikelihoodUnderLambda(const gsl_vector * variables)
{
	double rate=gsl_vector_get(variables,0);
	double ancestralstate=gsl_vector_get(variables,1);
	double lambda=gsl_vector_get(variables,2);
	int ntax=Matrix1->size1;
	gsl_matrix *VCV=gsl_matrix_calloc(ntax,ntax);
	gsl_vector *tipvariance=gsl_vector_calloc(ntax);
	gsl_vector *observedtips=gsl_vector_calloc(ntax);
	gsl_matrix_memcpy(VCV,Matrix1);
	gsl_vector_memcpy(tipvariance,Vector2);
	gsl_vector_memcpy(observedtips,Vector1);
	gsl_matrix *TransformedMatrix=brownie.ConvertVCVwithLambda(VCV,lambda);
	gsl_matrix *VCVfinal=gsl_matrix_calloc(ntax,ntax);
	gsl_matrix_add_constant (TransformedMatrix, -1.0*gsl_matrix_min (TransformedMatrix)); //replaces DeleteStem
		gsl_matrix_memcpy(VCVfinal,TransformedMatrix);
	for (int r=0;r<ntax;r++) {
        gsl_matrix_set(VCVfinal,r,r,(gsl_matrix_get(TransformedMatrix,r,r)+(gsl_vector_get(tipvariance,r)))); //replaces AddTipVarianceVectorToRateTImesVCV
    }
	
	//gsl_matrix_memcpy(VCVfinal,brownie.AddTipVarianceVectorToRateTimesVCV(brownie.DeleteStem(TransformedMatrix),tipvariance));
	gsl_vector *tipresiduals=gsl_vector_calloc(ntax);
	//cout<<"RateTimesVCV cell 0,0 ="<<gsl_matrix_get(RateTimesVCV,0,0)<<endl;
	//cout<<"VCVfinal cell 0,0 = "<<gsl_matrix_get(VCVfinal,0,0)<<endl;
	//cout<<"Size VCVfinal ="<<VCVfinal->size1<<" size of observedtips = "<<observedtips->size<<endl;
	// double ancestralstate=brownie.GetAncestralState(VCVfinal,observedtips);
	//cout<<"Ancestral state = "<<ancestralstate<<endl;
	gsl_vector_memcpy(tipresiduals,brownie.GetTipResiduals(observedtips,ancestralstate));
	double likelihood=(brownie.GetLScore(VCVfinal,tipresiduals,1)); //-lnL actually
	if (rate<0 || lambda<0) {
		likelihood=GSL_POSINF;
	}
	//cout<<"likelihood "<<likelihood<<endl;
	gsl_matrix_free(VCV);
	gsl_vector_free(tipvariance);
	gsl_vector_free(observedtips);
	gsl_matrix_free(TransformedMatrix);
	gsl_matrix_free(VCVfinal);
	gsl_vector_free(tipresiduals);
	return likelihood; //-lnL actually
}

double OptimizationFn::GetLikelihoodWithOptimizedTipVariance(const gsl_vector * variables) // where params points to the MatrixVectorVector struct //gsl_matrix * VCV, gsl_vector *tipvariance, gsl_vector * observedtips)
{
    //cout<<"Now in OptimizationFn::GetLikelihoodWithGivenTipVariance"<<endl;
    //cout<<"Size of variables vector = "<<variables->size<<endl;
    //gsl_matrix *VCV=((struct MatrixVectorVector *)params)->Matrix1;
    //gsl_vector *tipvariance=((struct MatrixVectorVector *)params)->Vector2;
    //gsl_vector *observedtips=((struct MatrixVectorVector *)params)->Vector1;
	double rate=gsl_vector_get(variables,0);
	double tipvar=gsl_vector_get(variables,1);
    //cout<<"Rate is "<<rate<<endl;
	int ntax=Matrix1->size1;
	//cout<<"Ntax = "<<ntax<<endl;
	gsl_matrix *VCV=gsl_matrix_calloc(ntax,ntax);
	gsl_vector *tipvariance=gsl_vector_calloc(ntax);
	gsl_vector *observedtips=gsl_vector_calloc(ntax);
	gsl_matrix_memcpy(VCV,Matrix1);
	gsl_vector_set_all(tipvariance,tipvar);
	gsl_vector_memcpy(observedtips,Vector1);
	// cout<<"Observed state in first taxon = "<<gsl_vector_get(observedtips,0)<<endl;
	// cout<<"Tip variance in first taxon = "<<gsl_vector_get(tipvariance,0)<<endl;
	gsl_matrix *RateTimesVCV=gsl_matrix_calloc(ntax,ntax);
	gsl_matrix_memcpy(RateTimesVCV, VCV);
	gsl_matrix_scale(RateTimesVCV,rate);
	gsl_matrix *VCVfinal=gsl_matrix_calloc(ntax,ntax);
	gsl_matrix_add_constant (RateTimesVCV, -1.0*gsl_matrix_min (RateTimesVCV)); //replaces DeleteStem
		gsl_matrix_memcpy(VCVfinal,RateTimesVCV);
	for (int r=0;r<ntax;r++) {
        gsl_matrix_set(VCVfinal,r,r,(gsl_matrix_get(RateTimesVCV,r,r)+(gsl_vector_get(tipvariance,r)))); //replaces AddTipVarianceVectorToRateTImesVCV
    }
	
	//gsl_matrix_memcpy(VCVfinal,brownie.AddTipVarianceVectorToRateTimesVCV(brownie.DeleteStem(RateTimesVCV),tipvariance));
	gsl_vector *tipresiduals=gsl_vector_calloc(ntax);
    //cout<<"RateTimesVCV cell 0,0 ="<<gsl_matrix_get(RateTimesVCV,0,0)<<endl;
    //cout<<"VCVfinal cell 0,0 = "<<gsl_matrix_get(VCVfinal,0,0)<<endl;
    //cout<<"Size VCVfinal ="<<VCVfinal->size1<<" size of observedtips = "<<observedtips->size<<endl;
	double ancestralstate=brownie.GetAncestralState(VCVfinal,observedtips);
    //cout<<"Ancestral state = "<<ancestralstate<<endl;
	gsl_vector_memcpy(tipresiduals,brownie.GetTipResiduals(observedtips,ancestralstate));
	double likelihood=(brownie.GetLScore(VCVfinal,tipresiduals,1)); //-lnL actually
	if (tipvar<0 || rate<0) {
		likelihood=GSL_POSINF;
	}
	//cout<<"rate = "<<rate<<" ancstate ="<<ancestralstate<<" tipvar = "<<tipvar<<" likelihood = "<<likelihood<<endl;
	gsl_matrix_free(VCV);
	gsl_vector_free(tipvariance);
	gsl_vector_free(observedtips);
	gsl_matrix_free(RateTimesVCV);
	gsl_matrix_free(VCVfinal);
	gsl_vector_free(tipresiduals);
	return likelihood; //-lnL actually
}


double OptimizationFn::GetLikelihoodWithOptimizedTipVariance_gsl(const gsl_vector * variables, void *obj)
{
	double temp;
	temp= ((OptimizationFn*)obj)->GetLikelihoodWithOptimizedTipVariance(variables);
	return temp;
}

gsl_vector * OptimizationFn::OptimizeRateWithOptimizedTipVariance()
{
	//gsl_vector * inputrate=gsl_vector_calloc(1);
	//gsl_vector_set(inputrate,0,1);
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	size_t np = 2;
	size_t iter = 0, i;
	int status;
	double size;
	/* Initial vertex size vector */
	ss = gsl_vector_alloc (np);
	/* Set all step sizes to .01 */ //Note that it was originally 1
	gsl_vector_set_all (ss, .001);
	/* Starting point */
	//cout<<"Now in OPtimizaRateWithGivenTipVariance in OptimizationFn"<<endl;
	x = gsl_vector_alloc (np);
	gsl_vector_set (x,0,1);
	gsl_vector_set (x,1,0.01);
	OptimizationFn *pt;
	pt=(this);
	double (*F)(const gsl_vector *, void *);
	F = &OptimizationFn::GetLikelihoodWithOptimizedTipVariance_gsl;
	/* Initialize method and iterate */
	gsl_multimin_function minex_func;
	minex_func.f=*F;
	minex_func.params=pt;
	minex_func.n = np;
	s = gsl_multimin_fminimizer_alloc (T, np);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
	do
	{
		cout<<"Now on iteration "<<iter<<endl;
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status!=0) { //0 Means it's a success
			printf ("error: %s\n", gsl_strerror (status));
			break;
		}
		size = gsl_multimin_fminimizer_size (s);
		//status = gsl_multimin_test_size (size, 1e-2);
		status = gsl_multimin_test_size (size, brownie.stoppingprecision); //since we want more precision
		if (status == GSL_SUCCESS)
		{
			//printf ("converged to minimum at\n");
		}
		//printf ("%5d ", iter);
		for (i = 0; i < np; i++)
		{
			//printf ("%10.3e ", gsl_vector_get (s->x, i));
		}
		//printf ("f() = %7.3f size = %.3f\n", s->fval, size);
	}
	while (status == GSL_CONTINUE && iter < maxiterations);
	gsl_vector * results=gsl_vector_calloc(np);
	gsl_vector_memcpy(results,s->x);
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);
	//cout<<gsl_vector_get(results,0)<<" "<<gsl_vector_get(results,1)<<endl;
	return results;
	
};


double OptimizationFn::GetLikelihoodUnderACDC_gsl( const gsl_vector * variables, void *obj) 
{
	double temp;
	temp= ((OptimizationFn*)obj)->GetLikelihoodUnderACDC(variables);
	return temp;
}

double OptimizationFn::GetLikelihoodUnderOU1_gsl( const gsl_vector * variables, void *obj) 
{
	double temp;
	temp= ((OptimizationFn*)obj)->GetLikelihoodUnderOU1(variables);
	return temp;
}

double OptimizationFn::GetLikelihoodUnderDelta_gsl( const gsl_vector * variables, void *obj) 
{
	double temp;
	temp= ((OptimizationFn*)obj)->GetLikelihoodUnderDelta(variables);
	return temp;
}

double OptimizationFn::GetLikelihoodUnderLambda_gsl( const gsl_vector * variables, void *obj) 
{
	double temp;
	temp= ((OptimizationFn*)obj)->GetLikelihoodUnderLambda(variables);
	return temp;
}

gsl_vector * OptimizationFn::GeneralOptimization(int ChosenModel)
{
	//Model 1=Optimize rate numerically under simple brownian motion, with given tip variance (all given tip variances can also be zero)
	//Model 2=Optimize rate and tip variance numerically.
	//Model 3=Optimize rate, OU mean, and attraction parameter, with given tip variance. Blomberg et al's d.
	//Model 4=Optimize rate, mean, and ACDC parameter, with given tip variance. Blomberg et al's g.
	double bestlikelihood=GSL_POSINF;
	int hitlimitscount=0;
	size_t np;
	if (ChosenModel==1) {
		np = 1;
	}
	else if (ChosenModel==2) {
		np = 2;
	}
	else if (ChosenModel==3) {
		np = 3;
	}
	else if (ChosenModel==4) {
		np = 3;
	}
	else if (ChosenModel==21) {
		np = 3;
	}
	else if (ChosenModel==22) {
		np = 3;
	}
	gsl_vector * results=gsl_vector_calloc(np);
	double startingancestralstatemean=brownie.GetAncestralState(Matrix1,Vector1);
	double startingratemean=brownie.EstimateRate(Matrix1,brownie.GetTipResiduals(Vector1,startingancestralstatemean));
	double estimates[randomstarts][np];
	double startingvalues[randomstarts][np];
	double likelihoods[randomstarts][1];
	if(detailedoutput==false) {
		brownie.ProgressBar(randomstarts);
	}
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
			gsl_vector_set (x,0,gsl_ran_exponential (r,startingratemean)); //starting rate
			startingvalues[startnum][0]=gsl_vector_get(x,0);
		}
		else if (ChosenModel==2) {
			gsl_vector_set (x,0,gsl_ran_exponential (r,startingratemean)); //starting rate
			startingvalues[startnum][0]=gsl_vector_get(x,0);
			gsl_vector_set (x,1,gsl_ran_exponential (r,startingratemean*gsl_matrix_get(Matrix1,0,0)/10)); // guess a mean 1/10 of the height of the VCV matrix under simple brownian motion
			startingvalues[startnum][1]=gsl_vector_get(x,1);
		}
		else if (ChosenModel==3) {
			gsl_vector_set (x,0,gsl_ran_exponential (r,startingratemean)); //starting rate
			startingvalues[startnum][0]=gsl_vector_get(x,0);
			gsl_vector_set (x,1,gsl_ran_exponential (r,startingancestralstatemean)); //starting ancestral state
			startingvalues[startnum][1]=gsl_vector_get(x,1);
			gsl_vector_set (x,2,gsl_ran_flat (r,0,1)); //starting value of d
			startingvalues[startnum][2]=gsl_vector_get(x,2);
		}
		else if (ChosenModel==4) {
			gsl_vector_set (x,0,gsl_ran_exponential (r,startingratemean)); //starting rate
			startingvalues[startnum][0]=gsl_vector_get(x,0);
			gsl_vector_set (x,1,gsl_ran_exponential (r,startingancestralstatemean)); //starting ancestral state
			startingvalues[startnum][1]=gsl_vector_get(x,1);
			gsl_vector_set (x,2,gsl_ran_exponential (r,1)); //starting value of g
			startingvalues[startnum][2]=gsl_vector_get(x,2);
		}
		else if (ChosenModel==21) {
			gsl_vector_set (x,0,gsl_ran_exponential (r,startingratemean)); //starting rate
			startingvalues[startnum][0]=gsl_vector_get(x,0);
			gsl_vector_set (x,1,gsl_ran_exponential (r,startingancestralstatemean)); //starting ancestral state
			startingvalues[startnum][1]=gsl_vector_get(x,1);
			gsl_vector_set (x,2,gsl_ran_exponential (r,1)); //starting value of delta
			startingvalues[startnum][2]=gsl_vector_get(x,2);
		}
		else if (ChosenModel==22) {
			gsl_vector_set (x,0,gsl_ran_exponential (r,startingratemean)); //starting rate
			startingvalues[startnum][0]=gsl_vector_get(x,0);
			gsl_vector_set (x,1,gsl_ran_exponential (r,startingancestralstatemean)); //starting ancestral state
			startingvalues[startnum][1]=gsl_vector_get(x,1);
			gsl_vector_set (x,2,gsl_ran_exponential (r,1)); //starting value of lambda
			startingvalues[startnum][2]=gsl_vector_get(x,2);
		}
		
		
		OptimizationFn *pt;
		pt=(this);
		double (*F)(const gsl_vector *, void *);
		if (ChosenModel==1) {
			F = &OptimizationFn::GetLikelihoodWithGivenTipVariance_gsl;
		}
		else if (ChosenModel==2) {
			F = &OptimizationFn::GetLikelihoodWithOptimizedTipVariance_gsl;
		}
		else if (ChosenModel==3) {
			F = &OptimizationFn::GetLikelihoodUnderOU1_gsl;
		}
		else if (ChosenModel==4) {
			F = &OptimizationFn::GetLikelihoodUnderACDC_gsl;
		}
		else if (ChosenModel==21) {
			F = &OptimizationFn::GetLikelihoodUnderDelta_gsl;
		}
		else if (ChosenModel==22) {
			F = &OptimizationFn::GetLikelihoodUnderLambda_gsl;
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
		brownie.message="Replicate ";
		brownie.message+=startnum+1;
		if (hitlimits) {
			brownie.message+=" **WARNING**";
		}
		brownie.message+="\n   NM iterations needed = ";
		int iterationsrequired=iter;
		brownie.message+=iterationsrequired;
		if (hitlimits) {
			brownie.message+=" **Max iterations hit; see WARNING below**";
		}
		brownie.message+="\n   LnL = ";
		char outputstring[60];
		sprintf(outputstring,"%60.45f",-1*(s->fval));
		brownie.message+=outputstring;
		// brownie.message+="\n   Rate = ";
		// brownie.message+=gsl_vector_get(s->x,0);
		if (detailedoutput) {
			brownie.PrintMessage();
		}
		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free (s);
		if(detailedoutput==false) {
			brownie.ProgressBar(0);
		}
	}
	if (hitlimitscount>0) {
		brownie.message="\n----------------------------------------------------------------------------\n WARNING: Out of ";
		brownie.message+=randomstarts;
		brownie.message+=" optimization starts, ";
		brownie.message+=hitlimitscount;
		if (hitlimitscount==1) {
			brownie.message+=" was ";
		}
		else {
			brownie.message+=" were ";
		}
		brownie.message+="stopped by hitting\n  the maximum # of iterations. This means that those replicates\n  may not even have hit the local maximum.\n\n  You can increase the maximum number of iterations or decrease the\n  precision with the NumOpt command. You could also consider\n  increasing the number of random starts using that same command.\n\n  If this happened on a small proportion of replicates, though,\n  or if the precision (below) is good enough, don't worry about it.\n----------------------------------------------------------------------------";
		if (detailedoutput) {
			brownie.PrintMessage();
		}
		else if ((randomstarts-hitlimitscount)<10 && (hitlimitscount/randomstarts)>.1) {
			brownie.PrintMessage();
		}
	}
	brownie.message="\n\nRate ML estimate = ";
	brownie.message+=gsl_vector_get(results,0);
	// brownie.message+="\nMean estimate across starts = ";
	// brownie.message+=gsl_stats_mean(rateestimates,1,randomstarts);
	
	if (detailedoutput) {
		brownie.PrintMessage();
	}
	// brownie.message="\nStandard deviation of estimates across ";
	//  brownie.message+=randomstarts;
	// brownie.message+=" starts = ";
	// brownie.message+=gsl_stats_sd(rateestimates,1,randomstarts);
	//  brownie.message+="\n[This is the precision of the rate estimate: numerical optimization does not give an exact value]";
	//  brownie.PrintMessage();
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
};

////////////////////////////////////////////////////////////////////////

gsl_vector * OptimizationFnMultiModel::GeneralOptimization(int ChosenModel)
{
	double bestlikelihood=GSL_POSINF;
	int hitlimitscount=0;
	size_t np;
	int numberofnontrivialmatrices=0;
	//We only want to adjust rates for matrices with some nonzero entries (adjusting rates on the others won't affect likelihood).
	//We assume that matrices are filled in order:	if Matrix3 has non-zero entries, then Matrices0, 1, and 2 should also have
	//non-zero entries.
	if (ChosenModel==5 || ChosenModel==12) {
		if (gsl_matrix_max(Matrix9) > 0) {
			numberofnontrivialmatrices=10;
		}
		else if (gsl_matrix_max(Matrix8) > 0) {
			numberofnontrivialmatrices=9;
		}
		else if (gsl_matrix_max(Matrix7) > 0) {
			numberofnontrivialmatrices=8;
		}
		else if (gsl_matrix_max(Matrix6) > 0) {
			numberofnontrivialmatrices=7;
		}
		else if (gsl_matrix_max(Matrix5) > 0) {
			numberofnontrivialmatrices=6;
		}
		else if (gsl_matrix_max(Matrix4) > 0) {
			numberofnontrivialmatrices=5;
		}
		else if (gsl_matrix_max(Matrix3) > 0) {
			numberofnontrivialmatrices=4;
		}
		else if (gsl_matrix_max(Matrix2) > 0) {
			numberofnontrivialmatrices=3;
		}
		else if (gsl_matrix_max(Matrix1) > 0) {
			numberofnontrivialmatrices=2;
		}
		else if (gsl_matrix_max(Matrix0) > 0) {
			numberofnontrivialmatrices=1;
		}
		else {
			brownie.message="All the VCV matrices had entries of zero. Are you sure you loaded a tree with branch lengths?";
			brownie.PrintMessage();
		}
		if (ChosenModel==5) {
			np = 1+numberofnontrivialmatrices;
		}
		if (ChosenModel==12) {
			np = 2+numberofnontrivialmatrices; //BM param, OU param, rootstate param, then all the optima. Matrix 0 has the tree, Matrices>=1 have start stop times
		}
	}
	if (ChosenModel==6) {
		np=3;
	}
	size_t npouterloop=np;
	if (ChosenModel==12) {
		npouterloop=2;
	}
	
	const gsl_rng_type * T_rng_type;
	gsl_rng * r_rng;
	gsl_rng_env_setup();
	T_rng_type = gsl_rng_default;
	r_rng = gsl_rng_alloc (T_rng_type);
	
	gsl_vector * results=gsl_vector_calloc(npouterloop);
	gsl_vector * resultsallparam=gsl_vector_calloc(np);
	fixedparams=gsl_vector_calloc(np);
	//cout<<"fixedparams->size = "<<fixedparams->size<<endl;
	int ntax=Matrix0->size1;
	gsl_matrix * CombinedVCV = gsl_matrix_calloc(ntax,ntax);
	if (ChosenModel==5) {
		gsl_matrix_add(CombinedVCV,Matrix0);
		gsl_matrix_add(CombinedVCV,Matrix1);
		gsl_matrix_add(CombinedVCV,Matrix2);
		gsl_matrix_add(CombinedVCV,Matrix3);
		gsl_matrix_add(CombinedVCV,Matrix4);
		gsl_matrix_add(CombinedVCV,Matrix5);
		gsl_matrix_add(CombinedVCV,Matrix6);
		gsl_matrix_add(CombinedVCV,Matrix7);
		gsl_matrix_add(CombinedVCV,Matrix8);
		gsl_matrix_add(CombinedVCV,Matrix9);
	}
	else if (ChosenModel==6) {
		gsl_matrix_add(CombinedVCV,Matrix0);
		gsl_matrix_add(CombinedVCV,Matrix1);	
	}
	else if (ChosenModel==12) {
		gsl_matrix_add(CombinedVCV,Matrix0);
	}
	double startingancestralstatemean=brownie.GetAncestralState(CombinedVCV,Vector1);
	double startingratemean=brownie.EstimateRate(CombinedVCV,brownie.GetTipResiduals(Vector1,startingancestralstatemean));
	double estimates[randomstarts][np];
	double startingvalues[randomstarts][npouterloop];
	double likelihoods[randomstarts][1];
	if(detailedoutput==false) {
		brownie.ProgressBar(randomstarts);
	}
	for (int startnum=0;startnum<randomstarts;startnum++) {
		const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
		gsl_multimin_fminimizer *s = NULL;
		gsl_vector *ss, *x;
		size_t iter = 0, i;
		int status;
		double size;
		bool hitlimits=false;
		/* Initial vertex size vector */
		ss = gsl_vector_alloc (npouterloop);
		gsl_vector_set_all (ss, stepsize);
		
		/* Starting point */
		x = gsl_vector_calloc (npouterloop);
		if (ChosenModel==5) {
			//first value is ancestral state, others are rates
			gsl_vector_set (x,0,startingancestralstatemean+gsl_ran_gaussian(r,startingratemean));
			startingvalues[startnum][0]=gsl_vector_get(x,0);
			for (int cell=1; cell<np;cell++) {
				gsl_vector_set (x,cell,gsl_ran_exponential (r,startingratemean));	
				startingvalues[startnum][cell]=gsl_vector_get(x,cell);
			}
		}
		if (ChosenModel==6) {
			//first value is ancestral state, others are rates
			gsl_vector_set (x,0,startingancestralstatemean+gsl_ran_gaussian(r,startingratemean));
			startingvalues[startnum][0]=gsl_vector_get(x,0);
			gsl_vector_set (x,1,gsl_ran_exponential (r,startingratemean));	
			startingvalues[startnum][1]=gsl_vector_get(x,1);
			gsl_vector_set (x,2,gsl_ran_exponential (r,startingratemean));	
			startingvalues[startnum][2]=gsl_vector_get(x,2);
		}
		if (ChosenModel==12) { //OUSM
			gsl_vector_set (x,0,gsl_ran_exponential (r,startingratemean)); //BM rate param
			startingvalues[startnum][0]=gsl_vector_get(x,0);
			gsl_vector_set (x,1,gsl_ran_exponential (r,1)); //OU attraction param. Might want to have a way of tuning this distribution	
			while ((gsl_vector_get(x,1)<minATTRACTIONPARAM) || (gsl_vector_get(x,1)>maxATTRACTIONPARAM)) {
				gsl_vector_set (x,1,gsl_ran_exponential (r,1)); //Don't want to bother with an illegal starting point
			}
			startingvalues[startnum][0]=gsl_vector_get(x,0);
			startingvalues[startnum][1]=gsl_vector_get(x,1);

			//for (int cell=2; cell<np;cell++) {
			//	gsl_vector_set (x,cell,startingancestralstatemean+gsl_ran_gaussian(r,startingratemean));	
			//	startingvalues[startnum][cell]=gsl_vector_get(x,cell);
			//}

		}
			
			OptimizationFnMultiModel *pt;
			pt=(this);
			double (*F)(const gsl_vector *, void *);
			if (ChosenModel==5) {
				F = &OptimizationFnMultiModel::GetLikelihoodWithGivenTipVarianceOneRatePerState_gsl;
			}
			else if (ChosenModel==6) {
				F = &OptimizationFnMultiModel::GetLikelihoodWithGivenTipVarianceDiffRateOnChange_gsl;
			}
			else if (ChosenModel==12) {
				//F = &OptimizationFnMultiModel::GetLikelihoodOUSM_gsl;
				//F = &OptimizationFnMultiModel::GetLikelihoodOUSM_OnlyVariablesAttractionRate_gsl;
				F = &OptimizationFnMultiModel::GetLikelihoodOUSM_AnalyticMeans_gsl;
			}
			gsl_multimin_function minex_func;
			minex_func.f=*F;
			minex_func.params=pt;
			minex_func.n = npouterloop;
			s = gsl_multimin_fminimizer_alloc (T, npouterloop);
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
				if(ChosenModel==12) {
					for (int i=0;i<npouterloop;i++) {
						gsl_vector_set(resultsallparam,i,gsl_vector_get(results,i));
					}
					for (int i=npouterloop;i<np;i++) {
						gsl_vector_set(resultsallparam,i,gsl_vector_get(fixedparams,i));
					}
				}
				else {
					gsl_vector_memcpy(resultsallparam,results);
				}
				bestlikelihood=s->fval;
			}
			for (int parameternumber=0; parameternumber<npouterloop; parameternumber++) {
				estimates[startnum][parameternumber]=gsl_vector_get(s->x,parameternumber);
			}
			for (int parameternumber=npouterloop; parameternumber<np; parameternumber++) {
				estimates[startnum][parameternumber]=gsl_vector_get(fixedparams,parameternumber);
			}
			if (iter==maxiterations) {
				hitlimits=true;
				hitlimitscount++;
			}
			brownie.message="Replicate ";
			brownie.message+=startnum+1;
			if (hitlimits) {
				brownie.message+=" **WARNING**";
			}
			brownie.message+="\n   NM iterations needed = ";
			int iterationsrequired=iter;
			brownie.message+=iterationsrequired;
			if (hitlimits) {
				brownie.message+=" **Max iterations hit; see WARNING below**";
			}
			brownie.message+="\n   LnL = ";
			char outputstring[60];
			sprintf(outputstring,"%60.45f",-1*(s->fval));
			brownie.message+=outputstring;
			if (ChosenModel==5) {
			brownie.message+="\n   RootVal = ";
			brownie.message+=gsl_vector_get(s->x,0);
				for (int cell=1;cell<np;cell++) {
					brownie.message+="\n   Rate in state ";
					brownie.message+=cell-1;
					brownie.message+=" = ";
					brownie.message+=gsl_vector_get(s->x,cell);
				}
			}
			else if (ChosenModel==6) {
			brownie.message+="\n   RootVal = ";
			brownie.message+=gsl_vector_get(s->x,0);
				brownie.message+="\n   Rate on branches with zero changes = ";
				brownie.message+=gsl_vector_get(s->x,1);
				brownie.message+="\n   Rate on branches with changes = ";
				brownie.message+=gsl_vector_get(s->x,2);
			}
			else if (ChosenModel==12) {
				brownie.message+="\n  BM Rate = ";
				brownie.message+=gsl_vector_get(s->x,0);
				brownie.message+="\n  OU attraction = ";
				brownie.message+=gsl_vector_get(s->x,1);
				brownie.message+="\n  Anc state = ";
				brownie.message+=gsl_vector_get(fixedparams,2);
				for (int cell=3;cell<np;cell++) {
					brownie.message+="\n  Mean value in state ";
					brownie.message+=cell-3;
					brownie.message+=" = ";
					brownie.message+=gsl_vector_get(fixedparams,cell);
				}
			}
			if (detailedoutput) {
				brownie.PrintMessage();
			}
			gsl_vector_free(x);
			gsl_vector_free(ss);
			gsl_multimin_fminimizer_free (s);
			if(detailedoutput==false) {
				brownie.ProgressBar(0);
			}
	}
	if (hitlimitscount>0) {
		brownie.message="\n----------------------------------------------------------------------------\n WARNING: Out of ";
		brownie.message+=randomstarts;
		brownie.message+=" optimization starts, ";
		brownie.message+=hitlimitscount;
		if (hitlimitscount==1) {
			brownie.message+=" was ";
		}
		else {
			brownie.message+=" were ";
		}
		brownie.message+="stopped by hitting\n  the maximum # of iterations. This means that those replicates\n  may not even have hit the local maximum.\n\n  You can increase the maximum number of iterations or decrease the\n  precision with the NumOpt command. You could also consider\n  increasing the number of random starts using that same command.\n\n  If this happened on a small proportion of replicates, though,\n  or if the precision (below) is good enough, don't worry about it.\n----------------------------------------------------------------------------";
		if (detailedoutput) {
			brownie.PrintMessage();
		}
		else if ((randomstarts-hitlimitscount)<10 && (hitlimitscount/randomstarts)>.1) {
			brownie.PrintMessage();
		}
	}
/*	brownie.message="\n\nRoot value estimate = ";
	brownie.message+=gsl_vector_get(results,0);
	for (int cell=1;cell<np;cell++) {
		brownie.message+="\n   Rate in state ";
		brownie.message+=cell-1;
		brownie.message+=" estimate = ";
		brownie.message+=gsl_vector_get(results,cell);
	}
*/	
	if (detailedoutput) {
		brownie.PrintMessage();
	}
// brownie.message="\nStandard deviation of estimates across ";
//  brownie.message+=randomstarts;
// brownie.message+=" starts = ";
// brownie.message+=gsl_stats_sd(rateestimates,1,randomstarts);
//  brownie.message+="\n[This is the precision of the rate estimate: numerical optimization does not give an exact value]";
//  brownie.PrintMessage();
	gsl_vector * finalvector;
	if (ChosenModel==5) {
		finalvector=gsl_vector_calloc(24); //change this if you change the max number of models
									   //finalvector=[lnL np ancstate rate1 rate2 rate3 ... rate10 sd(ancstate) sd(rate1) sd(rate2) ... sd(rate10)]
		gsl_vector_set(finalvector,0,bestlikelihood);
		gsl_vector_set(finalvector,1,np); //stores the number of optimized params
		for (int position=2; position<(np+2); position++) {
			gsl_vector_set(finalvector,position,gsl_vector_get(results,position-2));
			double paramestimate[randomstarts];
			for (int startnumber=0;startnumber<randomstarts;startnumber++) {
				paramestimate[startnumber]=estimates[startnumber][position-2];
			}
			gsl_vector_set(finalvector,position+11,gsl_stats_sd(paramestimate,1,randomstarts));
		}
		
	}
	else if (ChosenModel==6) {
		finalvector=gsl_vector_calloc(7);
		gsl_vector_set(finalvector,0,bestlikelihood);
		gsl_vector_set(finalvector,1,gsl_vector_get(results,0));
		gsl_vector_set(finalvector,2,gsl_vector_get(results,1));
		gsl_vector_set(finalvector,3,gsl_vector_get(results,2));
		for (int position=0; position<3; position++) {
			double paramestimate[randomstarts];
			for (int startnumber=0;startnumber<randomstarts;startnumber++) {
				paramestimate[startnumber]=estimates[startnumber][position];
			}
			gsl_vector_set(finalvector,position+4,gsl_stats_sd(paramestimate,1,randomstarts));
		}
	}
	else if (ChosenModel==12) {
		finalvector=gsl_vector_calloc(28); //change this if you change the max number of models
		gsl_vector_set(finalvector,0,bestlikelihood);
		gsl_vector_set(finalvector,1,np); //stores the number of optimized params
		for (int position=2; position<(np+2); position++) {
			gsl_vector_set(finalvector,position,gsl_vector_get(resultsallparam,position-2));
			double paramestimate[randomstarts];
			for (int startnumber=0;startnumber<randomstarts;startnumber++) {
				paramestimate[startnumber]=estimates[startnumber][position-2];
			}
			gsl_vector_set(finalvector,position+13,gsl_stats_sd(paramestimate,1,randomstarts));
		}
	}
//cout<<"Final vector: "<<endl;
//for (int finalvectorposition=0;finalvectorposition<finalvector->size;finalvectorposition++) {
//	cout<<finalvectorposition<<": "<<gsl_vector_get(finalvector,finalvectorposition)<<endl;
//		}
	return finalvector;
}




double OptimizationFnMultiModel::GetLikelihoodWithGivenTipVarianceOneRatePerState(const gsl_vector * variables) {
	double ancstate=gsl_vector_get(variables,0);
    int ntax=Matrix0->size1;
	gsl_matrix *VCVtotal=gsl_matrix_calloc(ntax,ntax);
	gsl_vector *tipvariance=gsl_vector_calloc(ntax);
	gsl_vector *observedtips=gsl_vector_calloc(ntax);
	gsl_vector_memcpy(tipvariance,Vector2);
	gsl_vector_memcpy(observedtips,Vector1);
	int numberofmodels=-1+(variables->size);
	
	gsl_matrix *RateTimesVCV=gsl_matrix_calloc(ntax,ntax);
	gsl_matrix_memcpy(RateTimesVCV, Matrix0);
	gsl_matrix_scale(RateTimesVCV,gsl_vector_get(variables,1));
	gsl_matrix_add(VCVtotal,RateTimesVCV);
	//here's where an eval function or a 3-d matrix structure would come in handy
	if (numberofmodels>1) {
		gsl_matrix_memcpy(RateTimesVCV, Matrix1);
		gsl_matrix_scale(RateTimesVCV,gsl_vector_get(variables,2));
		gsl_matrix_add(VCVtotal,RateTimesVCV);
	}
	if (numberofmodels>2) {
		gsl_matrix_memcpy(RateTimesVCV, Matrix2);
		gsl_matrix_scale(RateTimesVCV,gsl_vector_get(variables,3));
		gsl_matrix_add(VCVtotal,RateTimesVCV);
	}	
	if (numberofmodels>3) {
		gsl_matrix_memcpy(RateTimesVCV, Matrix3);
		gsl_matrix_scale(RateTimesVCV,gsl_vector_get(variables,4));
		gsl_matrix_add(VCVtotal,RateTimesVCV);
	}		
	if (numberofmodels>4) {
		gsl_matrix_memcpy(RateTimesVCV, Matrix4);
		gsl_matrix_scale(RateTimesVCV,gsl_vector_get(variables,5));
		gsl_matrix_add(VCVtotal,RateTimesVCV);
	}	
	if (numberofmodels>5) {
		gsl_matrix_memcpy(RateTimesVCV, Matrix5);
		gsl_matrix_scale(RateTimesVCV,gsl_vector_get(variables,6));
		gsl_matrix_add(VCVtotal,RateTimesVCV);
	}	
	if (numberofmodels>6) {
		gsl_matrix_memcpy(RateTimesVCV, Matrix6);
		gsl_matrix_scale(RateTimesVCV,gsl_vector_get(variables,7));
		gsl_matrix_add(VCVtotal,RateTimesVCV);
	}	
	if (numberofmodels>7) {
		gsl_matrix_memcpy(RateTimesVCV, Matrix7);
		gsl_matrix_scale(RateTimesVCV,gsl_vector_get(variables,8));
		gsl_matrix_add(VCVtotal,RateTimesVCV);
	}	
	if (numberofmodels>8) {
		gsl_matrix_memcpy(RateTimesVCV, Matrix8);
		gsl_matrix_scale(RateTimesVCV,gsl_vector_get(variables,9));
		gsl_matrix_add(VCVtotal,RateTimesVCV);
	}	
	if (numberofmodels>9) {
		gsl_matrix_memcpy(RateTimesVCV, Matrix9);
		gsl_matrix_scale(RateTimesVCV,gsl_vector_get(variables,10));
		gsl_matrix_add(VCVtotal,RateTimesVCV);
	}	
	
	gsl_matrix *VCVfinal=gsl_matrix_calloc(ntax,ntax);
	gsl_matrix_add_constant (VCVtotal, -1.0*gsl_matrix_min (VCVtotal)); //replaces DeleteStem
	gsl_matrix_memcpy(VCVfinal,VCVtotal);
	for (int r=0;r<ntax;r++) {
        gsl_matrix_set(VCVfinal,r,r,(gsl_matrix_get(VCVtotal,r,r)+(gsl_vector_get(tipvariance,r)))); //replaces AddTipVarianceVectorToRateTImesVCV
    }
	
	//gsl_matrix_memcpy(VCVfinal,brownie.AddTipVarianceVectorToRateTimesVCV(brownie.DeleteStem(VCVtotal),tipvariance));
	gsl_vector *tipresiduals=gsl_vector_calloc(ntax);
	gsl_vector_memcpy(tipresiduals,brownie.GetTipResiduals(observedtips,ancstate));
	double likelihood=(brownie.GetLScore(VCVfinal,tipresiduals,1)); //-lnL actually
	for (int cell=1;cell<=numberofmodels;cell++) {
		if(gsl_vector_get(variables,cell)<0) {
			likelihood=GSL_POSINF;
		}
	}
	if (gsl_vector_min(tipvariance)<0) {
		likelihood=GSL_POSINF;
	}
	gsl_vector_free(observedtips);
	gsl_vector_free(tipvariance);
	gsl_matrix_free(VCVtotal);
	gsl_matrix_free(RateTimesVCV);
	gsl_matrix_free(VCVfinal);
	gsl_vector_free(tipresiduals);
	return likelihood; //-lnL actually
}

double OptimizationFnMultiModel::GetLikelihoodWithGivenTipVarianceOneRatePerState_gsl( const gsl_vector * variables, void *obj) {
	double temp;
	temp= ((OptimizationFnMultiModel*)obj)->GetLikelihoodWithGivenTipVarianceOneRatePerState(variables);
	return temp;
	
}

gsl_vector * OptimizationFnMultiModel::OptimizeRateWithGivenTipVarianceOneRatePerState() {
	double bestlikelihood=GSL_POSINF;
	int hitlimitscount=0;
	//We only want to adjust rates for matrices with some nonzero entries (adjusting rates on the others won't affect likelihood).
	//We assume that matrices are filled in order:	if Matrix3 has non-zero entries, then Matrices0, 1, and 2 should also have
	//non-zero entries.
	int numberofnontrivialmatrices=0;
	if (gsl_matrix_max(Matrix9) > 0) {
		numberofnontrivialmatrices=10;
	}
	else if (gsl_matrix_max(Matrix8) > 0) {
		numberofnontrivialmatrices=9;
	}
	else if (gsl_matrix_max(Matrix7) > 0) {
		numberofnontrivialmatrices=8;
	}
	else if (gsl_matrix_max(Matrix6) > 0) {
		numberofnontrivialmatrices=7;
	}
	else if (gsl_matrix_max(Matrix5) > 0) {
		numberofnontrivialmatrices=6;
	}
	else if (gsl_matrix_max(Matrix4) > 0) {
		numberofnontrivialmatrices=5;
	}
	else if (gsl_matrix_max(Matrix3) > 0) {
		numberofnontrivialmatrices=4;
	}
	else if (gsl_matrix_max(Matrix2) > 0) {
		numberofnontrivialmatrices=3;
	}
	else if (gsl_matrix_max(Matrix1) > 0) {
		numberofnontrivialmatrices=2;
	}
	else if (gsl_matrix_max(Matrix0) > 0) {
		numberofnontrivialmatrices=1;
	}
	else {
		brownie.message="All the VCV matrices had entries of zero. Are you sure you loaded a tree with branch lengths?";
		brownie.PrintMessage();
	}
	size_t np = 1+numberofnontrivialmatrices;
	gsl_vector * results=gsl_vector_calloc(np);
	const gsl_rng_type * T_rng_type;
	gsl_rng * r_rng;
	gsl_rng_env_setup();
	T_rng_type = gsl_rng_default;
	r_rng = gsl_rng_alloc (T_rng_type);
	int ntax=Matrix0->size1;
	gsl_matrix * CombinedVCV = gsl_matrix_calloc(ntax,ntax);
	gsl_matrix_add(CombinedVCV,Matrix0);
	gsl_matrix_add(CombinedVCV,Matrix1);
	gsl_matrix_add(CombinedVCV,Matrix2);
	gsl_matrix_add(CombinedVCV,Matrix3);
	gsl_matrix_add(CombinedVCV,Matrix4);
	gsl_matrix_add(CombinedVCV,Matrix5);
	gsl_matrix_add(CombinedVCV,Matrix6);
	gsl_matrix_add(CombinedVCV,Matrix7);
	gsl_matrix_add(CombinedVCV,Matrix8);
	gsl_matrix_add(CombinedVCV,Matrix9);
	
	double ancstatestart=brownie.GetAncestralState(CombinedVCV,Vector1);
	double ratestart=brownie.EstimateRate(CombinedVCV,brownie.GetTipResiduals(Vector1,ancstatestart));
	gsl_matrix * estimates=gsl_matrix_calloc(randomstarts,np);
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
		x = gsl_vector_alloc (np); 
		//first value is ancestral state, others are rates
		gsl_vector_set (x,0,ancstatestart+gsl_ran_gaussian(r,ratestart));
		for (int cell=1; cell<np;cell++) {
			gsl_vector_set (x,cell,gsl_ran_exponential (r,ratestart));			
		}
		OptimizationFnMultiModel *pt;
		pt=(this);
		double (*F)(const gsl_vector *, void *);
		F = &OptimizationFnMultiModel::GetLikelihoodWithGivenTipVarianceOneRatePerState_gsl;
		/* Initialize method and iterate */
		gsl_multimin_function minex_func;
		minex_func.f=*F;
		minex_func.params=pt;
		minex_func.n = np;
		s = gsl_multimin_fminimizer_alloc (T, np);
		gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
		do
		{
			// cout<<"Now on iteration "<<iter<<endl;
			iter++;
			status = gsl_multimin_fminimizer_iterate(s);
			if (status!=0) { //0 Means it's a success in c++
				printf ("error: %s\n", gsl_strerror (status));
				break;
			}
			size = gsl_multimin_fminimizer_size (s);
			status = gsl_multimin_test_size (size, stoppingprecision); //since we want more precision
			if (status == GSL_SUCCESS)
			{
				//printf ("converged to minimum at\n");
			}
			//printf ("%5d ", iter);
			//	for (i = 0; i < np; i++)
			//	{
			//printf ("%10.3e ", gsl_vector_get (s->x, i));
			//	}
			//printf ("f() = %7.3f size = %.3f\n", s->fval, size);
		}
		while (status == GSL_CONTINUE && iter < maxiterations);
		if (s->fval<bestlikelihood) {
			gsl_vector_memcpy(results,s->x);
			bestlikelihood=s->fval;
		}
		for (int cell=0; cell<np; cell++) {
			gsl_matrix_set(estimates,startnum,cell,gsl_vector_get(s->x,cell));
		}
		if (iter==maxiterations) {
			hitlimits=true;
			hitlimitscount++;
		}
		brownie.message="Replicate ";
		brownie.message+=startnum+1;
		if (hitlimits) {
			brownie.message+=" **WARNING**";
		}
		//brownie.message+="\n   Starting value = ";
		//brownie.message+=startingvalue;
		brownie.message+="\n   NM iterations needed = ";
		int iterationsrequired=iter;
		brownie.message+=iterationsrequired;
		if (hitlimits) {
			brownie.message+=" **Max iterations hit; see WARNING below**";
		}
		brownie.message+="\n   LnL = ";
		//brownie.message+=s->fval;
		char outputstring[60];
		sprintf(outputstring,"%60.45f",-1*(s->fval));
		brownie.message+=outputstring;
		brownie.message+="\n   RootVal = ";
		brownie.message+=gsl_vector_get(s->x,0);
		for (int cell=1;cell<np;cell++) {
			brownie.message+="\n   Rate in state ";
			brownie.message+=cell-1;
			brownie.message+=" = ";
			brownie.message+=gsl_vector_get(s->x,cell);
		}
		if (detailedoutput) {
			brownie.PrintMessage();
		}
		//cout<<"Rep "<<startnum+1<<" Iter "<<iter<<" LnL "<<s->fval<<" Rate "<<gsl_vector_get(s->x,0)<<endl;
		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free (s);
	}
	if (hitlimitscount>0) {
		brownie.message="\n----------------------------------------------------------------------------\n WARNING: Out of ";
		brownie.message+=randomstarts;
		brownie.message+=" optimization starts, ";
		brownie.message+=hitlimitscount;
		if (hitlimitscount==1) {
			brownie.message+=" was ";
		}
		else {
			brownie.message+=" were ";
		}
		brownie.message+="stopped by hitting\n  the maximum # of iterations. This means that those replicates\n  may not even have hit the local maximum.\n\n  You can increase the maximum number of iterations or decrease the\n  precision with the NumOpt command. You could also consider\n  increasing the number of random starts using that same command.\n\n  If this happened on a small proportion of replicates, though,\n  or if the precision (below) is good enough, don't worry about it.\n----------------------------------------------------------------------------";
		if (detailedoutput) {
			brownie.PrintMessage();
		}
		else if ((randomstarts-hitlimitscount)<10 && (hitlimitscount/randomstarts)>.1) {
			brownie.PrintMessage();
		}
	}
	brownie.message="\n\nRoot value estimate = ";
	brownie.message+=gsl_vector_get(results,0);
	for (int cell=1;cell<np;cell++) {
		brownie.message+="\n   Rate in state ";
		brownie.message+=cell-1;
		brownie.message+=" estimate = ";
		brownie.message+=gsl_vector_get(results,cell);
	}
	
	//	brownie.message+="\nMean estimate across starts = ";
	//	brownie.message+=gsl_stats_mean(rateestimates,1,randomstarts);
	
	if (detailedoutput) {
		brownie.PrintMessage();
	}
	//	brownie.message="\nStandard deviation of estimates across ";
	//	brownie.message+=randomstarts;
	//	brownie.message+=" starts = ";
	//	brownie.message+=gsl_stats_sd(rateestimates,1,randomstarts);
	//	brownie.message+="\n[This is the precision of the rate estimate: numerical optimization does not give an exact value]";
	//	brownie.PrintMessage();
	return results;
	
}

double OptimizationFnMultiModel::GetLikelihoodWithGivenTipVarianceDiffRateOnChange(const gsl_vector * variables) {
	double ancstate=gsl_vector_get(variables,0);
    int ntax=Matrix0->size1;
	gsl_matrix *VCVtotal=gsl_matrix_calloc(ntax,ntax);
	gsl_vector *tipvariance=gsl_vector_calloc(ntax);
	gsl_vector *observedtips=gsl_vector_calloc(ntax);
	gsl_vector_memcpy(tipvariance,Vector2);
	gsl_vector_memcpy(observedtips,Vector1);
	int numberofmodels=-1+(variables->size);
	gsl_matrix *RateTimesVCV=gsl_matrix_calloc(ntax,ntax);
	gsl_matrix_memcpy(RateTimesVCV, Matrix0);
	gsl_matrix_scale(RateTimesVCV,gsl_vector_get(variables,1));
	gsl_matrix_add(VCVtotal,RateTimesVCV);
	gsl_matrix_memcpy(RateTimesVCV, Matrix1);
	gsl_matrix_scale(RateTimesVCV,gsl_vector_get(variables,2));
	gsl_matrix_add(VCVtotal,RateTimesVCV);
	gsl_matrix *VCVfinal=gsl_matrix_calloc(ntax,ntax);
	gsl_matrix_add_constant (VCVtotal, -1.0*gsl_matrix_min (VCVtotal)); //replaces DeleteStem
	gsl_matrix_memcpy(VCVfinal,VCVtotal);
	for (int r=0;r<ntax;r++) {
        gsl_matrix_set(VCVfinal,r,r,(gsl_matrix_get(VCVtotal,r,r)+(gsl_vector_get(tipvariance,r)))); //replaces AddTipVarianceVectorToRateTImesVCV
    }
	
	//VCVfinal=brownie.AddTipVarianceVectorToRateTimesVCV(brownie.DeleteStem(VCVtotal),tipvariance);
	//gsl_matrix_memcpy(VCVfinal,brownie.AddTipVarianceVectorToRateTimesVCV(brownie.DeleteStem(VCVtotal),tipvariance));
	gsl_vector *tipresiduals=gsl_vector_calloc(ntax);
	gsl_vector_memcpy(tipresiduals,brownie.GetTipResiduals(observedtips,ancstate));
	double likelihood=(brownie.GetLScore(VCVfinal,tipresiduals,1)); //-lnL actually
	if(gsl_vector_get(variables,1)<0 || gsl_vector_get(variables,2)<0) {
		likelihood=GSL_POSINF;
	}
	if (gsl_vector_min(tipvariance)<0) {
		likelihood=GSL_POSINF;
	}
	gsl_matrix_free(VCVtotal);
	gsl_vector_free(tipvariance);
	gsl_vector_free(observedtips);
	gsl_matrix_free(RateTimesVCV);
	gsl_matrix_free(VCVfinal);
	gsl_vector_free(tipresiduals);
	return likelihood; //-lnL actually
}

double OptimizationFnMultiModel::GetLikelihoodWithGivenTipVarianceDiffRateOnChange_gsl( const gsl_vector * variables, void *obj) {
	double temp;
	temp= ((OptimizationFnMultiModel*)obj)->GetLikelihoodWithGivenTipVarianceDiffRateOnChange(variables);
	return temp;
}

double OptimizationFnMultiModel::GetLikelihoodOUSM_gsl( const gsl_vector * variables, void *obj) {
	double temp;
	temp= ((OptimizationFnMultiModel*)obj)->GetLikelihoodOUSM(variables);
	return temp;
}

double OptimizationFnMultiModel::GetLikelihoodOUSM(const gsl_vector * variables) {
	//For OU models, Matrix1 consists of time between root and MRCA of two taxa (so, a typical Brownie VCV).
	//To convert this to an OU VCV (as in Butler & King, equation A5), have to use this matrix to get T-s (entry in the given matrix)
	//and time s (which is total root-tip length minus Matrix1(i,j))
	//Each subsequent matrix will have times for a different model (ie, OU mean). Each may not be square: it will have one 
	//row per taxon, and odd columns will have the start time for an interval spent with that model and even columns will have the end
	//time for the corresponding interval. These matrices can then be used to pull out times to transform using the OU model,
	//as in Butler & King A7.
	//Note that unlike Butler & King 2004 and Hansen 1997, Brownie allows multiple models per branch 
	
	double rate=gsl_vector_get(variables,0);
	double attraction=gsl_vector_get(variables,1);
	int numberofmeans=-3+(variables->size);
	double rootmean=gsl_vector_get(variables,2);
	int ntax=Matrix0->size1;
	gsl_matrix *VCVtotal=gsl_matrix_calloc(ntax,ntax);
	gsl_vector *tipvariance=gsl_vector_calloc(ntax);
	gsl_vector *observedtips=gsl_vector_calloc(ntax);
	gsl_vector *expectedtips=gsl_vector_calloc(ntax);
	gsl_matrix *W_BK_A7=gsl_matrix_calloc(ntax,numberofmeans+1); //W matrix based on equation A7 of Butler and King
		gsl_vector_memcpy(observedtips,Vector1);
		gsl_vector_memcpy(tipvariance,Vector2);
		double likelihood;
		if (attraction<0) {
			likelihood=GSL_POSINF;
		}
		else if (rate<0) {
			likelihood=GSL_POSINF;
		}
		else if (gsl_vector_min(tipvariance)<0) {
			likelihood=GSL_POSINF;
		}
		else if (attraction>maxATTRACTIONPARAM) {
			likelihood=GSL_POSINF;
			if(detailedoutput) {
				brownie.message="Warning: OU attraction parameter of ";
				brownie.message+=attraction;
				brownie.message+=" is greater than the maximum allowed, ";
				brownie.message+=maxATTRACTIONPARAM;
				brownie.PrintMessage();
			}
		}
		else if (attraction<minATTRACTIONPARAM) {
			likelihood=GSL_POSINF;
			if (detailedoutput) {
				brownie.message="Warning: OU attraction parameter of ";
				brownie.message+=attraction;
				brownie.message+=" is less than the minimum allowed, ";
				brownie.message+=minATTRACTIONPARAM;
				brownie.PrintMessage();
			}
		}
		else {
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
			//cout<<"attraction="<<attraction<<" roottotiptime="<<roottotiptime;
//			gsl_error_handler_t * previoushandler=//gsl_set_error_handler_off();
			double exptonegalphaT=0;
/*			int status=brownie.browniesafe_gsl_sf_exp_e10_e(-1.0*attraction*roottotiptime,exptonegalphaTresult);
			if (status) {
				if (status == GSL_ERANGE) {
					cout<<"Underflow (or overflow) -- attraction parameter too big"<<endl;
				}
			}*/
			exptonegalphaT=brownie.browniesafe_gsl_sf_exp(-1.0*attraction*roottotiptime);
			//cout<<" exptonegalphaT="<<exptonegalphaT<<endl;
//			gsl_set_error_handler (previoushandler);
			for (int rowtaxon=0;rowtaxon<ntax;rowtaxon++) {
				for (int coltaxon=0;coltaxon<ntax;coltaxon++) {
					gsl_matrix_set(ScaledVCV,rowtaxon,coltaxon,(0.5*rate/attraction)*(brownie.browniesafe_gsl_sf_exp(-2.0*attraction*(roottotiptime-gsl_matrix_get(BranchingTimes,rowtaxon,coltaxon))))*(1.0-(brownie.browniesafe_gsl_sf_exp(-2.0*attraction*(gsl_matrix_get(BranchingTimes,rowtaxon,coltaxon))))));
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,0,exptonegalphaT);
				if (numberofmeans>0) { //Here's where we do Butler and King A7
					double runningtotal=0;
					int chosencolumn=0;
					while (chosencolumn<Matrix1->size2) {
						runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix1,rowtaxon,chosencolumn)); //if we don't have entries, we'll be taking e^0-e^0=0
						chosencolumn++;
						runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix1,rowtaxon,chosencolumn));
						chosencolumn++;
					}
					gsl_matrix_set(W_BK_A7,rowtaxon,1,exptonegalphaT*runningtotal);
				}
				if (numberofmeans>1) { 
					double runningtotal=0;
					int chosencolumn=0;
					while (chosencolumn<Matrix2->size2) {
						runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix2,rowtaxon,chosencolumn)); 
						chosencolumn++;
						runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix2,rowtaxon,chosencolumn));
						chosencolumn++;
					}
					gsl_matrix_set(W_BK_A7,rowtaxon,2,exptonegalphaT*runningtotal);
				}
				if (numberofmeans>2) { 
					double runningtotal=0;
					int chosencolumn=0;
					while (chosencolumn<Matrix3->size2) {
						runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix3,rowtaxon,chosencolumn)); 
						chosencolumn++;
						runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix3,rowtaxon,chosencolumn));
						chosencolumn++;
					}
					gsl_matrix_set(W_BK_A7,rowtaxon,3,exptonegalphaT*runningtotal);
				}
				if (numberofmeans>3) { 
					double runningtotal=0;
					int chosencolumn=0;
					while (chosencolumn<Matrix4->size2) {
						runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix4,rowtaxon,chosencolumn)); 
						chosencolumn++;
						runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix4,rowtaxon,chosencolumn));
						chosencolumn++;
					}
					gsl_matrix_set(W_BK_A7,rowtaxon,4,exptonegalphaT*runningtotal);
				}
				if (numberofmeans>4) { 
					double runningtotal=0;
					int chosencolumn=0;
					while (chosencolumn<Matrix5->size2) {
						runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix5,rowtaxon,chosencolumn)); 
						chosencolumn++;
						runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix5,rowtaxon,chosencolumn));
						chosencolumn++;
					}
					gsl_matrix_set(W_BK_A7,rowtaxon,5,exptonegalphaT*runningtotal);
				}
				if (numberofmeans>5) { 
					double runningtotal=0;
					int chosencolumn=0;
					while (chosencolumn<Matrix6->size2) {
						runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix6,rowtaxon,chosencolumn)); 
						chosencolumn++;
						runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix6,rowtaxon,chosencolumn));
						chosencolumn++;
					}
					gsl_matrix_set(W_BK_A7,rowtaxon,6,exptonegalphaT*runningtotal);
				}
				if (numberofmeans>6) { 
					double runningtotal=0;
					int chosencolumn=0;
					while (chosencolumn<Matrix7->size2) {
						runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix7,rowtaxon,chosencolumn)); 
						chosencolumn++;
						runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix7,rowtaxon,chosencolumn));
						chosencolumn++;
					}
					gsl_matrix_set(W_BK_A7,rowtaxon,7,exptonegalphaT*runningtotal);
				}
				if (numberofmeans>7) { 
					double runningtotal=0;
					int chosencolumn=0;
					while (chosencolumn<Matrix8->size2) {
						runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix8,rowtaxon,chosencolumn)); 
						chosencolumn++;
						runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix8,rowtaxon,chosencolumn));
						chosencolumn++;
					}
					gsl_matrix_set(W_BK_A7,rowtaxon,8,exptonegalphaT*runningtotal);
				}
			}
			
			gsl_matrix *VCVfinal=gsl_matrix_calloc(ntax,ntax);
			gsl_matrix_add_constant (ScaledVCV, -1.0*gsl_matrix_min (ScaledVCV)); //replaces DeleteStem
			gsl_matrix_memcpy(VCVfinal,ScaledVCV);
			for (int r=0;r<ntax;r++) {
				gsl_matrix_set(VCVfinal,r,r,(gsl_matrix_get(ScaledVCV,r,r)+(gsl_vector_get(tipvariance,r)))); //replaces AddTipVarianceVectorToRateTImesVCV
			}			
			//gsl_matrix_memcpy(VCVfinal,brownie.AddTipVarianceVectorToRateTimesVCV(brownie.DeleteStem(ScaledVCV),tipvariance));
			gsl_vector *tipresiduals=gsl_vector_calloc(ntax);
			gsl_vector *tipexpectations=gsl_vector_calloc(ntax);
			gsl_vector * OUmeans=gsl_vector_calloc(numberofmeans+1);
			gsl_vector_set(OUmeans,0,rootmean);
			for (int position=1;position<=numberofmeans;position++) {
				gsl_vector_set(OUmeans,position,gsl_vector_get(variables,position+2));
			}
			gsl_blas_dgemv (CblasNoTrans,1, W_BK_A7, OUmeans,0, tipexpectations); 
			gsl_vector_memcpy(tipresiduals,observedtips);
			gsl_vector_sub(tipresiduals,tipexpectations);
			likelihood=(brownie.GetLScore(VCVfinal,tipresiduals,1)); //-lnL actually
			gsl_matrix_free (VCVfinal);
			gsl_vector_free (tipresiduals);
			gsl_vector_free (tipexpectations);
			gsl_vector_free (OUmeans);
			gsl_matrix_free(BranchingTimes);
			gsl_matrix_free(ScaledVCV);

		}
		//if (detailedoutput) {
		//	brownie.message="lnL = ";
		//	brownie.message+=likelihood;
		//	brownie.PrintMessage();
		//}
				gsl_matrix_free(VCVtotal);
		gsl_vector_free(tipvariance);
		gsl_vector_free(observedtips);
		gsl_vector_free(expectedtips);
		gsl_matrix_free(W_BK_A7);
		return likelihood; //-lnL actually
}

double OptimizationFnMultiModel::GetLikelihoodOUSM_FixedAttractionRate_gsl( const gsl_vector * variables, void *obj) {
	double temp;
	temp= ((OptimizationFnMultiModel*)obj)->GetLikelihoodOUSM_FixedAttractionRate(variables);
	return temp;
}

double OptimizationFnMultiModel::GetLikelihoodOUSM_FixedAttractionRate(const gsl_vector * variables) {
	//For OU models, Matrix1 consists of time between root and MRCA of two taxa (so, a typical Brownie VCV).
	//To convert this to an OU VCV (as in Butler & King, equation A5), have to use this matrix to get T-s (entry in the given matrix)
	//and time s (which is total root-tip length minus Matrix1(i,j))
	//Each subsequent matrix will have times for a different model (ie, OU mean). Each may not be square: it will have one 
	//row per taxon, and odd columns will have the start time for an interval spent with that model and even columns will have the end
	//time for the corresponding interval. These matrices can then be used to pull out times to transform using the OU model,
	//as in Butler & King A7.
	//Note that unlike Butler & King 2004 and Hansen 1997, Brownie allows multiple models per branch 
	gsl_vector *variablesoriginalpos=gsl_vector_calloc(2+variables->size);
	//cout<<"fixedparams = ";
	//for (int i=0;i<fixedparams->size;i++) {
	//	cout<<gsl_vector_get(fixedparams,i)<<"\t";
	//}
	//cout<<endl;
	for (int i=0;i<2;i++) {
		gsl_vector_set(variablesoriginalpos,i,gsl_vector_get(fixedparams,i));
	}
	for (int i=2;i<variablesoriginalpos->size;i++) {
		gsl_vector_set(variablesoriginalpos,i,gsl_vector_get(variables,i-2));
	}
	double rate=gsl_vector_get(variablesoriginalpos,0);
	double attraction=gsl_vector_get(variablesoriginalpos,1);
	int numberofmeans=-3+(variablesoriginalpos->size);
	double rootmean=gsl_vector_get(variablesoriginalpos,2);
	int ntax=Matrix0->size1;
	gsl_matrix *VCVtotal=gsl_matrix_calloc(ntax,ntax);
	gsl_vector *tipvariance=gsl_vector_calloc(ntax);
	gsl_vector *observedtips=gsl_vector_calloc(ntax);
	gsl_vector *expectedtips=gsl_vector_calloc(ntax);
	gsl_matrix *W_BK_A7=gsl_matrix_calloc(ntax,numberofmeans+1); //W matrix based on equation A7 of Butler and King
	gsl_vector_memcpy(observedtips,Vector1);
	gsl_vector_memcpy(tipvariance,Vector2);
	double likelihood;
	if (attraction<0) {
		likelihood=GSL_POSINF;
	}
	else if (rate<0) {
		likelihood=GSL_POSINF;
	}
	else if (gsl_vector_min(tipvariance)<0) {
		likelihood=GSL_POSINF;
	}
	else if (attraction>maxATTRACTIONPARAM) {
		likelihood=GSL_POSINF;
		if(detailedoutput) {
			brownie.message="Warning: OU attraction parameter of ";
			brownie.message+=attraction;
			brownie.message+=" is greater than the maximum allowed, ";
			brownie.message+=maxATTRACTIONPARAM;
			brownie.PrintMessage();
		}
	}
	else if (attraction<minATTRACTIONPARAM) {
		likelihood=GSL_POSINF;
		if (detailedoutput) {
			brownie.message="Warning: OU attraction parameter of ";
			brownie.message+=attraction;
			brownie.message+=" is less than the minimum allowed, ";
			brownie.message+=minATTRACTIONPARAM;
			brownie.PrintMessage();
		}
	}
	else {
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
			//cout<<"attraction="<<attraction<<" roottotiptime="<<roottotiptime;
//			gsl_error_handler_t * previoushandler=//gsl_set_error_handler_off();
		double exptonegalphaT=0;
		/*			int status=brownie.browniesafe_gsl_sf_exp_e10_e(-1.0*attraction*roottotiptime,exptonegalphaTresult);
		if (status) {
			if (status == GSL_ERANGE) {
				cout<<"Underflow (or overflow) -- attraction parameter too big"<<endl;
			}
		}*/
		exptonegalphaT=brownie.browniesafe_gsl_sf_exp(-1.0*attraction*roottotiptime);
			//cout<<" exptonegalphaT="<<exptonegalphaT<<endl;
//			gsl_set_error_handler (previoushandler);
		for (int rowtaxon=0;rowtaxon<ntax;rowtaxon++) {
			for (int coltaxon=0;coltaxon<ntax;coltaxon++) {
				gsl_matrix_set(ScaledVCV,rowtaxon,coltaxon,(0.5*rate/attraction)*(brownie.browniesafe_gsl_sf_exp(-2.0*attraction*(roottotiptime-gsl_matrix_get(BranchingTimes,rowtaxon,coltaxon))))*(1.0-(brownie.browniesafe_gsl_sf_exp(-2.0*attraction*(gsl_matrix_get(BranchingTimes,rowtaxon,coltaxon))))));
			}
			gsl_matrix_set(W_BK_A7,rowtaxon,0,exptonegalphaT);
			if (numberofmeans>0) { //Here's where we do Butler and King A7
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix1->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix1,rowtaxon,chosencolumn)); //if we don't have entries, we'll be taking e^0-e^0=0
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix1,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,1,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>1) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix2->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix2,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix2,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,2,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>2) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix3->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix3,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix3,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,3,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>3) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix4->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix4,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix4,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,4,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>4) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix5->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix5,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix5,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,5,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>5) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix6->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix6,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix6,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,6,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>6) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix7->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix7,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix7,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,7,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>7) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix8->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix8,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix8,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,8,exptonegalphaT*runningtotal);
			}
		}
		
		gsl_matrix *VCVfinal=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix_add_constant (ScaledVCV, -1.0*gsl_matrix_min (ScaledVCV)); //replaces DeleteStem
		gsl_matrix_memcpy(VCVfinal,ScaledVCV);
		for (int r=0;r<ntax;r++) {
			gsl_matrix_set(VCVfinal,r,r,(gsl_matrix_get(ScaledVCV,r,r)+(gsl_vector_get(tipvariance,r)))); //replaces AddTipVarianceVectorToRateTImesVCV
		}			
			//gsl_matrix_memcpy(VCVfinal,brownie.AddTipVarianceVectorToRateTimesVCV(brownie.DeleteStem(ScaledVCV),tipvariance));
		gsl_vector *tipresiduals=gsl_vector_calloc(ntax);
		gsl_vector *tipexpectations=gsl_vector_calloc(ntax);
		gsl_vector * OUmeans=gsl_vector_calloc(numberofmeans+1);
		gsl_vector_set(OUmeans,0,rootmean);
		for (int position=1;position<=numberofmeans;position++) {
			gsl_vector_set(OUmeans,position,gsl_vector_get(variablesoriginalpos,position+2));
		}
		gsl_blas_dgemv (CblasNoTrans,1, W_BK_A7, OUmeans,0, tipexpectations); 
		gsl_vector_memcpy(tipresiduals,observedtips);
		gsl_vector_sub(tipresiduals,tipexpectations);
		likelihood=(brownie.GetLScore(VCVfinal,tipresiduals,1)); //-lnL actually
		if(detailedoutput) {
			char outputstring[30];
			sprintf(outputstring,"%30.28f",likelihood);
			cout<<outputstring<<"\t";
			cout<<gsl_vector_get(fixedparams,0)<<"\t";
			cout<<gsl_vector_get(fixedparams,1)<<"\t";
			for (int i=0;i<variables->size;i++) {
				cout<<gsl_vector_get(variables,i)<<"\t";
			}
			cout<<endl;
		}
		gsl_matrix_free (VCVfinal);
		gsl_vector_free (tipresiduals);
		gsl_vector_free (tipexpectations);
		gsl_vector_free (OUmeans);
		gsl_matrix_free(BranchingTimes);
		gsl_matrix_free(ScaledVCV);
		
	}
		//if (detailedoutput) {
		//	brownie.message="lnL = ";
		//	brownie.message+=likelihood;
		//	brownie.PrintMessage();
		//}
				gsl_matrix_free(VCVtotal);
gsl_vector_free(tipvariance);
gsl_vector_free(observedtips);
gsl_vector_free(expectedtips);
gsl_matrix_free(W_BK_A7);
gsl_vector_free(variablesoriginalpos);
return likelihood; //-lnL actually
}

double OptimizationFnMultiModel::GetLikelihoodOUSM_OnlyVariablesAttractionRate_gsl( const gsl_vector * variables, void *obj) {
	double temp;
	temp= ((OptimizationFnMultiModel*)obj)->GetLikelihoodOUSM_OnlyVariablesAttractionRate(variables);
	return temp;
}

double OptimizationFnMultiModel::GetLikelihoodOUSM_OnlyVariablesAttractionRate(const gsl_vector * variables) {
	//optimizing attraction, rate, and OU means does not seem to work well (slight difference in attraction params can lead to very different optimal OU means; the search algorithm seems to get stuck on local optima with very weird OU means and pretty good attraction and rate
	//though Butler and King A8 can be used to get OU means, it doesn't seem to work, so we still have to optimize them. Basic idea here is to have an outer loop optimizing the attraction and rate parameters, with the likelihood for a given pair of them
	//being calculated by finding the best OU means with the attraction and rate params fixed.
	gsl_vector *variablesoriginalpos=gsl_vector_calloc(fixedparams->size); //fixed params also includes entries for rate and attraction, though these can actually vary
	for (int i=0;i<2;i++) {
		gsl_vector_set(fixedparams,i,gsl_vector_get(variables,i));
	}
	double likelihood=GSL_POSINF;
	gsl_vector * results=gsl_vector_calloc(-2+fixedparams->size);
	double bestlikelihood=GSL_POSINF;
	int hitlimitscount=0;
	double startingancestralstatemean=brownie.GetAncestralState(Matrix0,Vector1);
	double startingratemean=brownie.EstimateRate(Matrix0,brownie.GetTipResiduals(Vector1,startingancestralstatemean));
	double localstepsize=0.001*GSL_MAX(fabs(gsl_vector_max(Vector1)),fabs(gsl_vector_min(Vector1)));
	for (int startnum=0;startnum<randomstarts;startnum++) {
		const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
		gsl_multimin_fminimizer *s = NULL;
		gsl_vector *ss, *x;
		size_t iter = 0, i;
		int status;
		double size;
		bool hitlimits=false;
		ss = gsl_vector_alloc (-2+fixedparams->size);
		gsl_vector_set_all (ss, localstepsize); //because the initial value can be too large for some datasets
		//cout<<"fixedparams = ";
		//for (int i=0;i<fixedparams->size;i++) {
		//	cout<<gsl_vector_get(fixedparams,i)<<"\t";
		//}
		//cout<<endl;
		
		/* Starting point */
		x = gsl_vector_calloc (-2+fixedparams->size);
			for (int cell=0; cell<-2+fixedparams->size;cell++) {
				gsl_vector_set (x,cell,startingancestralstatemean+gsl_ran_gaussian(r,startingratemean));
				if(detailedoutput) {
					cout<<"Starting ancstate "<<cell<<" "<<gsl_vector_get(x,cell)<<endl;
				}
				
			}
		OptimizationFnMultiModel *pt;
		pt=(this);
		double (*F)(const gsl_vector *, void *);
		F = &OptimizationFnMultiModel::GetLikelihoodOUSM_FixedAttractionRate_gsl;
		gsl_multimin_function minex_func;
		minex_func.f=*F;
		minex_func.params=pt;
		minex_func.n = -2+fixedparams->size;
		s = gsl_multimin_fminimizer_alloc (T, -2+fixedparams->size);
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
		}
		while (status == GSL_CONTINUE && iter < maxiterations);
		if (s->fval<bestlikelihood) {
			gsl_vector_memcpy(results,s->x);
			bestlikelihood=s->fval;
		}
		int iterationsrequired=iter;
		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free (s);
		
	}
	for (int i=0;i<results->size;i++) {
		gsl_vector_set(fixedparams,i+2,gsl_vector_get(results,i));
	}
	return bestlikelihood;
}


double OptimizationFnMultiModel::GetLikelihoodOUSM_AnalyticMeans_gsl( const gsl_vector * variables, void *obj) {
	double temp;
	temp= ((OptimizationFnMultiModel*)obj)->GetLikelihoodOUSM_AnalyticMeans(variables);
	return temp;
}

double OptimizationFnMultiModel::GetLikelihoodOUSM_AnalyticMeans(const gsl_vector * variables) {
	//For OU models, Matrix1 consists of time between root and MRCA of two taxa (so, a typical Brownie VCV).
	//To convert this to an OU VCV (as in Butler & King, equation A5), have to use this matrix to get T-s (entry in the given matrix)
	//and time s (which is total root-tip length minus Matrix1(i,j))
	//Each subsequent matrix will have times for a different model (ie, OU mean). Each may not be square: it will have one 
	//row per taxon, and odd columns will have the start time for an interval spent with that model and even columns will have the end
	//time for the corresponding interval. These matrices can then be used to pull out times to transform using the OU model,
	//as in Butler & King A7.
	//Note that unlike Butler & King 2004 and Hansen 1997, Brownie allows multiple models per branch 
	
	double rate=gsl_vector_get(variables,0);
	double attraction=gsl_vector_get(variables,1);
	int numberofmeans=-3+(fixedparams->size); //does not include ancstate
	int ntax=Matrix0->size1;
	gsl_matrix *VCVtotal=gsl_matrix_calloc(ntax,ntax);
	gsl_vector *tipvariance=gsl_vector_calloc(ntax);
	gsl_vector *observedtips=gsl_vector_calloc(ntax);
	gsl_vector *expectedtips=gsl_vector_calloc(ntax);
	gsl_matrix *W_BK_A7=gsl_matrix_calloc(ntax,numberofmeans+1); //W matrix based on equation A7 of Butler and King
	gsl_vector_memcpy(observedtips,Vector1);
	gsl_vector_memcpy(tipvariance,Vector2);
	double likelihood;
	if (attraction<0) {
		likelihood=GSL_POSINF;
	}
	else if (rate<0) {
		likelihood=GSL_POSINF;
	}
	else if (gsl_vector_min(tipvariance)<0) {
		likelihood=GSL_POSINF;
	}
	else if (attraction>maxATTRACTIONPARAM) {
		likelihood=GSL_POSINF;
		if(detailedoutput) {
			brownie.message="Warning: OU attraction parameter of ";
			brownie.message+=attraction;
			brownie.message+=" is greater than the maximum allowed, ";
			brownie.message+=maxATTRACTIONPARAM;
			brownie.PrintMessage();
		}
	}
	else if (attraction<minATTRACTIONPARAM) {
		likelihood=GSL_POSINF;
		if (detailedoutput) {
			brownie.message="Warning: OU attraction parameter of ";
			brownie.message+=attraction;
			brownie.message+=" is less than the minimum allowed, ";
			brownie.message+=minATTRACTIONPARAM;
			brownie.PrintMessage();
		}
	}
	else {
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
			//cout<<"attraction="<<attraction<<" roottotiptime="<<roottotiptime;
//			gsl_error_handler_t * previoushandler=//gsl_set_error_handler_off();
		double exptonegalphaT=0;
		/*			int status=brownie.browniesafe_gsl_sf_exp_e10_e(-1.0*attraction*roottotiptime,exptonegalphaTresult);
		if (status) {
			if (status == GSL_ERANGE) {
				cout<<"Underflow (or overflow) -- attraction parameter too big"<<endl;
			}
		}*/
		exptonegalphaT=brownie.browniesafe_gsl_sf_exp(-1.0*attraction*roottotiptime);
			//cout<<" exptonegalphaT="<<exptonegalphaT<<endl;
//			gsl_set_error_handler (previoushandler);
		for (int rowtaxon=0;rowtaxon<ntax;rowtaxon++) {
			for (int coltaxon=0;coltaxon<ntax;coltaxon++) {
				gsl_matrix_set(ScaledVCV,rowtaxon,coltaxon,(0.5*rate/attraction)*(brownie.browniesafe_gsl_sf_exp(-2.0*attraction*(roottotiptime-gsl_matrix_get(BranchingTimes,rowtaxon,coltaxon))))*(1.0-(brownie.browniesafe_gsl_sf_exp(-2.0*attraction*(gsl_matrix_get(BranchingTimes,rowtaxon,coltaxon))))));
			}
			gsl_matrix_set(W_BK_A7,rowtaxon,0,exptonegalphaT);
			if (numberofmeans>0) { //Here's where we do Butler and King A7
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix1->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix1,rowtaxon,chosencolumn)); //if we don't have entries, we'll be taking e^0-e^0=0
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix1,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,1,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>1) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix2->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix2,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix2,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,2,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>2) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix3->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix3,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix3,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,3,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>3) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix4->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix4,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix4,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,4,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>4) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix5->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix5,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix5,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,5,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>5) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix6->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix6,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix6,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,6,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>6) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix7->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix7,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix7,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,7,exptonegalphaT*runningtotal);
			}
			if (numberofmeans>7) { 
				double runningtotal=0;
				int chosencolumn=0;
				while (chosencolumn<Matrix8->size2) {
					runningtotal+=brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix8,rowtaxon,chosencolumn)); 
					chosencolumn++;
					runningtotal+=-1.0*brownie.browniesafe_gsl_sf_exp(attraction*gsl_matrix_get(Matrix8,rowtaxon,chosencolumn));
					chosencolumn++;
				}
				gsl_matrix_set(W_BK_A7,rowtaxon,8,exptonegalphaT*runningtotal);
			}
		}
		
		gsl_matrix *VCVfinal=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix_add_constant (ScaledVCV, -1.0*gsl_matrix_min (ScaledVCV)); //replaces DeleteStem
		gsl_matrix_memcpy(VCVfinal,ScaledVCV);
		for (int r=0;r<ntax;r++) {
			gsl_matrix_set(VCVfinal,r,r,(gsl_matrix_get(ScaledVCV,r,r)+(gsl_vector_get(tipvariance,r)))); //replaces AddTipVarianceVectorToRateTImesVCV
		}			
			//gsl_matrix_memcpy(VCVfinal,brownie.AddTipVarianceVectorToRateTimesVCV(brownie.DeleteStem(ScaledVCV),tipvariance));
		gsl_vector *tipresiduals=gsl_vector_calloc(ntax);
		gsl_vector *tipexpectations=gsl_vector_calloc(ntax);
		
		//Here we analytically get the means, following the glssoln approach of OUCH
		gsl_matrix *ScaledVCVsansRate=gsl_matrix_calloc(ntax,ntax);
		gsl_matrix_memcpy(ScaledVCVsansRate,ScaledVCV); //Note that this means it won't work with nonzero fixed tip variance, as we've not added taht to ScaledVCV
	//	cout<<"Initial ScaledVCVsansRate"<<endl;
	//	brownie.PrintMatrix(ScaledVCVsansRate);
		gsl_matrix_scale(ScaledVCVsansRate,1.0/rate);
	//	cout<<"ScaledVCVsansRate after dividing by rate"<<endl;
	//	brownie.PrintMatrix(ScaledVCVsansRate);
		gsl_linalg_cholesky_decomp (ScaledVCVsansRate);
	//	cout<<"Chol(ScaledVCVsansRate)"<<endl;
	//	brownie.PrintMatrix(ScaledVCVsansRate);
		for(int i=0;i<ntax;i++) {
			for (int j=0;j<i;j++) {
				gsl_matrix_set(ScaledVCVsansRate,i,j,0.0);
			}
		}
	//	cout<<"Chol(ScaledVCVsansRate) deleted lower triangle"<<endl;
	//	brownie.PrintMatrix(ScaledVCVsansRate);
		gsl_matrix *vh = gsl_matrix_calloc(ScaledVCVsansRate->size2,ScaledVCVsansRate->size1);	
		gsl_matrix_transpose_memcpy (vh,ScaledVCVsansRate);
	//	cout<<"vh=(Chol(ScaledVCVsansRate))^T"<<endl;
	//	brownie.PrintMatrix(vh);
		//vh=transpose(chol(Vtilde))
		//gsl_permutation * p = gsl_permutation_alloc (ntax);
		//int signum;
		//gsl_matrix *Winvertstart=gsl_matrix_calloc(W_BK_A7->size1,W_BK_A7->size2);
		//gsl_matrix *Winvert=gsl_matrix_calloc(W_BK_A7->size1,W_BK_A7->size2);
		//gsl_matrix_memcpy (Winvertstart, W_BK_A7);
		//gsl_linalg_LU_decomp (Winvertstart,p, &signum);
		//gsl_linalg_LU_invert (Winvertstart,p, Winvert);
		//gsl_permutation_free(p);
		//gsl_matrix_free (Winvertstart);
		gsl_matrix *vhOverW = gsl_matrix_calloc(W_BK_A7->size1,W_BK_A7->size2);
		gsl_matrix *Usvd = gsl_matrix_calloc(W_BK_A7->size1,W_BK_A7->size2);
		gsl_matrix *Vsvd = gsl_matrix_calloc(W_BK_A7->size2,W_BK_A7->size2);
		gsl_vector *Dsvd = gsl_vector_calloc(W_BK_A7->size2);
		gsl_vector *workspace = gsl_vector_calloc(W_BK_A7->size2);
		//[U,D,v]=svd(vh\w)
		//A=vh\w
		/* From Matlab documentation:
			If A is an m-by-n matrix with m ~= n and B is a column vector with m components, or a matrix with several such columns,
			then X = A\B is the solution in the least squares sense to the under- or overdetermined system of equations AX = B. 
			The effective rank, k, of A, is determined from the QR decomposition with pivoting (see "Algorithm" for details).
			A solution X is computed which has at most k nonzero components per column. 
		 */
		gsl_matrix *QR=gsl_matrix_calloc(vh->size1,vh->size2);
		gsl_matrix_memcpy(QR,vh);
		gsl_vector *tau=gsl_vector_calloc(GSL_MIN(vh->size1,vh->size2));
		gsl_linalg_QR_decomp (QR,tau);
		/* from GSL documentation
			Function: int gsl_linalg_QR_decomp (gsl_matrix * A, gsl_vector * tau)
			
			This function factorizes the M-by-N matrix A into the QR decomposition A = Q R. 
			On output the diagonal and upper triangular part of the input matrix contain the matrix R. 
			The vector tau and the columns of the lower triangular part of the matrix A contain the 
			Householder coefficients and Householder vectors which encode the orthogonal matrix Q. 
			The vector tau must be of length k=\min(M,N). The matrix Q is related to these components by,
			Q = Q_k ... Q_2 Q_1 where Q_i = I - \tau_i v_i v_i^T and v_i is the Householder vector
			v_i = (0,...,1,A(i+1,i),A(i+2,i),...,A(m,i)). This is the same storage scheme as used by lapack. 
		*/
		
		
		/* from GSL documentation
			Ñ Function: int gsl_linalg_QR_solve (const gsl_matrix * QR, const gsl_vector * tau, const gsl_vector * b, gsl_vector * x)
			
			This function solves the square system A x = b using the QR decomposition of A into (QR, tau) given by gsl_linalg_QR_decomp. 
		*/
		//Since gsl can't solve A x = B (where B is  matrix), we'll do A x = B[all row,col 0], A x = B[all row, col 1], etc., or for these variables, vhOverW[all row, selected col]=vh\w[all row,selected col]
		
		for (int col=0;col<vhOverW->size2;col++) {
			gsl_vector *b=gsl_vector_calloc(vhOverW->size1);
			for (int row=0;row<vhOverW->size1;row++) {
				gsl_vector_set(b,row,gsl_matrix_get(W_BK_A7,row,col));
			}
			gsl_vector *x=gsl_vector_calloc(vhOverW->size1);
			gsl_linalg_QR_solve (QR, tau, b, x);
			for (int row=0;row<vhOverW->size1;row++) {
				gsl_matrix_set(vhOverW,row,col,gsl_vector_get(x,row));
			}
			gsl_vector_free(b);
			gsl_vector_free(x);
		}
	
		
		//store vhOverW in Usvd, as in gsl, to solve A=USV' (R is X=UDV'), object A is used to store object U		
		gsl_matrix_memcpy(Usvd,vhOverW);
		gsl_linalg_SV_decomp (Usvd, Vsvd,Dsvd, workspace);
		//double tol=GSL_MIN(stoppingprecision,1e-12); //A bit hacky
		double tol=0.0001;
		double workingtol=tol;
		//for (double workingtol=tol;workingtol<1;workingtol=workingtol*5) {
			int rcount=0;
			vector<int> bigenoughvalues;
			for (int i=0;i<Dsvd->size;i++) {
				if (gsl_vector_get(Dsvd,i)>workingtol) {
					rcount++;
					bigenoughvalues.push_back(i);
				}
			}
		//	cout<<"Tolerance is "<<workingtol<<endl<<"rcount="<<rcount<<endl;
			if (W_BK_A7->size1==0 || W_BK_A7->size2==0 || rcount==0) {
				likelihood=GSL_POSINF;
			}
			else {
			gsl_matrix *Usvdshrunk = gsl_matrix_calloc(W_BK_A7->size1,rcount);
			gsl_matrix *Vsvdshrunk = gsl_matrix_calloc(W_BK_A7->size2,rcount);
			gsl_matrix *Smatrix = gsl_matrix_calloc(rcount,rcount);
			for (int i=0;i<bigenoughvalues.size();i++) {
				for (int j=0;j<W_BK_A7->size1;j++) {
					gsl_matrix_set(Usvdshrunk,j,i,gsl_matrix_get(Usvd,j,bigenoughvalues[i]));
				}
				for (int j=0;j<W_BK_A7->size2;j++) {
					gsl_matrix_set(Vsvdshrunk,j,i,gsl_matrix_get(Vsvd,j,bigenoughvalues[i]));
				}
				gsl_matrix_set(Smatrix,i,i,1.0/gsl_vector_get(Dsvd,bigenoughvalues[i]));
			}
		//OUmeansMat=(Vsvdshrunk*(Smatrix*(Usvdshrunk)'))*(vh\tips)
		//vhOverTips
		//MatA=Smatrix*(Usvdshrunk)'
		//MatB=Vsvdshrunk*MatA
		//OUmeansMat=MatB*vhOverTips
			gsl_vector * vhOverTips = gsl_vector_calloc(W_BK_A7->size1);
		//gsl_vector *residual=gsl_vector_calloc(W_BK_A7->size1);
			gsl_linalg_QR_solve (QR, tau, observedtips, vhOverTips);
			gsl_matrix *MatA=gsl_matrix_calloc(rcount,W_BK_A7->size1);
			gsl_matrix *MatB=gsl_matrix_calloc(W_BK_A7->size2,W_BK_A7->size1);
			gsl_vector *OUmeans = gsl_vector_calloc(W_BK_A7->size2);
			gsl_blas_dgemm (CblasNoTrans,CblasTrans, 1.0, Smatrix, Usvdshrunk, 0.0, MatA);
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Vsvdshrunk, MatA, 0.0, MatB);
			gsl_blas_dgemv (CblasNoTrans, 1.0, MatB, vhOverTips, 1.0, OUmeans);
		//	cout<<"numberofmeans="<<numberofmeans<<" fixedparams="<<fixedparams->size<<" OUmeans="<<OUmeans->size<<" W_BK_A7 #col="<<W_BK_A7->size2<<endl;
		//	cout<<"attraction="<<attraction<<" rate="<<rate<<endl;
			/*			cout<<"W_BK_A7 = "<<endl;
			brownie.PrintMatrix(W_BK_A7);
			cout<<"vh="<<endl;
			brownie.PrintMatrix(vh);
			cout<<"vhOverW="<<endl;
			brownie.PrintMatrix(vhOverW);
			cout<<"Vsvd="<<endl;
			brownie.PrintMatrix(Vsvd);
			cout<<"Usvd="<<endl;
			brownie.PrintMatrix(Usvd);
			cout<<"Dsvd="<<endl;
			brownie.PrintVector(Dsvd);
			cout<<"Usvdshrunk="<<endl;
			brownie.PrintMatrix(Usvdshrunk);
			cout<<"Vsvdshrunk="<<endl;
			brownie.PrintMatrix(Vsvdshrunk);
			cout<<"Smatrix="<<endl;
			brownie.PrintMatrix(Smatrix);
			cout<<"vhOverTips="<<endl;
			brownie.PrintVector(vhOverTips);
			*/
		//	cout<<"OUmeans="<<endl;
		//	brownie.PrintVector(OUmeans);
			gsl_matrix_free(MatA);
			gsl_matrix_free(MatB);
			gsl_matrix_free(Smatrix);
			gsl_matrix_free(Vsvdshrunk);
			gsl_matrix_free(Usvdshrunk);
			gsl_vector_free(vhOverTips);
			
			
			for (int i=0;i<numberofmeans+1;i++) {
				gsl_vector_set(fixedparams,i+2,gsl_vector_get(OUmeans,i));
			}
			
			
			gsl_blas_dgemv (CblasNoTrans,1, W_BK_A7, OUmeans,0, tipexpectations); 
			gsl_vector_memcpy(tipresiduals,observedtips);
			gsl_vector_sub(tipresiduals,tipexpectations);
			likelihood=(brownie.GetLScore(VCVfinal,tipresiduals,1)); //-lnL actually
			gsl_vector_free (OUmeans);
			
		//}
			gsl_matrix_free (VCVfinal);
			gsl_vector_free (tipresiduals);
			gsl_vector_free (tipexpectations);
			gsl_matrix_free(BranchingTimes);
			gsl_matrix_free(ScaledVCV);
			gsl_matrix_free(QR);
			gsl_vector_free(tau);
			gsl_matrix_free(Vsvd);
			gsl_matrix_free(Usvd);
			gsl_matrix_free(vh);
			gsl_matrix_free(ScaledVCVsansRate);
			gsl_vector_free(workspace);
			gsl_vector_free(Dsvd);
		}
	}
		//if (detailedoutput) {
		//	brownie.message="lnL = ";
		//	brownie.message+=likelihood;
		//	brownie.PrintMessage();
		// }
	
	gsl_matrix_free(VCVtotal);
	gsl_vector_free(tipvariance);
	gsl_vector_free(observedtips);
	gsl_vector_free(expectedtips);
	gsl_matrix_free(W_BK_A7);
	return likelihood; //-lnL actually
}




////////////////////////

//double OptimizationFnDiscreteRate::GetLikelihoodOnTree_gsl( const gsl_vector * variables, void *obj) {
//	double temp;
//	temp= ((OptimizationFnDiscreteRate*)obj)->GetLikelihoodOnTree(variables);
//	return temp;
//}

//double OptimizationFnDiscreteRate::GetLikelihoodOnTree(const gsl_vector * variables) {
//	int numberofrates=variables->size;
//	int numberofstates=Matrix0->size1;
	
//	return likelihood;
//}




LindyFn::LindyFn(int Inmaxiterations, double Instoppingprecision, int Inrandomstarts, double Instepsize, bool Indetailedoutput,TreesBlock* intrees, TaxaBlock* intaxa, AssumptionsBlock* inassumptions, CharactersBlock* incharacters) :brownie()
{
	trees=intrees;
    taxa=intaxa;
    assumptions=inassumptions;
    characters=incharacters;		
    maxiterations=Inmaxiterations;
    stoppingprecision=Instoppingprecision;
    randomstarts=Inrandomstarts;
    stepsize=Instepsize;
    detailedoutput=Indetailedoutput;
	// cout<<"First entry in Matrix1 is "<<gsl_matrix_get(Matrix1,0,0)<<endl;
}



double LindyFn::GetLikelihoodUnderLindy2_gsl( const gsl_vector * variables, void *obj) 
{
	double temp;
	temp= ((LindyFn*)obj)->GetLikelihoodUnderLindy2(variables);
	return temp;
}

double LindyFn::GetLikelihoodUnderLindy2(const gsl_vector * variables)
{
	double rateA=gsl_vector_get(variables,0);
	double rateB=gsl_vector_get(variables,1);
	double likelihood=(brownie.CalculateDiscreteLindy2(rateA,rateB));
	return likelihood;
}

double LindyFn::GetLikelihoodUnderLindy1_gsl( const gsl_vector * variables, void *obj) 
{
	double temp;
	temp= ((LindyFn*)obj)->GetLikelihoodUnderLindy1(variables);
	return temp;
}

double LindyFn::GetLikelihoodUnderLindy1(const gsl_vector * variables)
{
	double rateA=gsl_vector_get(variables,0);
	double likelihood=(brownie.CalculateDiscreteLindy1(rateA));
	return likelihood;
}

gsl_vector * LindyFn::GeneralOptimization(int ChosenModel)
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
		LindyFn *pt;
		pt=(this);
		double (*F)(const gsl_vector *, void *);
		if (ChosenModel==1) {
			F = &LindyFn::GetLikelihoodUnderLindy1_gsl;
		}
		else if (ChosenModel==2) {
			F = &LindyFn::GetLikelihoodUnderLindy2_gsl;
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
		brownie.message="Replicate ";
		brownie.message+=startnum+1;
		if (hitlimits) {
			brownie.message+=" **WARNING**";
		}
		brownie.message+="\n   NM iterations needed = ";
		int iterationsrequired=iter;
		brownie.message+=iterationsrequired;
		if (hitlimits) {
			brownie.message+=" **Max iterations hit; see WARNING below**";
		}
		brownie.message+="\n   LnL = ";
		char outputstring[60];
		sprintf(outputstring,"%60.45f",-1.0*(s->fval));
		brownie.message+=outputstring;
		// brownie.message+="\n   Rate = ";
		// brownie.message+=gsl_vector_get(s->x,0);
		if (detailedoutput) {
			brownie.PrintMessage();
		}
		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free (s);
		if(detailedoutput==false) {
			brownie.ProgressBar(0);
		}
	}
	if (hitlimitscount>0) {
		brownie.message="\n----------------------------------------------------------------------------\n WARNING: Out of ";
		brownie.message+=randomstarts;
		brownie.message+=" optimization starts, ";
		brownie.message+=hitlimitscount;
		if (hitlimitscount==1) {
			brownie.message+=" was ";
		}
		else {
			brownie.message+=" were ";
		}
		brownie.message+="stopped by hitting\n  the maximum # of iterations. This means that those replicates\n  may not even have hit the local maximum.\n\n  You can increase the maximum number of iterations or decrease the\n  precision with the NumOpt command. You could also consider\n  increasing the number of random starts using that same command.\n\n  If this happened on a small proportion of replicates, though,\n  or if the precision (below) is good enough, don't worry about it.\n----------------------------------------------------------------------------";
		if (detailedoutput) {
			brownie.PrintMessage();
		}
		else if ((randomstarts-hitlimitscount)<10 && (hitlimitscount/randomstarts)>.1) {
			brownie.PrintMessage();
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
};




