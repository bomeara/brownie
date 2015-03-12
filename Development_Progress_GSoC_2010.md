# Google Summer of Code Project #

---

## Objective ##
The main goal of this project will be the integration of important phylogenetic tools for ancestral state reconstruction and rate testing (via Brownie) into a commonly used scripting environment (R).  As R is a cross-platform scripting environment for which there are currently external packages to read-in and manipulate genetic and phylogenetic data.  R is an ideal place to port the functions contained in Brownie due to its robust linkage to C/C++ (innately and through external packages) and because of its potential to act as a staging area for all steps in phylogenetic data analysis, from reading in genetic data to running analyses to plotting the results.


## Project Executors ##
  * [Conrad Stack ](http://www.conradstack.info)
  * [Luke Harmon ](http://www.webpages.uidaho.edu/~lukeh/)
  * [Brian O'Meara ](http://www.brianomeara.info/)

## Final thoughts ##
Anyone please feel free to contact me (conrad.stack@gmail) with any suggestions for improvement or general comments.
<br><br>
Also check the <a href='http://www.conradstack.info'>project blog</a> for exciting technical updates and random thoughts, displaying in near-real time!<br>
<br>
<br><br>



<h1>Development progress (updated weekly)</h1>
Note:  The initial project proposal can be seen <a href='http://docs.google.com/Doc?docid=0AcvZ2biv57wNZHpnZmI0c18xMDJjNHo2Y3JmaA&hl=en'>here</a>.  The development progress described below will roughly follow this, but may deviate sharply if our strategy changes regarding how the implementation should be done.<br>
<br>
<hr />

<h3>05/05/2010</h3>
<ul><li>Tasks:<br>
</li></ul><blockquote>- Created this page</blockquote>

<hr />

<h3>05/17/2010</h3>
<ul><li>Tasks:<br>
<ul><li>Setup build environment for brownie and related R-package on windows, linux, and mac boxes<br>
</li><li>Background reading <code>[1]</code>
</li><li>Setup <a href='http://www.conradstack.info'>project blog</a>
</li><li>Step through brownie to get familiar with project structure and coding style. <code>[2]</code></li></ul></li></ul>


<ul><li>Links:<br>
<ol><li><a href='http://www.amazon.com/Computational-Molecular-Evolution-Oxford-Ecology/dp/0198567022/ref=sr_1_1?ie=UTF8&s=books&qid=1273644680&sr=8-1'>Computational Molecular Evolution</a>
</li><li><a href='http://bodegaphylo.wikispot.org/Morphological_Diversification_and_Rates_of_Evolution'>Brownie Tutorial</a>
<hr /></li></ol></li></ul>


<h3>05/24/2010</h3>
<ul><li>Add Rpack to brownie source branch<br>
<ul><li>Add configuration script which copies current Brownie build into Rpack, builds R package<br>
</li></ul></li><li>Add phylobase, Rcpp requirements to Rpack<br>
<ul><li>R functions which checks for singleton nodes and convert them to SIMMAP format (supported by Brownie)<br>
</li><li>R functions which take sufficient data (tree, taxsets, data, characters) and outputs them to a Brownie-compatible nexus format<br>
<ul><li>NOTE: This suite of data will be streamed to Brownie<br>
</li></ul></li></ul></li><li>Stubs in Brownie to calculate confidence intervals<br>
<ul><li>Discuss with Brian how these should be calculated<br>
</li><li>Determine input/output types<br>
</li></ul></li><li>Build Brownie source under R compilation process with Rpack subdir<br>
<ul><li>Add new configuration requirements<br>
</li><li>More streamlined Makefile.in and Makefile.win for compiling Brownie source.<br>
<hr /></li></ul></li></ul>


<h3>05/31/2010</h3>
<ul><li>R Brownie compilation:<br>
<ul><li>Complete Brownie compilation through R (excluding newer confidence interval tests)<br>
</li><li>Hard code tests (and test files) which run basic brownie functions (testing for runtime errors)<br>
</li></ul></li><li>R-Brownie interface:<br>
<ul><li>R functions which take sufficient data (tree, taxsets, data, characters) and outputs them to a Brownie-compatible nexus format<br>
</li><li>Pass information from phylobase to C++<br>
</li><li>Pass tree, data from C++ to phylobase<br>
<ul><li>Conversion might be able to be done by Rcpp using as() function, but could also be done in R<br>
</li></ul></li><li>Stubs for desired functions in R<br>
<ul><li>Associated stubs for documentation<br>
</li></ul></li></ul></li><li>Brownie additions:<br>
<ul><li>Begin extensive runtime tests of confidence interval tests</li></ul></li></ul>

<hr />


<h3>06/07/2010</h3>
<ul><li>R Brownie compilation:<br>
<ul><li>Hard code tests for confidence interval finding.<br>
</li></ul></li><li>R-Brownie interface:<br>
<ul><li>Define stubs for functions that will need to included in R package<br>
<ul><li>Discuss final list with Luke, Brian<br>
</li><li>Initial ideas:<br>
<ul><li>Read in or specify character data<br>
</li><li>Read in from csv or Excel file<br>
</li><li>Read in from Nexus file<br>
</li><li>Read in or specify set of trees<br>
</li><li>Read in set of posterior trees (with weights!) from BEAST, MrBayes, etc.<br>
</li><li>Read in or specify directly the nexus files, SIMMAP files to be used<br>
</li><li>Collate all necessary elements into temporary nexus file<br>
</li><li>Specify rate model<br>
</li><li>Specify output file<br>
</li><li>Check for necessary elements (taxa, characters, trees, assumptions)<br>
</li><li>Run ratetest (censored, uncensored)<br>
</li><li>Run ancestral state reconstruction<br>
</li><li>Generic functions (S4 style)<br>
</li></ul></li></ul></li><li>Generate documentation files (latex)<br>
</li></ul></li><li>Brownie additions<br>
<ul><li>Finish any debugging that is necessary<br>
<hr /></li></ul></li></ul>



<h3>06/14/2010</h3>
<ul><li>Debugging of confidence interval code:<br>
<ul><li>Confidence interval comparisons with other phylogenetic packages (e.g. HYPHY, BEAST)<br>
</li></ul></li><li>R-Brownie interface:<br>
<ul><li>C++ classes for handling data passage between R package and libBrownie<br>
</li><li>Begin R coding to pass data to C++ interface and retrieve results<br>
</li><li>Unit tests for each major function<br>
</li><li>Document (latex) as progress is made on each function<br>
</li></ul></li><li>R front-end<br>
<ul><li>Plan out list of methods and classes to add.<br>
</li><li>Plan out which features of ape and/or phylobase to use<br>
<hr /></li></ul></li></ul>




<h3>06/21/2010</h3>
<ul><li>Debugging of confidence interval code:<br>
<ul><li>Confidence interval comparisons with other phylogenetic packages (e.g. HYPHY, BEAST)<br>
</li></ul></li><li>R-Brownie interface:<br>
<ul><li>Thorough testing of the interface functions<br>
</li><li>Choose examples to eventually become vignettes<br>
</li></ul></li><li>R front-end<br>
<ul><li>Begin writing method and class code<br>
</li><li>Start latex documentation for each function<br>
<hr /></li></ul></li></ul>



<h3>06/28/2010</h3>
<ul><li>Debugging of confidence interval code:<br>
<ul><li>Confidence interval comparisons with other phylogenetic packages (e.g. HYPHY, BEAST)<br>
</li></ul></li><li>R-Brownie interface:<br>
<ul><li>Should be nearly complete, but will still have bugs<br>
</li></ul></li><li>R front-end<br>
<ul><li>Continue writing method and class code<br>
</li><li>Continue documentation<br>
</li><li>Begin vignettes<br>
</li></ul></li><li>Start and complete formal mid-term evaluation<br>
<hr /></li></ul>


<h3>07/5/2010</h3>
<ul><li>Start/Finish prototype of confidence interval additions<br>
</li><li>Streamline brownie interface<br>
<ul><li>Add more needed functions to dlInterface class<br>
</li><li>Rework package to utilize Rcpp Modules (remove wrapper from dlInterface)<br>
</li></ul></li><li>Add stubs for ALL R functions that will be included in first release<br>
</li><li>Started midterm review<br>
</li><li>Work on Latex documentation<br>
<hr /></li></ul>


