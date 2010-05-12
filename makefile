###################################################
#
# Makefile for Brownie
# Creator [Xcode -> Makefile Ver: 3.01]
# Created: [Fri Mar 17 10:05:58 2006]
#
###################################################

#note that version.h is only updated when makefiledebug is run (to make sure that the regular makefile works when the source is not in subversion)
#
# Macros
#

CC = /usr/bin/g++
CC_OPTIONS = -m32 -fexceptions -O0 -Wno-deprecated -fpermissive -g
#added these as link options
#LNK_OPTIONS = -t -L/usr/local/lib/ -lgsl -lgslcblas -lm -L./gtp.0.15_Modified/nexus_parser/ -lnp -L./gtp.0.15_Modified/my_structures/ -lmy_structures 
LNK_OPTIONS = -m32 -t -L/usr/local/lib/ -lgsl -lgslcblas -lm


#
# INCLUDE directories for Brownie
#

#INCLUDE = -I./ncl-2.0/src/ -I./TreeLib/ -I. -I./TreeLib/gport/ -I./gtp.0.15_Modified/  -I./gtp.0.15_Modified/nexus_parser/ -I./gtp.0.15_Modified/my_structures/
INCLUDE = -I./ncl-2.0/src/ -I./TreeLib/ -I. -I./TreeLib/gport/

#INCLUDE1= -I./gtp.0.15_Modified/nexus_parser/ 
#INCLUDE2= -I./gtp.0.15_Modified/my_structures/

#LIB2 = -I./gtp.0.15_Modified/my_structures/
#LIB = -I./gtp.0.15_Modified/nexus_parser/


#
# Build Brownie
#

Brownie : \
		brownie.o\
		cdfvectorholder.o\
#		gtp.o\
#		lex.yy.o\
#		main.o\
#		mySmallTreeLib.o\
#		y.tab.o\
#		yywrap.o\
#		my_hash.o\
#		my_hash_fixed.o\
#		my_hgraph.o\
#		my_queue.o\
#		my_slist.o\
#		my_structures.o\
#		my_vector.o\
		allelesblock.o\
		assumptionsblock.o\
		charactersblock.o\
		charactersblock2.o\
		datablock.o\
		discretedatum.o\
		discretematrix.o\
		distancedatum.o\
		distancesblock.o\
		nexus.o\
		nexusblock.o\
		nexustoken.o\
		nxsdate.o\
		nxsstring.o\
		setreader.o\
		taxablock.o\
		treesblock.o\
		xnexus.o\
		gport.o\
		gtree.o\
		ntree.o\
		stree.o\
		containingtree.o\
		Parse.o\
		tokeniser.o\
		treedrawer.o\
		TreeLib.o\
		treeorder.o\
		treereader.o\
		treewriter.o\
		lcaquery.o\
		quartet.o\
		optimizationfn.o

	$(CC) $(LNK_OPTIONS) \
		brownie.o\
		cdfvectorholder.o\
#		gtp.o\
#		lex.yy.o\
#		main.o\
#		mySmallTreeLib.o\
#		y.tab.o\
#		yywrap.o\
#		my_hash.o\
#		my_hash_fixed.o\
#		my_hgraph.o\
#		my_queue.o\
#		my_slist.o\
#		my_structures.o\
#		my_vector.o\
		allelesblock.o\
		assumptionsblock.o\
		charactersblock.o\
		charactersblock2.o\
		datablock.o\
		discretedatum.o\
		discretematrix.o\
		distancedatum.o\
		distancesblock.o\
		nexus.o\
		nexusblock.o\
		nexustoken.o\
		nxsdate.o\
		nxsstring.o\
		setreader.o\
		taxablock.o\
		treesblock.o\
		xnexus.o\
		gport.o\
		gtree.o\
		ntree.o\
		stree.o\
		containingtree.o\		
		Parse.o\
		tokeniser.o\
		treedrawer.o\
		TreeLib.o\
		treeorder.o\
		treereader.o\
		treewriter.o\
		lcaquery.o\
		quartet.o\
		optimizationfn.o\
		Brownie

#
# Build the parts of Brownie
#


# Item # 1 -- brownie --
brownie.o : brownie.cpp cdfvectorholder.o
	$(CC) $(CC_OPTIONS) brownie.cpp -c $(INCLUDE) -o brownie.o 

#don't optimize cdf vector holder for speed (which takes a loooooooooong time compiling, probably doesn't improve speed); instead optimize it for size
cdfvectorholder.o : cdfvectorholder.cpp
	$(CC) $(CC_OPTIONS) -fexceptions cdfvectorholder.cpp -c $(INCLUDE) -o cdfvectorholder.o 

#gtp.o: 
#	cd gtp.0.15_Modified; make clean; make; mv gtp.o ../gtp.o

#my_hash.o:
#	cd gtp.0.15_Modified/my_structures; make

#my_hash_fixed.o:
#	cd gtp.0.15_Modified/my_structures; make

#my_hgraph.o:
#	cd gtp.0.15_Modified/my_structures; make

#my_queue.o:
#	cd gtp.0.15_Modified/my_structures; make

#my_slist.o:
#	cd gtp.0.15_Modified/my_structures; make

#my_structures.o:
#	cd gtp.0.15_Modified/my_structures; make

#my_vector.o:
#	cd gtp.0.15_Modified/my_structures; make

#lex.yy.o:
#	cd gtp.0.15_Modified/nexus_parser; make

#main.o:
#	cd gtp.0.15_Modified/nexus_parser; make

#mySmallTreeLib.o:
#	cd gtp.0.15_Modified/nexus_parser; make

#y.tab.o:
#	cd gtp.0.15_Modified/nexus_parser; make

#yywrap.o:
#	cd gtp.0.15_Modified/nexus_parser; make
	
# Item #1b -- optimizationfn --
optimizationfn.o : optimizationfn.cpp
	$(CC) $(CC_OPTIONS) optimizationfn.cpp -c $(INCLUDE) -o optimizationfn.o

# Item # 2 -- allelesblock --
allelesblock.o : ./ncl-2.0/src/allelesblock.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/allelesblock.cpp -c $(INCLUDE) -o allelesblock.o


# Item # 3 -- assumptionsblock --
assumptionsblock.o : ./ncl-2.0/src/assumptionsblock.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/assumptionsblock.cpp -c $(INCLUDE) -o assumptionsblock.o


# Item # 4 -- charactersblock --
charactersblock.o : ./ncl-2.0/src/charactersblock.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/charactersblock.cpp -c $(INCLUDE) -o charactersblock.o


# Item # 4b -- charactersbloc2k --
charactersblock2.o : ./ncl-2.0/src/charactersblock2.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/charactersblock2.cpp -c $(INCLUDE) -o charactersblock2.o


# Item # 5 -- datablock --
datablock.o : ./ncl-2.0/src/datablock.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/datablock.cpp -c $(INCLUDE) -o datablock.o


# Item # 6 -- discretedatum --
discretedatum.o : ./ncl-2.0/src/discretedatum.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/discretedatum.cpp -c $(INCLUDE) -o discretedatum.o


# Item # 7 -- discretematrix --
discretematrix.o : ./ncl-2.0/src/discretematrix.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/discretematrix.cpp -c $(INCLUDE) -o discretematrix.o


# Item # 8 -- distancedatum --
distancedatum.o : ./ncl-2.0/src/distancedatum.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/distancedatum.cpp -c $(INCLUDE) -o distancedatum.o


# Item # 9 -- distancesblock --
distancesblock.o : ./ncl-2.0/src/distancesblock.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/distancesblock.cpp -c $(INCLUDE) -o distancesblock.o


# Item # 10 -- nexus --
nexus.o : ./ncl-2.0/src/nexus.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/nexus.cpp -c $(INCLUDE) -o nexus.o


# Item # 11 -- nexusblock --
nexusblock.o : ./ncl-2.0/src/nexusblock.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/nexusblock.cpp -c $(INCLUDE) -o nexusblock.o


# Item # 12 -- nexustoken --
nexustoken.o : ./ncl-2.0/src/nexustoken.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/nexustoken.cpp -c $(INCLUDE) -o nexustoken.o


# Item # 13 -- nxsdate --
nxsdate.o : ./ncl-2.0/src/nxsdate.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/nxsdate.cpp -c $(INCLUDE) -o nxsdate.o


# Item # 14 -- nxsstring --
nxsstring.o : ./ncl-2.0/src/nxsstring.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/nxsstring.cpp -c $(INCLUDE) -o nxsstring.o


# Item # 15 -- setreader --
setreader.o : ./ncl-2.0/src/setreader.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/setreader.cpp -c $(INCLUDE) -o setreader.o


# Item # 16 -- taxablock --
taxablock.o : ./ncl-2.0/src/taxablock.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/taxablock.cpp -c $(INCLUDE) -o taxablock.o


# Item # 17 -- treesblock --
treesblock.o : ./ncl-2.0/src/treesblock.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/treesblock.cpp -c $(INCLUDE) -o treesblock.o


# Item # 18 -- xnexus --
xnexus.o : ./ncl-2.0/src/xnexus.cpp
	$(CC) $(CC_OPTIONS) ./ncl-2.0/src/xnexus.cpp -c $(INCLUDE) -o xnexus.o


# Item # 19 -- gport --
gport.o : ./TreeLib/gport/gport.cpp
	$(CC) $(CC_OPTIONS) ./TreeLib/gport/gport.cpp -c $(INCLUDE) -o gport.o


# Item # 20 -- gtree --
gtree.o : ./TreeLib/gtree.cpp
	$(CC) $(CC_OPTIONS) ./TreeLib/gtree.cpp -c $(INCLUDE) -o gtree.o


# Item # 21 -- ntree --
ntree.o : ./TreeLib/ntree.cpp
	$(CC) $(CC_OPTIONS) ./TreeLib/ntree.cpp -c $(INCLUDE) -o ntree.o

# Item # 21b -- stree --
stree.o : ./TreeLib/stree.cpp
	$(CC) $(CC_OPTIONS) ./TreeLib/stree.cpp -c $(INCLUDE) -o stree.o

# Item # 21c -- containingtree --
containingtree.o : ./TreeLib/containingtree.cpp
	$(CC) $(CC_OPTIONS) ./TreeLib/containingtree.cpp -c $(INCLUDE) -o containingtree.o

# Item # 22 -- Parse --
Parse.o : ./TreeLib/Parse.cpp
	$(CC) $(CC_OPTIONS) ./TreeLib/Parse.cpp -c $(INCLUDE) -o Parse.o


# Item # 23 -- tokeniser --
tokeniser.o : ./TreeLib/tokeniser.cpp
	$(CC) $(CC_OPTIONS) ./TreeLib/tokeniser.cpp -c $(INCLUDE) -o tokeniser.o


# Item # 24 -- treedrawer --
treedrawer.o : ./TreeLib/treedrawer.cpp
	$(CC) $(CC_OPTIONS) ./TreeLib/treedrawer.cpp -c $(INCLUDE) -o treedrawer.o


# Item # 25 -- TreeLib --
TreeLib.o : ./TreeLib/TreeLib.cpp
	$(CC) $(CC_OPTIONS) ./TreeLib/TreeLib.cpp -c $(INCLUDE) -o TreeLib.o


# Item # 26 -- treeorder --
treeorder.o : ./TreeLib/treeorder.cpp
	$(CC) $(CC_OPTIONS) ./TreeLib/treeorder.cpp -c $(INCLUDE) -o treeorder.o


# Item # 27 -- treereader --
treereader.o : ./TreeLib/treereader.cpp
	$(CC) $(CC_OPTIONS) ./TreeLib/treereader.cpp -c $(INCLUDE) -o treereader.o


# Item # 28 -- treewriter --
treewriter.o : ./TreeLib/treewriter.cpp
	$(CC) $(CC_OPTIONS) ./TreeLib/treewriter.cpp -c $(INCLUDE) -o treewriter.o


# Item # 29 -- lcaquery --
lcaquery.o:  ./TreeLib/lcaquery.cpp
	$(CC) $(CC_OPTIONS) ./TreeLib/lcaquery.cpp -c $(INCLUDE) -o lcaquery.o

# Item # 30 -- quartet --
quartet.o: ./TreeLib/quartet.cpp 
	$(CC) $(CC_OPTIONS) ./TreeLib/quartet.cpp -c $(INCLUDE) -o quartet.o


# FINAL ITEM -- BROWNIE
Brownie : brownie.o cdfvectorholder.o optimizationfn.o allelesblock.o assumptionsblock.o charactersblock.o charactersblock2.o datablock.o discretedatum.o discretematrix.o distancedatum.o distancesblock.o nexus.o nexusblock.o nexustoken.o nxsdate.o nxsstring.o setreader.o taxablock.o treesblock.o xnexus.o gport.o gtree.o ntree.o stree.o containingtree.o Parse.o tokeniser.o treedrawer.o TreeLib.o treeorder.o treereader.o treewriter.o lcaquery.o quartet.o
	$(CC) brownie.o cdfvectorholder.o optimizationfn.o allelesblock.o assumptionsblock.o charactersblock.o charactersblock2.o datablock.o discretedatum.o discretematrix.o distancedatum.o distancesblock.o nexus.o nexusblock.o nexustoken.o nxsdate.o nxsstring.o setreader.o taxablock.o treesblock.o xnexus.o gport.o gtree.o ntree.o stree.o containingtree.o Parse.o tokeniser.o treedrawer.o TreeLib.o treeorder.o treereader.o treewriter.o lcaquery.o quartet.o $(LNK_OPTIONS) -o brownie
	@echo ""
	@chmod a+x brownie
	@echo "brownie has now been compiled. yippee."

macintel : brownie.o cdfvectorholder.o optimizationfn.o allelesblock.o assumptionsblock.o charactersblock.o charactersblock2.o datablock.o discretedatum.o discretematrix.o distancedatum.o distancesblock.o nexus.o nexusblock.o nexustoken.o nxsdate.o nxsstring.o setreader.o taxablock.o treesblock.o xnexus.o gport.o gtree.o ntree.o stree.o containingtree.o Parse.o tokeniser.o treedrawer.o TreeLib.o treeorder.o treereader.o treewriter.o lcaquery.o quartet.o
	$(CC) -arch i386 brownie.o cdfvectorholder.o optimizationfn.o allelesblock.o assumptionsblock.o charactersblock.o charactersblock2.o datablock.o discretedatum.o discretematrix.o distancedatum.o distancesblock.o nexus.o nexusblock.o nexustoken.o nxsdate.o nxsstring.o setreader.o taxablock.o treesblock.o xnexus.o gport.o gtree.o ntree.o stree.o containingtree.o Parse.o tokeniser.o treedrawer.o TreeLib.o treeorder.o treereader.o treewriter.o lcaquery.o quartet.o $(LNK_OPTIONS) -o brownie
	@echo ""
	@chmod a+x brownie
	@echo "brownie has now been compiled. yippee."

macppc : brownie.o cdfvectorholder.o optimizationfn.o allelesblock.o assumptionsblock.o charactersblock.o charactersblock2.o datablock.o discretedatum.o discretematrix.o distancedatum.o distancesblock.o nexus.o nexusblock.o nexustoken.o nxsdate.o nxsstring.o setreader.o taxablock.o treesblock.o xnexus.o gport.o gtree.o ntree.o stree.o containingtree.o Parse.o tokeniser.o treedrawer.o TreeLib.o treeorder.o treereader.o treewriter.o lcaquery.o quartet.o
	$(CC) -arch ppc brownie.o cdfvectorholder.o optimizationfn.o allelesblock.o assumptionsblock.o charactersblock.o charactersblock2.o datablock.o discretedatum.o discretematrix.o distancedatum.o distancesblock.o nexus.o nexusblock.o nexustoken.o nxsdate.o nxsstring.o setreader.o taxablock.o treesblock.o xnexus.o gport.o gtree.o ntree.o stree.o containingtree.o Parse.o tokeniser.o treedrawer.o TreeLib.o treeorder.o treereader.o treewriter.o lcaquery.o quartet.o $(LNK_OPTIONS) -o brownie
	@echo ""
	@chmod a+x brownie
	@echo "brownie has now been compiled. yippee."

macppc64 : brownie.o cdfvectorholder.o optimizationfn.o allelesblock.o assumptionsblock.o charactersblock.o charactersblock2.o datablock.o discretedatum.o discretematrix.o distancedatum.o distancesblock.o nexus.o nexusblock.o nexustoken.o nxsdate.o nxsstring.o setreader.o taxablock.o treesblock.o xnexus.o gport.o gtree.o ntree.o stree.o containingtree.o Parse.o tokeniser.o treedrawer.o TreeLib.o treeorder.o treereader.o treewriter.o lcaquery.o quartet.o
	$(CC) -arch ppc64 brownie.o cdfvectorholder.o optimizationfn.o allelesblock.o assumptionsblock.o charactersblock.o charactersblock2.o datablock.o discretedatum.o discretematrix.o distancedatum.o distancesblock.o nexus.o nexusblock.o nexustoken.o nxsdate.o nxsstring.o setreader.o taxablock.o treesblock.o xnexus.o gport.o gtree.o ntree.o stree.o containingtree.o Parse.o tokeniser.o treedrawer.o TreeLib.o treeorder.o treereader.o treewriter.o lcaquery.o quartet.o $(LNK_OPTIONS) -o brownie
	@echo ""
	@chmod a+x brownie
	@echo "brownie has now been compiled. yippee."

clean: 
	-rm *.o


##### END RUN ####
