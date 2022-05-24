

COMPILFLAGS = -D_UNIX -D_BIT32 
COPTIMIZE = -O3

MAINDIR = .

CCC = g++


# BOOSTDIR = /Users/Shared/boost_1_78_0

XCSP3PDIR=$(MAINDIR)/XCSP3-CPP-Parser


BIN=$(MAINDIR)/bin
SRC=$(MAINDIR)/src/lib
MOD=$(MAINDIR)/examples
OBJ=$(MAINDIR)/src/obj
INC=$(MAINDIR)/src/include
DOC=$(MAINDIR)/doc
TCL=$(MAINDIR)/tools/tclap/include
XINC=$(XCSP3PDIR)/include
XSRC=$(XCSP3PDIR)/src
XOBJ=$(XCSP3PDIR)/obj

MODELS = $(wildcard $(MOD)/src/*.cpp)
BINS = $(patsubst $(MOD)/src/%, $(BIN)/%, $(MODELS:.cpp=))


PINCSRC = $(wildcard $(INC)/*.hpp)
PLIBSRC = $(wildcard $(SRC)/*.cpp)
PLIBAUX = $(PLIBSRC:.cpp=.o)
PLIBOBJ = $(patsubst $(SRC)/%, $(OBJ)/%, $(PLIBAUX))
XLIBSRC = $(wildcard $(XSRC)/*.cc)
XLIBAUX = $(XLIBSRC:.cc=.o)
XLIBOBJ = $(patsubst $(XSRC)/%, $(XOBJ)/%, $(XLIBAUX))


CFLAGS = -Wall -std=c++11 -I$(INC) -I$(TCL) -I$(XINC) `xml2-config --cflags`
LFLAGS = -L$(OBJ)


## Compile options
%.o:			CFLAGS +=$(COPTIMIZE)  $(COMPILFLAGS) #-ggdb -D DEBUG
%.op:			CFLAGS +=$(COPTIMIZE) -pg -ggdb -D NDEBUG
%.od:			CFLAGS +=-O0 -ggdb -D DEBUG -D INVARIANTS #-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
%.or:			CFLAGS +=$(COPTIMIZE) -D NDEBUG
%.oc:                   CFLAGS +=-O0 -fprofile-arcs -ftest-coverage -ggdb -D DEBUG


#------------------------------------------------------------
#  make all      : to compile the examples.
#------------------------------------------------------------


all: lib xcsp3 $(BINS)
	
xcsp3: $(XCSP3PDIR)/lib/libparserxcsp3core.a
	make MistralXCSP

$(XCSP3PDIR)/lib/libparserxcsp3core.a:
	cd $(XCSP3PDIR)/samples/; make lib


# The library
lib: $(PLIBOBJ) $(PUTIOBJ)
$(OBJ)/%.o:  $(SRC)/%.cpp $(INC)/%.hpp
	@echo 'compile '$<
	$(CCC) $(CFLAGS) -c $< -o $@ 

# The examples
$(BIN)/%: $(MOD)/obj/%.o $(PLIBOBJ)
	@echo 'link '$<
	$(CCC) $(CFLAGS) $(PLIBOBJ) -L $(XCSP3PDIR)/lib `xml2-config --libs` -lparserxcsp3core $< -lm -o $@

$(MOD)/obj/%.o: $(MOD)/src/%.cpp
	@echo 'compile '$<
	$(CCC) $(CFLAGS) -c $< -o $@ 
	
# MistralXCSP: $(MOD)/obj/MistralXCSP.o $(XCSP3PDIR)/lib/libparserxcsp3core.a
# 	@echo 'link '$<
# 	$(CCC) $(CFLAGS) $(PLIBOBJ) -L $(XCSP3PDIR)/lib `xml2-config --libs` -lparserxcsp3core $< -lm -o $(BIN)/$@
	
# Examples, one at a time
%: $(MOD)/obj/%.o $(PLIBOBJ) 
	@echo 'link '$<	
	$(CCC) $(CFLAGS) $(PLIBOBJ) -L $(XCSP3PDIR)/lib `xml2-config --libs` -lparserxcsp3core $< -lm -o $(BIN)/$@ 


DATE := $(shell date '+%y-%m-%d')

clean : 
	cd fz;	make clean;
	cd XCSP3-CPP-Parser/samples ; make clean; 
	rm -rf $(OBJ)/*.o $(OBJ)/*.a $(SRC)/*~ $(MOD)/obj/*.o $(MOD)/src/*~ $(MOD)/src/*/*~ $(INC)/*~ $(UTI)/*~  *~ $(BIN)/* $(DOC)/*~ ./fzn-mistral fz/mistral-fzn

archive: 
	@echo Export Mistral version 2.0.$(DATE)
	rm -rf Mistral-2.0.$(DATE).tar.gz
	rm -rf Mistral-2.0.$(DATE)
	mkdir Mistral-2.0.$(DATE)
	git archive master --format=tar | tar -x -C ./Mistral-2.0.$(DATE)
	mkdir ./Mistral-2.0.$(DATE)/globals
	#cp -rf ./fz/mznlib/*.mzn ./Mistral-2.0.$(DATE)/globals 
	tar -cvzf Mistral-2.0.$(DATE).tar.gz Mistral-2.0.$(DATE)

test: bin/unit_test
	cd fz; make; 
	bin/unit_test
	cd fz; ./test_fcts.sh passed
