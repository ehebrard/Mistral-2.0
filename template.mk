
MAINDIR ?= .

COPTIMIZE ?= -O3

CCC = g++
#CCC       = g++-5

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


default: xcsp3 flatzinc 
#default: flatzinc 

flatzinc:
	cd fz; make
	cp fz/mistral-fz ./bin/fzn-mistral
	
xcsp3: $(XCSP3PDIR)
	cd $(XCSP3PDIR)/samples/; make
	make MistralXCSP

$(XCSP3PDIR):
	git clone https://github.com/xcsp3team/XCSP3-CPP-Parser.git

#fz/mistral-fz: fz/mistral-fzn
#	cp fz/mistral-fzn ./bin

#fz/mistral-fzn: 
#	cd fz; make

all: lib $(BINS) flatzinc

# The library
lib: $(PLIBOBJ) $(PUTIOBJ)
$(OBJ)/%.o:  $(SRC)/%.cpp $(INC)/%.hpp
	@echo 'compile '$<
	$(CCC) $(CFLAGS) -c $< -o $@ 

# The examples
$(BIN)/%: $(MOD)/obj/%.o $(PLIBOBJ)
	@echo 'link '$<
	$(CCC) $(CFLAGS) $(PLIBOBJ) $(XLIBOBJ) $< -lm -o $@

$(MOD)/obj/%.o: $(MOD)/src/%.cpp
	@echo 'compile '$<
	$(CCC) $(CFLAGS) -c $< -o $@ 

# Examples, one at a time
%: $(MOD)/obj/%.o $(PLIBOBJ) $(XLIBOBJ)
	@echo 'link '$<	
	$(CCC) $(CFLAGS) $(PLIBOBJ) $(XLIBOBJ) `xml2-config --libs` $< -lm -o $(BIN)/$@ 

