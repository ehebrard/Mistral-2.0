

COMPILFLAGS = -D_UNIX -D_BIT32 
COPTIMIZE = -O3

MAINDIR = .

CCC = g++
#CCC = g++-5

# BOOSTDIR = /Users/Shared/boost_1_78_0

XCSP3DIR=$(MAINDIR)/XCSP3-CPP-Parser


BIN=$(MAINDIR)/bin
SRC=$(MAINDIR)/src/lib
MOD=$(MAINDIR)/examples
OBJ=$(MAINDIR)/src/obj
INC=$(MAINDIR)/src/include
DOC=$(MAINDIR)/doc
TCL=$(MAINDIR)/tools/tclap/include
XINC=$(XCSP3DIR)/include
XSRC=$(XCSP3DIR)/src
XOBJ=$(XCSP3DIR)/obj
XLIB=$(XCSP3DIR)/lib

XLIBFLAG = `xml2-config --libs`
XINCFLAG = `xml2-config --cflags`  -I$(XINC)
XCPP_FILES := $(wildcard $(XSRC)/*.cc)
XOBJ_FILES := $(addprefix $(XOBJ)/,$(notdir $(XCPP_FILES:.cc=.o)))


CFLAGS = -Wall -std=c++11 -I$(INC) -I$(TCL) 


MODELS = $(wildcard $(MOD)/src/*.cpp)
BINS = $(patsubst $(MOD)/src/%, $(BIN)/%, $(MODELS:.cpp=))


PINCSRC = $(wildcard $(INC)/*.hpp)
PLIBSRC = $(wildcard $(SRC)/*.cpp)
PLIBAUX = $(PLIBSRC:.cpp=.o)
PLIBOBJ = $(patsubst $(SRC)/%, $(OBJ)/%, $(PLIBAUX))
XLIBSRC = $(wildcard $(XSRC)/*.cc)
XLIBAUX = $(XLIBSRC:.cc=.o)
XLIBOBJ = $(patsubst $(XSRC)/%, $(XOBJ)/%, $(XLIBAUX))



LFLAGS = -L$(OBJ)


#------------------------------------------------------------
#  make all      : to compile the examples.
#------------------------------------------------------------


all: lib xcsp3 $(BINS)
	
# xcsp3: $(XCSP3DIR)/lib/libparserxcsp3core.a
# 	make MistralXCSP
#
# # $(XCSP3DIR)/lib/libparserxcsp3core.a:
# # 	cd $(XCSP3DIR)/samples/; make lib


# The library
lib: $(PLIBOBJ) $(PUTIOBJ)
$(OBJ)/%.o:  $(SRC)/%.cpp $(INC)/%.hpp
	@echo 'compile '$<
	$(CCC) $(CFLAGS) -c $< -o $@ 

# The examples
$(BIN)/%: $(MOD)/obj/%.o $(PLIBOBJ)
	@echo 'link '$<
	$(CCC) $(CFLAGS) $(PLIBOBJ) $< -lm -o $@

$(MOD)/obj/%.o: $(MOD)/src/%.cpp
	@echo 'compile '$<
	$(CCC) $(CFLAGS) -c $< -o $@ 
	
# MistralXCSP: $(MOD)/obj/MistralXCSP.o $(XCSP3DIR)/lib/libparserxcsp3core.a
# 	@echo 'link '$<
# 	$(CCC) $(CFLAGS) $(PLIBOBJ) -L $(XCSP3DIR)/lib `xml2-config --libs` -lparserxcsp3core $< -lm -o $(BIN)/$@
	
# Examples, one at a time
%: $(MOD)/obj/%.o $(PLIBOBJ) 
	@echo 'link '$<	
	$(CCC) $(CFLAGS) $(PLIBOBJ) $< -lm -o $(BIN)/$@ 
	
	
xcsp3: $(BIN)/MistralXCSP

testlib: $(XCSP3DIR)/samples/main.cc
	$(CCC) $(CFLAGS) $(XINCFLAG) -o $(BIN)/testlib main.cc -L $(XLIB) $(XLIBFLAG) -lparserxcsp3core

$(BIN)/MistralXCSP: $(MOD)/obj/MistralXCSP.o $(PLIBOBJ) $(XLIB)/libparserxcsp3core.a
	@echo 'link '$<
	$(CCC) $(CFLAGS) $(XINCFLAG) $(PLIBOBJ) -L $(XLIB) $(XLIBFLAG) -lparserxcsp3core $< -lm -o $@

$(MOD)/obj/MistralXCSP.o: $(MOD)/src/MistralXCSP.cpp
	@echo 'compile '$<
	$(CCC) $(CFLAGS) $(XINCFLAG) -c $< -o $@ 
	
$(XOBJ)/%.o: $(XSRC)/%.cc
	@mkdir -p $(XOBJ)
	$(CCC) $(CFLAGS) $(XINCFLAG) -c -o $@ $<

$(XLIB)/libparserxcsp3core.a: $(XOBJ_FILES) 
	@mkdir -p $(XLIB)
	ar rcsv $(XLIB)/libparserxcsp3core.a $(XOBJ_FILES)


# xclean:
# 	rm -rf $(XOBJ)/*.o $(XLIB)/* $(BIN)/testlib
#


DATE := $(shell date '+%y-%m-%d')

clean : 
	cd fz;	make clean;
	# if [ -d "./XCSP3-CPP-Parser/samples" ]; then cd XCSP3-CPP-Parser/samples ; make clean; echo exists; fi
	rm -rf $(XOBJ)/*.o $(XLIB)/* $(BIN)/testlib $(OBJ)/*.o $(OBJ)/*.a $(SRC)/*~ $(MOD)/obj/*.o $(MOD)/src/*~ $(MOD)/src/*/*~ $(INC)/*~ $(UTI)/*~  *~ $(BIN)/* $(DOC)/*~ ./fzn-mistral fz/mistral-fzn

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
	
	

