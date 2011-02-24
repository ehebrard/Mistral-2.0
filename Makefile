
DATE := $(shell date '+%y-%m-%d')


OPTFLAGS = -O3 #-m32
#OPTFLAGS = -g 

#COMPILFLAGS = -Wall -D_UNIX -D_BIT32 -D_DEBUG_SEARCH #-D_DEBUG_AC -D_DEBUG_PROPAG #-D_DELTA 
COMPILFLAGS = -Wall -D_UNIX -D_BIT32 -DNDEBUG #-D_DEBUG_SEARCH #-D_DEBUG_NOGOOD #-D_CHRONOLOGICAL #-D_DEBUG_AC -D_DEBUG_PROPAG

CCC = g++ $(OPTFLAGS) $(COMPILFLAGS)

BIN=bin
SRC=src/lib
MOD=examples
OBJ=src/obj
INC=src/include
DOC=doc

CFLAGS = -I$(INC) 
LFLAGS = -L$(OBJ)


MODELS = $(wildcard $(MOD)/src/*.cpp)
BINS = $(patsubst $(MOD)/src/%, $(BIN)/%, $(MODELS:.cpp=))


PINCSRC = $(wildcard $(INC)/*.hpp)
PLIBSRC = $(wildcard $(SRC)/*.cpp)
PLIBAUX = $(PLIBSRC:.cpp=.o)
PLIBOBJ = $(patsubst $(SRC)/%, $(OBJ)/%, $(PLIBAUX))


#------------------------------------------------------------
#  make all      : to compile the examples.
#------------------------------------------------------------

all: lib $(PLIBOBJ) $(BINS)

clean : 
	rm -rf $(OBJ)/*.o $(OBJ)/*.a $(SRC)/*~ $(MOD)/obj/*.o $(MOD)/src/*~ $(MOD)/src/*/*~ $(INC)/*~ *~ $(BIN)/* $(DOC)/*~

# The library
lib: $(PLIBOBJ)
$(OBJ)/%.o:  $(SRC)/%.cpp $(INC)/%.hpp
	@echo 'compile '$<
	$(CCC) $(CFLAGS) -c $< -o $@ 

# The examples
$(BIN)/%: $(MOD)/obj/%.o $(PLIBOBJ)
	@echo 'link '$<
	$(CCC) $(CFLAGS)   $(PLIBOBJ) $< -lm -o $@

$(MOD)/obj/%.o: $(MOD)/src/%.cpp
	@echo 'compile '$<
	$(CCC) $(CFLAGS) -c $< -o $@ 

# Examples, one at a time
%: $(MOD)/obj/%.o $(PLIBOBJ)
	@echo 'link '$<	
	$(CCC) $(CFLAGS)   $(PLIBOBJ) $< -lm -o $(BIN)/$@ 

release: 
	@echo Export Mistral version 2.0.$(DATE)
	rm -rf Mistral-2.0.$(DATE).bz2
	rm -rf Mistral-2.0.$(DATE)
	mkdir Mistral-2.0.$(DATE)
	git archive master --format=tar | tar -x -C ./Mistral-2.0.$(DATE)
	mkdir ./Mistral-2.0.$(DATE)/bin
	mkdir ./Mistral-2.0.$(DATE)/examples/obj
	mkdir ./Mistral-2.0.$(DATE)/src/obj
	tar -cvjf Mistral-2.0.$(DATE).bz2 Mistral-2.0.$(DATE)
	scp Mistral-2.0.$(DATE).bz2 4c60.ucc.ie:/home/ehebrard/tmp/
