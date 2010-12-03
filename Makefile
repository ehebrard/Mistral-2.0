

OPTFLAGS = -O3 #-m32
#OPTFLAGS = -g 

COMPILFLAGS = -Wall -D_UNIX -D_BIT32 #-D_STATIC_CAST #-D_DEBUG_SEARCH -D_DEBUG_AC -D_DEBUG_PROPAG #-D_DELTA 

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