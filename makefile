# Name of your program:

PROGRAM1 = createResponse_noempty_pp
PROGRAM2 = unfoldData_noempty_pp
PROGRAM3 = unfoldDataUncertainties_noempty_pp
PROGRAM4 = makeRawDeltaPhi
PROGRAM5 = makeHerwig_hist
PROGRAM6 = createEfficiencies

# Your source file:
SRC1 = createResponse_noempty_pp.cxx
SRC2 = unfoldData_noempty_pp.cxx
SRC3 = unfoldDataUncertainties_noempty_pp.cxx
SRC4 = makeRawDeltaPhi.cxx
SRC5 = makeHerwig_hist.cxx
SRC6 = createEfficiencies.cxx

CLASS = dijetfinder.cc
# Compiler:
CXX = clang++

# ROOT flags:
ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs)

ROOUNFOLD_SRC = /Users/daniel/sPHENIX/ppg08/roounfold/RooUnfold-master/src/

# RooUnfold install location (adjust if needed)
# Common Homebrew paths:
ROOUNFOLD_INC = /usr/local/include
ROOUNFOLD_LIB = /usr/local/lib

# Add the RooUnfold library:
ROOUNFOLD_LIBS = -lRooUnfold

# Final compiler flags:
CXXFLAGS = -O2 -Wall $(ROOTCFLAGS) -I$(ROOUNFOLD_INC) -I$(ROOUNFOLD_SRC)
LDFLAGS  = $(ROOTLIBS) -L$(ROOUNFOLD_LIB) $(ROOUNFOLD_LIBS)

OBJ1 = dijetfinder.o createResponse_noempty_pp.o
OBJ2 = dijetfinder.o unfoldData_noempty_pp.o
OBJ3 = dijetfinder.o unfoldDataUncertainties_noempty_pp.o
OBJ4 = dijetfinder.o makeRawDeltaPhi.o
OBJ5 = dijetfinder.o makeHerwig_hist.o	
OBJ6 = dijetfinder.o createEfficiencies.o

# Build rule:

all: $(PROGRAM1) $(PROGRAM2) $(PROGRAM3) $(PROGRAM4) $(PROGRAM5) $(PROGRAM6)

$(PROGRAM1): $(OBJ1)
	$(CXX) $(OBJ1) -o $(PROGRAM1) $(LDFLAGS)

$(PROGRAM2): $(OBJ2)
	$(CXX) $(OBJ2) -o $(PROGRAM2) $(LDFLAGS)

$(PROGRAM3): $(OBJ3)
	$(CXX) $(OBJ3) -o $(PROGRAM3) $(LDFLAGS)

$(PROGRAM4): $(OBJ4)
	$(CXX) $(OBJ4) -o $(PROGRAM4) $(LDFLAGS)

$(PROGRAM5): $(OBJ5)
	$(CXX) $(OBJ5) -o $(PROGRAM5) $(LDFLAGS)

$(PROGRAM6): $(OBJ6)
	$(CXX) $(OBJ6) -o $(PROGRAM6) $(LDFLAGS)

# Compile dijetfinder
dijetfinder.o: $(CLASS)
	$(CXX) $(CXXFLAGS) -c $(CLASS) -o dijetfinder.o

# Compile main source
createResponse_noempty_pp.o: $(SRC1)
	$(CXX) $(CXXFLAGS) -c $(SRC1) -o createResponse_noempty_pp.o

unfoldData_noempty_pp.o: $(SRC2)
	$(CXX) $(CXXFLAGS) -c $(SRC2) -o unfoldData_noempty_pp.o

unfoldDataUncertainties_noempty_pp.o: $(SRC3)
	$(CXX) $(CXXFLAGS) -c $(SRC3) -o unfoldDataUncertainties_noempty_pp.o

makeRawDeltaPhi.o: $(SRC4)
	$(CXX) $(CXXFLAGS) -c $(SRC4) -o makeRawDeltaPhi.o

makeHerwig_hist.o: $(SRC5)
	$(CXX) $(CXXFLAGS) -c $(SRC5) -o makeHerwig_hist.o

createEfficiencies.o: $(SRC6)
	$(CXX) $(CXXFLAGS) -c $(SRC6) -o createEfficiencies.o

# Cleanup:
clean:
	rm -f $(OBJ1) $(PROGRAM1) $(OBJ2) $(PROGRAM2) $(OBJ3) $(PROGRAM3) $(OBJ4) $(PROGRAM4) $(OBJ5) $(PROGRAM5) $(OBJ6) $(PROGRAM6) *.o
