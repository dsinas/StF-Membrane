EXEC := stfMem

SDIR := src
SRC := $(SDIR)/main.cpp

LDIR := src/fftsrc
LSRC := $(LDIR)/fftw++.cc

OBJ := $(LSRC:.cc=.o) $(SRC:.cpp=.o)

CXX := g++
CXXFLAGS := -O2 -Wall
LIB := -lfftw3 -lm

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ $(LIB)

#$(LDIR)%.o: %.cc
$(SDIR)%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LIB)

.PHONY: clean

clean:
	@rm -f $(SDIR)/*.o
	@rm -f $(LDIR)/*.o
	@echo "clean obj done!"

clean_data:
	@rm -rf data
	@echo "clean data done!"
