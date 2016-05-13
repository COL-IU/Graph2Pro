# Makefile
#
#objects := $(patsubst %.cpp, %.o, $(wildcard *.cpp))
objects0 :=  smallapp.o seq.o 
objects1 :=  soap.o FastG.o seqgraph.o
objects2 :=  Graph2Pro.o
objects3 :=  GraphPep2Pro.o

all: DBGraph2Pro DBGraphPep2Pro

CXX := g++

CXXFLAGs :=

BUILD_TYPE :=

ifeq ($(BUILD_TYPE), debug)
	CXXFLAGS += -g
else
	CXXFLAGS += -O2 -w
endif

ifeq ($(pg), 1)
	CXXFLAGS += -pg
else
	CXXFLAGS +=
endif

INC := -I ./

DBGraph2Pro: DBGraph2Pro_main.cpp $(objects0) $(objects1) $(objects2)
	$(CXX) $(CXXFLAGS) -o DBGraph2Pro DBGraph2Pro_main.cpp $(objects0) $(objects1) $(objects2) $(INC) 

DBGraphPep2Pro: DBGraphPep2Pro_main.cpp $(objects0) $(objects1) $(objects3)
	$(CXX) $(CXXFLAGS) -o DBGraphPep2Pro DBGraphPep2Pro_main.cpp $(objects0) $(objects1) $(objects3) $(INC) 

$(objects0) : %.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(INC)
$(objects1) : %.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(INC)
$(objects2) : %.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(INC)
$(objects3) : %.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(INC)

clean:
	$(RM) *.o DBGraph2Pro DBGraphPep2Pro
