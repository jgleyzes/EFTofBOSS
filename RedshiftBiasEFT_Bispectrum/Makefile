CXX = g++
CXXFlags = -Wall -Wextra -pedantic -g -02

EXEC =	RedshiftBiasEFT 
SRC = $(wildcard *.cpp)
OBJ = $(SRC: .cpp=.o) 
STD = -std=c++11
GSL = -lm -lgsl -lgslcblas
CUBA = -lcuba -lm
FFTW = -lfftw3


all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $^ -o $@ $(STD) $(GSL) $(CUBA) $(FFTW)

%.o: %.cpp
	$(CXX) $(CXXFlags) -o $@ -c $^ 

.PHONY: clean

clean:
	rm $(EXEC) *~

