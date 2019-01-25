#
# compiler
#
CC        = g++ 
#CC       	  = clang++
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG -g -o speed -fopenmp
#CXXDEBUG		= -Wall -Wextra -std=c++11 
#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/utils


Eigen            = -I/opt/local -I/usr/include/eigen3

.PHONY: manipulator

TARGET = manipulator

all: $(TARGET)

manipulator: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(Eigen) manipulator.cc -o manipulator

#simulate: 
#	$(CC) $(CXXFLAGS) $(SCOTSINC) simulate.cc -o simulate

clean:
	rm manipulator

