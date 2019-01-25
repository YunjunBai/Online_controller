#
# compiler
#
CC        = g++ 
#CC       	  = clang++
CXXFLAGS 		= -Wall -Wextra -std=c++11 -fopenmp -O3 -DNDEBUG -g -o speed 
					
#CXXDEBUG		= -Wall -Wextra -std=c++11 
#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/utils


Eigen            = -I/opt/local -I/usr/include/eigen3

.PHONY: apartment

TARGET = apartment

all: $(TARGET)

apartment: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(Eigen) apartment.cc -o apartment

#simulate: 
#	$(CC) $(CXXFLAGS) $(SCOTSINC) simulate.cc -o simulate

clean:
	rm apartment

