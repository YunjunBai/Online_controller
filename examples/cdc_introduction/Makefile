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


Eigen            =  -I/opt/local -I/usr/include/Eigen

.PHONY: vehicle

TARGET = vehicle

all: $(TARGET)

vehicle: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(Eigen) vehicle.cc -o vehicle

#simulate: 
#	$(CC) $(CXXFLAGS) $(SCOTSINC) simulate.cc -o simulate

clean:
	rm vehicle

