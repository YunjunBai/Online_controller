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

.PHONY: quadcopter_12d_test

TARGET = quadcopter_12d_test

all: $(TARGET)

quadcopter_12d_test: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(Eigen) quadcopter_12d_test.cc -o quadcopter_12d_test

#simulate: 
#	$(CC) $(CXXFLAGS) $(SCOTSINC) simulate.cc -o simulate

clean:
	rm quadcopter_12d_test

