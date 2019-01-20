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


Eigen            = -I /opt/local

.PHONY: quadcopter_12d

TARGET = quadcopter_12d

all: $(TARGET)

quadcopter_12d: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(Eigen) quadcopter_12d.cc -o quadcopter_12d

#simulate: 
#	$(CC) $(CXXFLAGS) $(SCOTSINC) simulate.cc -o simulate

clean:
	rm quadcopter_12d

