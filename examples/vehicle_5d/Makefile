#
# compiler
#
CC        = g++ 
#CC       	  = clang++
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG -g -o speed
#CXXDEBUG		= -Wall -Wextra -std=c++11 
#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/utils


Eigen            = -I/opt/local

.PHONY: vehicle_5d

TARGET = vehicle_5d

all: $(TARGET)

vehicle_5d: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(Eigen) vehicle_5d.cc -o vehicle_5d

#simulate: 
#	$(CC) $(CXXFLAGS) $(SCOTSINC) simulate.cc -o simulate

clean:
	rm vehicle_5d

