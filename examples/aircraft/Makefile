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

.PHONY: aircraft

TARGET = aircraft

all: $(TARGET)

aircraft: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(Eigen) aircraft.cc -o aircraft

#simulate: 
#	$(CC) $(CXXFLAGS) $(SCOTSINC) simulate.cc -o simulate

clean:
	rm aircraft

