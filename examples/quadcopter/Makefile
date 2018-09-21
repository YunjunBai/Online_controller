#
# compiler
#
CC        = g++
#CC       	  = clang++
CXXFLAGS 		= -g -Wall -Wextra -std=c++11
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG


#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/utils


.PHONY: quadcopter 

TARGET = quadcopter 

all: $(TARGET)

quadcopter: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) quadcopter.cc -o quadcopter


clean:
	rm quadcopter 
