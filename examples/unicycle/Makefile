#
# compiler
#
CC        = g++
#CC       	  = clang++
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG -g
#CXXDEBUG		= -Wall -Wextra -std=c++11 
#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/utils

.PHONY: unicycle

TARGET = unicycle

all: $(TARGET)

unicycle: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) unicycle.cc -o unicycle

#simulate: 
#	$(CC) $(CXXFLAGS) $(SCOTSINC) simulate.cc -o simulate

clean:
	rm unicycle 
