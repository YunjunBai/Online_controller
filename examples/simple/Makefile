#
# compiler
#
CC        = g++
#CC       	  = clang++
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG  -g
CXXDEBUG		= -g -Wall -Wextra -std=c++11 
#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/utils

.PHONY: simple simple_2

TARGET = simple simple_2

all: $(TARGET)

simple: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) simple.cc -o simple

simple_2: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) simple_2.cc -o simple_2

clean:
	rm simple simple_2
