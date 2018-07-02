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

.PHONY: test

TARGET = test

all: $(TARGET)

test: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) test.cc -o test

#simulate: 
#	$(CC) $(CXXFLAGS) $(SCOTSINC) simulate.cc -o simulate

clean:
	rm test 
