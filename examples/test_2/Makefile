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

.PHONY: test_2

TARGET = test_2

all: $(TARGET)

test_2: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) test_2.cc -o test_2

#simulate: 
#	$(CC) $(CXXFLAGS) $(SCOTSINC) simulate.cc -o simulate

clean:
	rm test_2
