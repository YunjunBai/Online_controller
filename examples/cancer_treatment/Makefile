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


.PHONY: cancer_treatment 

TARGET = cancer_treatment 

all: $(TARGET)

cancer_treatment: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) cancer_treatment.cc -o cancer_treatment


clean:
	rm cancer_treatment 
