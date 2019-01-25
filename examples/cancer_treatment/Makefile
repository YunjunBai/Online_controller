

#
# compiler
#
CC        = g++ 
#CC       	  = clang++
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG -g -fopenmp
#CXXDEBUG		= -Wall -Wextra -std=c++11 
#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/utils


Eigen            = -I/opt/local -I/usr/include/eigen3

.PHONY: cancer_treatment

TARGET = cancer_treatment

all: $(TARGET)

cancer_treatment: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(Eigen) cancer_treatment.cc -o cancer_treatment

#simulate: 
#	$(CC) $(CXXFLAGS) $(SCOTSINC) simulate.cc -o simulate

clean:
	rm cancer_treatment

