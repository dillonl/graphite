# =========================================
# MOSAIK Banded Smith-Waterman Makefile
# (c) 2009 Michael Stromberg & Wan-Ping Lee
# =========================================

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= smithwaterman.cpp BandedSmithWaterman.cpp SmithWatermanGotoh.cpp Repeats.cpp LeftAlign.cpp IndelAllele.cpp
OBJECTS= $(SOURCES:.cpp=.o) disorder.o
OBJECTS_NO_MAIN= disorder.o BandedSmithWaterman.o SmithWatermanGotoh.o Repeats.o LeftAlign.o IndelAllele.o

# ----------------
# compiler options
# ----------------

CFLAGS=-Wall -O3
LDFLAGS=-Wl,-s
#CFLAGS=-g
PROGRAM=smithwaterman
LIBS=

all: $(PROGRAM) sw.o

.PHONY: all

disorder.o: disorder.c disorder.h
	g++ -c -o disorder.o disorder.c

sw.o: $(OBJECTS_NO_MAIN)
	ld -r $(OBJECTS_NO_MAIN) -o sw.o

$(PROGRAM): $(OBJECTS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(LDFLAGS) $(CFLAGS) -o $@ $^ $(LIBS)

.PHONY: clean

clean:
	@echo "Cleaning up."
	@rm -f *.o $(PROGRAM) *~
