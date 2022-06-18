CC = gcc
CXX = g++
LIBS = -lrt
CFLAGS = --std=c++11 -W -Wall -Wno-sign-compare -O2 -s -pipe -mmmx -msse -msse2 -msse3 -mcmodel=large
TARGET = ./all

all: shortestPath SOR SRH BR IPM Statistic 
shortestPath: monitor.h shortestPath.cpp
	$(CXX) $(CFLAGS) -o shortestPath shortestPath.cpp $(LIBS) -D WATCH_MEM

SOR: monitor.h SOR.cpp
	$(CXX) $(CFLAGS) -o SOR SOR.cpp $(LIBS) -D WATCH_MEM

SRH: monitor.h SRH.cpp
	$(CXX) $(CFLAGS) -o SRH SRH.cpp $(LIBS) -D WATCH_MEM

BR: monitor.h BR.cpp
	$(CXX) $(CFLAGS) -o BR BR.cpp $(LIBS) -D WATCH_MEM

IPM: monitor.h IPM.cpp
	$(CXX) $(CFLAGS) -o IPM IPM.cpp $(LIBS) -D WATCH_MEM

Statistic: monitor.h Statistic.cpp
	$(CXX) $(CFLAGS) -o Statistic Statistic.cpp $(LIBS) -D WATCH_MEM

.PHONY: clean
clean:
		-@rm *.o *.gcno *~ 2> /dev/null || true
		-@rm SOR shortestPath SRH Statistic 2> /dev/null || true
