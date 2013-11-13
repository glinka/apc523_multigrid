SRCS=main.cc
OBJECTS=$(SRCS:.cc=.o)

CXX = g++

CXXFLAGS = -g -Wall -std=c++0x #-O3

all: jacobi

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $<

jacobi: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^ > ./.depend

clean:
	$(RM) *.o 

include .depend
