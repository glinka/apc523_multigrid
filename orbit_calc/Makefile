SRCS=main.cc
OBJECTS=$(SRCS:.cc=.o)

CXX = g++

CXXFLAGS = -g -Wall -std=c++0x -O3

all: integrate

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $<

integrate: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^ > ./.depend

clean:
	$(RM) *.o 

include .depend
