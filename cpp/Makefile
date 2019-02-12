CXX=g++ 
RM=rm -f
CPPFLAGS=-g -std=c++11
LDFLAGS=-g
LDLIBS=

SRCS=AIC.cpp SF.cpp PSP.cpp main.cpp
OBJS=$(subst .cc,.o,$(SRCS))

all: main

main: $(OBJS)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) -o main.out $(OBJS) $(LDLIBS) 

clean:
	$(RM) $(OBJS)

