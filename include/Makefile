CXX=g++ 
RM=rm -f
CPPFLAGS=-g -std=c++11 -I.

SRCS=AIC.cpp SF.cpp PSP.cpp main.cpp
OBJS=$(subst .cc,.o,$(SRCS))

all: main

main: $(OBJS)
	$(CXX) $(CPPFLAGS) -o main.out $(OBJS) 

%.o: %.cpp
	$(CXX) $(CPPFLAGS) -c $<

clean:
	$(RM) $(OBJS)

