CXX=g++ 
RM=rm -rf
CPPFLAGS=-g -std=c++11 -I.

SRCS=
OBJS=$(subst .cc,.o,$(SRCS))

all: demo

demo: $(OBJS)
	$(CXX) $(CPPFLAGS) -o demo.out $(OBJS) 

clean:
	$(RM) $(OBJS) *.dSYM

