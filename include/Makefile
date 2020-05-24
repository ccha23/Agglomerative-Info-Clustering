CXX=g++ 
RM=rm -rf
CPPFLAGS=-g -std=c++11 -I.

gaussian_SRCS=gaussian_demo.cpp 
hypergraph_SRCS=hypergraph_demo.cpp

gaussian_OBJS=$(subst .cpp,.o,$(gaussian_SRCS))
hypergraph_OBJS=$(subst .cpp,.o,$(hypergraph_SRCS))

all: gaussian_demo hypergraph_demo

gaussian_demo: $(gaussian_OBJS)
	$(CXX) $(CPPFLAGS) -o gaussian_demo.out $(gaussian_OBJS) 

hypergraph_demo: $(hypergraph_OBJS)
	$(CXX) $(CPPFLAGS) -o hypergraph_demo.out $(hypergraph_OBJS) 

clean:
	$(RM) $(gaussian_OBJS) $(hypergraph_OBJS) gaussian_demo.out hypergraph_demo.out

