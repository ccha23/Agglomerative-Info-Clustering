CXX=g++ 
RM=rm -rf
CPPFLAGS=-g -std=c++11 -I./include

gaussian_SRCS=gaussian_demo.cpp 
hypergraph_SRCS=hypergraph_demo.cpp
fls_SRCS=fls_demo.cpp

gaussian_OBJS=$(subst .cpp,.o,$(gaussian_SRCS))
hypergraph_OBJS=$(subst .cpp,.o,$(hypergraph_SRCS))
fls_OBJS=$(subst .cpp,.o,$(fls_SRCS))

all: gaussian_demo hypergraph_demo fls_demo

gaussian_demo: $(gaussian_OBJS)
	$(CXX) $(CPPFLAGS) -o gaussian_demo.out $(gaussian_OBJS) 

hypergraph_demo: $(hypergraph_OBJS)
	$(CXX) $(CPPFLAGS) -o hypergraph_demo.out $(hypergraph_OBJS) 

fls_demo: $(fls_OBJS)
	$(CXX) $(CPPFLAGS) -o fls_demo.out $(fls_OBJS) -lpari

clean:
	$(RM) $(gaussian_OBJS) $(hypergraph_OBJS) $(fls_OBJS) gaussian_demo.out hypergraph_demo.out fls_demo.out