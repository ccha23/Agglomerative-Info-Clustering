The source code requires Eigen C++ library, which is included in the subfolder Eigen and also available at
http://eigen.tuxfamily.org/

gaussian_demo.cpp is an example that performs the exact and approximate clustering algorithm on an AWGN model.
Compile using
g++ -g -std=c++11 -I. -o gaussian_demo.out gaussian_demo.cpp

hypergraph_demo.cpp is an example that performs the exact and approximate clustering algorithm on a hypergraphical source model.
Compile using
g++ -g -std=c++11 -I. -o hypergraph_demo.out hypergraph_demo.cpp