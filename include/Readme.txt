The source code requires Eigen C++ library, which is included in the subfolder Eigen and also available at
http://eigen.tuxfamily.org/

demo.cpp is an example that performs the exact and approximate clustering algorithm on an AWGN model.
Compile using
g++ -g -std=c++11 -I. -o demo.out demo.cpp