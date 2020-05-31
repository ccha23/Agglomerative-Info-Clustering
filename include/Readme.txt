The source code requires Eigen C++ library, which is included in the subfolder Eigen and also available at
http://eigen.tuxfamily.org/

gaussian_demo.cpp is an example that performs the exact and approximate clustering algorithm on an AWGN model. Run using 
make gaussian_demo; ./gaussian_demo.out

hypergraph_demo.cpp is an example that performs the exact and approximate clustering algorithm on a hypergraphical source model. Run using
make hypergraph_demo; ./hypergraph_demo.out

fls_demo.cpp is an example that performs the exact and approximate clustering algorithm on a finite linear source model. It requires the C library pari, which can be installed on ubuntu using
apt-get install libpari-dev
Run demo using
make fls_demo; ./fls_demo.out
