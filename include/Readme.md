The source code requires Eigen C++ library, which is included in the subfolder Eigen and also available at
http://eigen.tuxfamily.org/

cart_demo.cpp is an example that performs the exact and approximate clustering algorithm on an AWGN model.
Compile using cmake:

after installing opencv, use to following command to compile and run the demo

```
mkdir build 
cd build 
cmake .. 
make 
./AIC ../../data/test.csv 
```


``` 
test.csv
1, 2, 3, 4
1, 2, 3, 4
1, 2, 3, 4
1, 2, 6, 5
1, 2, 6, 5
1, 2, 6, 5
1, 2, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
4, 3, 6, 5
```

``` 
expected output
-- Configuring done
-- Generating done
-- Build files have been written to: /Users/handasontam/Documents/OpenSource/Agglomerative-Info-Clustering/include/build
Scanning dependencies of target AIC
[ 50%] Building CXX object CMakeFiles/AIC.dir/cart_demo.cpp.o
[100%] Linking CXX executable AIC
[100%] Built target AIC
1  2  3  4  
1  2  3  4  
1  2  3  4  
1  2  6  5  
1  2  6  5  
1  2  6  5  
1  2  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
4  3  6  5  
Info-clustering by CL tree approximation:
for B = {1 }
cart entropy: 0.174818

for B = {0 }
cart entropy: 1.57336

for B = {1 0 }
cart entropy: 1.57336

for B = {2 }
cart entropy: 0.786681

for B = {0 }
cart entropy: 1.57336

for B = {2 0 }
cart entropy: 1.57336

for B = {2 }
cart entropy: 0.786681

for B = {1 }
cart entropy: 0.174818

for B = {2 1 }
cart entropy: 0.174818

for B = {3 }
cart entropy: 0.0874089

for B = {0 }
cart entropy: 1.57336

for B = {3 0 }
cart entropy: 1.57336

for B = {3 }
cart entropy: 0.0874089

for B = {1 }
cart entropy: 0.174818

for B = {3 1 }
cart entropy: 0.174818

for B = {3 }
cart entropy: 0.0874089

for B = {2 }
cart entropy: 0.786681

for B = {3 2 }
cart entropy: 0.907068

critical values : [ 0.786681 0.786681 0.0874089 ]
partition at threshold 0.786681:[ [ 0 ] [ 1 ] [ 2 ] [ 3 ] ]
partition at threshold 0.786681:[ [ 0 ] [ 1 2 ] [ 3 ] ]
partition at threshold 0.0874089:[ [ 0 1 2 ] [ 3 ] ]
Agglomerative info-clustering:
for B = {0 }
cart entropy: 1.57336

for B = {1 }
cart entropy: 0.174818

for B = {2 }
cart entropy: 0.786681

for B = {3 }
cart entropy: 0.0874089

for B = {0 1 }
cart entropy: 1.57336

for B = {2 }
cart entropy: 0.786681

for B = {0 2 }
cart entropy: 1.57336

for B = {0 1 2 }
cart entropy: 2.125

for B = {2 }
cart entropy: 0.786681

for B = {0 2 }
cart entropy: 1.57336

for B = {0 1 2 }
cart entropy: 2.125

for B = {3 }
cart entropy: 0.0874089

for B = {0 3 }
cart entropy: 1.57336

for B = {0 1 3 }
cart entropy: 2.125

for B = {0 1 2 3 }
cart entropy: 2.125

for B = {3 }
cart entropy: 0.0874089

for B = {2 3 }
cart entropy: 0.907068

for B = {2 0 3 }
cart entropy: 1.57336

for B = {2 0 1 3 }
cart entropy: 2.125

for B = {3 }
cart entropy: 0.0874089

for B = {0 3 }
cart entropy: 1.57336

for B = {0 2 3 }
cart entropy: 1.57336

for B = {0 2 1 3 }
cart entropy: 2.125

for B = {0 2 }
cart entropy: 1.57336

for B = {1 }
cart entropy: 0.174818

for B = {3 }
cart entropy: 0.0874089

for B = {0 2 1 }
cart entropy: 2.125

for B = {3 }
cart entropy: 0.0874089

for B = {0 2 3 }
cart entropy: 1.57336

for B = {0 2 1 3 }
cart entropy: 2.125

for B = {3 }
cart entropy: 0.0874089

for B = {0 2 3 }
cart entropy: 1.57336

for B = {0 2 1 3 }
cart entropy: 2.125

for B = {0 2 3 }
cart entropy: 1.57336

for B = {1 }
cart entropy: 0.174818

for B = {0 2 3 1 }
cart entropy: 2.125

critical values : [ 0.786681 0.0874089 -0.376819 ]
partition at threshold 0.786681:[ [ 0 ] [ 1 ] [ 2 ] [ 3 ] ]
partition at threshold 0.0874089:[ [ 0 2 ] [ 1 ] [ 3 ] ]
partition at threshold -0.376819:[ [ 0 2 3 ] [ 1 ] ]
```


``` 
# for test2.csv
./AIC ../../data/test2.csv

2  2  4  4  
2  2  4  4  
2  2  4  4  
2  2  5  5  
2  2  5  5  
2  2  5  5  
2  2  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
3  3  5  5  
Info-clustering by CL tree approximation:
for B = {1 }
cart entropy: 0.174818

for B = {0 }
cart entropy: 0.174818

for B = {1 0 }
cart entropy: 0.174818

for B = {2 }
cart entropy: 0.0874089

for B = {0 }
cart entropy: 0.174818

for B = {2 0 }
cart entropy: 0.174818

for B = {2 }
cart entropy: 0.0874089

for B = {1 }
cart entropy: 0.174818

for B = {2 1 }
cart entropy: 0.174818

for B = {3 }
cart entropy: 0.0874089

for B = {0 }
cart entropy: 0.174818

for B = {3 0 }
cart entropy: 0.174818

for B = {3 }
cart entropy: 0.0874089

for B = {1 }
cart entropy: 0.174818

for B = {3 1 }
cart entropy: 0.174818

for B = {3 }
cart entropy: 0.0874089

for B = {2 }
cart entropy: 0.0874089

for B = {3 2 }
cart entropy: 0.207796

critical values : [ 0.174818 0.0874089 ]
partition at threshold 0.174818:[ [ 0 ] [ 1 ] [ 2 ] [ 3 ] ]
partition at threshold 0.0874089:[ [ 0 1 ] [ 2 ] [ 3 ] ]
Agglomerative info-clustering:
for B = {0 }
cart entropy: 0.174818

for B = {1 }
cart entropy: 0.174818

for B = {2 }
cart entropy: 0.0874089

for B = {3 }
cart entropy: 0.0874089

for B = {0 1 }
cart entropy: 0.174818

for B = {2 }
cart entropy: 0.0874089

for B = {0 2 }
cart entropy: 0.174818

for B = {0 1 2 }
cart entropy: 0.238947

for B = {2 }
cart entropy: 0.0874089

for B = {1 2 }
cart entropy: 0.174818

for B = {1 0 2 }
cart entropy: 0.238947

for B = {2 }
cart entropy: 0.0874089

for B = {1 2 }
cart entropy: 0.174818

for B = {1 0 2 }
cart entropy: 0.238947

for B = {3 }
cart entropy: 0.0874089

for B = {0 3 }
cart entropy: 0.174818

for B = {0 1 3 }
cart entropy: 0.238947

for B = {0 1 2 3 }
cart entropy: 0.238947

for B = {3 }
cart entropy: 0.0874089

for B = {1 3 }
cart entropy: 0.174818

for B = {1 0 3 }
cart entropy: 0.238947

for B = {1 0 2 3 }
cart entropy: 0.238947

for B = {3 }
cart entropy: 0.0874089

for B = {1 3 }
cart entropy: 0.174818

for B = {1 0 3 }
cart entropy: 0.238947

for B = {1 0 2 3 }
cart entropy: 0.238947

for B = {0 1 }
cart entropy: 0.174818

for B = {2 }
cart entropy: 0.0874089

for B = {3 }
cart entropy: 0.0874089

for B = {0 1 2 }
cart entropy: 0.238947

for B = {3 }
cart entropy: 0.0874089

for B = {0 1 3 }
cart entropy: 0.238947

for B = {0 1 2 3 }
cart entropy: 0.238947

for B = {3 }
cart entropy: 0.0874089

for B = {2 3 }
cart entropy: 0.207796

for B = {2 0 1 3 }
cart entropy: 0.238947

for B = {3 }
cart entropy: 0.0874089

for B = {0 1 3 }
cart entropy: 0.238947

for B = {0 1 2 3 }
cart entropy: 0.238947

critical values : [ 0.174818 0.0553443 ]
partition at threshold 0.174818:[ [ 0 ] [ 1 ] [ 2 ] [ 3 ] ]
partition at threshold 0.0553443:[ [ 0 1 ] [ 2 ] [ 3 ] ]

```
