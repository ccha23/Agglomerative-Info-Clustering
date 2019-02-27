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
...
...
...
...
critical values : [ 0.299173 0.235086 0.110689 ]
partition at threshold 0.299173:[ [ 0 ] [ 1 ] [ 2 ] [ 3 ] ]
partition at threshold 0.235086:[ [ 0 ] [ 1 ] [ 2 3 ] ]
partition at threshold 0.110689:[ [ 0 2 3 ] [ 1 ] ]
```


``` 
# for test2.csv
./AIC ../../data/test2.csv

...
...
...
...
critical values : [ 0.174818 0.0874089 0.0232797 ]
partition at threshold 0.174818:[ [ 0 ] [ 1 ] [ 2 ] [ 3 ] ]
partition at threshold 0.0874089:[ [ 0 1 ] [ 2 ] [ 3 ] ]
partition at threshold 0.0232797:[ [ 0 1 ] [ 2 3 ] ]


```
