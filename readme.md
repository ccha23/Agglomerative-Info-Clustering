[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ccha23/Agglomerative-Info-Clustering/master?urlpath=lab/tree/demo.ipynb) 

The above binder link launches a notebook that demonstrates the agglomerative info-clustering algorithm (AIC) implemented in C++. It is run using the [xeus-cling](https://xeus-cling.readthedocs.io/en/latest/) C++11 jupyter kernel.

- [Jupyter notebook preview](http://nbviewer.jupyter.org/github/ccha23/Agglomerative-Info-Clustering/blob/master/demo.ipynb)
- [API documentation](https://ccha23.github.io/Agglomerative-Info-Clustering)
- [Paper](https://www.overleaf.com/read/zkkyyxfvbmsc)

## Installation


Simply download the `include` folder and add it to the include path. 
- The source code requires [Eigen C++ library](http://eigen.tuxfamily.org), which is included in `include` folder. 
- To use the finite linear source model, you would need to install [`libpari-dev`](https://pari.math.u-bordeaux.fr). (See the next section.)

## Examples

Some examples from the jupyter notebook are also included as `.cpp` files as follows. See the jupyter notebook for more detailed explanations:

- gaussian_demo.cpp performs the exact and approximate clustering algorithm for a jointly gaussian source model. Run using 
```
make gaussian_demo; ./gaussian_demo.out
```

- hypergraph_demo.cpp is an example for hypergraphical source model. Run using
```
make hypergraph_demo; ./hypergraph_demo.out
```

- fls_demo.cpp is an example for clustering a finite linear source model. It requires the C library pari, which can be installed on ubuntu using
```
apt-get install libpari-dev
``` 
Run the using
```
make fls_demo; ./fls_demo.out
```