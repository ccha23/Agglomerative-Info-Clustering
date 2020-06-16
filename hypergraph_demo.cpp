#include <IC/AIC> // for exact info-clustering solution
#include <IC/ChowLiu> // for approximate solution
#include <IC/hypergraph> // for hypergraphical source model

int main() {

    // define the incidence matrix
	Eigen::MatrixXd M(6, 4); // 6 vertices and 4 edges
	M << 1, 1, 0, 0,
		 1, 1, 0, 0,
		 1, 0, 0, 0,
		 0, 0, 1, 0,
		 0, 0, 1, 0,
		 0, 0, 0, 1;
    
    std::cout << "Incidence matrix of the (weighted) hypergraph: " 
              << std::endl << "M= \n" << M << std::endl;
    IC::HypergraphEntropy h(M);
    
    // compute pairwise mutual information
    std::cout << std::endl << "Pairwise mutual information: " << std::endl;
    for (size_t i=0; i<h.size(); i++)
        for (size_t j=i+1; j<h.size(); j++)
        std::cout << "I(Z" << i << "; Z" << j << ") = " 
                  << h.mutual_info(i,j) << std::endl;    
    
    // generate exact info-clustering solution
    std::cout << std::endl << "Exact solution: " << std::endl;
    IC::AIC psp(h);
    std::cout << psp.getCoarsestPartition() 
              << std::endl; 
    while (psp.agglomerate()) 
            std::cout << "agglomerates to " << psp.getCoarsestPartition() 
                  << " at critical value " << psp.getCriticalValues().back() 
                  << std::endl;
    
    // generate approximate info-clustering solution
    std::cout << std::endl << "Approximate solution: " << std::endl;
    IC::CL cl(h);
    std::vector<double> cl_gamma = cl.getCriticalValues();
    std::cout << "critical values : " << cl_gamma << std::endl;
    for (double gamma : cl.getCriticalValues()) {
        std::cout << "partition at threshold " << gamma 
                  << ":" << cl.getPartition(gamma) << std::endl;
    }   
}