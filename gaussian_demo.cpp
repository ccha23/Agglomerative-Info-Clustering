#include <IC/AIC> // for exact info-clustering solution
#include <IC/ChowLiu> // for approximate solution
#include <IC/gaussian> // for gaussian source model

int main() {
    // define the linear mappings
    Eigen::MatrixXd M0(1,4), M1(1,4), M2(1,4);
    M0 << 1,1,0,0;
    M1 << 1,0,1,0;
    M2 << 0,0,0,1;
    std::cout << "Linear mappings: " 
              << std::endl << "M0= \n" << M0 << std::endl
              << std::endl << "M1= \n" << M1 << std::endl
              << std::endl << "M2= \n" << M2 << std::endl;
    std::vector<Eigen::MatrixXd> LinearMappings = {M0,M1,M2};
    IC::GaussianEntropy h(LinearMappings);
    
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