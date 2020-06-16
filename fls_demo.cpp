#include <IC/AIC> // for exact info-clustering solution
#include <IC/ChowLiu> // for approximate solution
#include <IC/FLS> // for finite linear source model

int main() {
    // define the linear mappings
    int pField = 2;
    int kField = 1;

    size_t m = 3;
    size_t n0 = 1, n1 = 1, n2 = 1, n3 = 0;
    Eigen::MatrixXd M0(n0,m);
    Eigen::MatrixXd M1(n1,m);
    Eigen::MatrixXd M2(n2,m);
    Eigen::MatrixXd M3(n3,m);
    M0 << 1,0,0;
    M1 << 0,1,0;  
    M2 << 1,1,0;

    std::cout << "Linear mappings: " 
      << std::endl << "M0= \n" << M0 << std::endl
      << std::endl << "M1= \n" << M1 << std::endl
      << std::endl << "M2= \n" << M2 << std::endl
      << std::endl << "M3= \n" << M3 << std::endl;

    pari_init(1000000,2);

    // compute pairwise mutual information
    std::vector<Eigen::MatrixXd> LinearMappings = {M0,M1,M2,M3};
    IC::FiniteLinearEntropy h(LinearMappings,pField,kField);
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
    
    pari_close();  
}