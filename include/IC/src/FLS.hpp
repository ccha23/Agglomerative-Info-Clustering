#ifndef _FLS_h_included_
#define _FLS_h_included_

#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include <pari/pari.h>

// using namespace std;
// using namespace Eigen;

namespace IC { 
    long grank(Eigen::MatrixXd A, int p, int k)
    {
      GEN p1;	  /* vec */
      size_t n_cols = A.cols();
      size_t n_rows = A.rows();
      p1 = cgetg(n_cols+1, t_MAT); // p1 = cgetg(4, t_MAT);
      for(size_t i = 0; i < n_cols; i++){ 
         gel(p1, i+1) = cgetg(n_rows+1, t_COL);
      }


      for(size_t i = 0; i < n_rows; i++){ 
          for(size_t j = 0; j < n_cols; j++){ 
              gcoeff(p1, i+1, j+1) = stoi(A(i,j));
          }
      }


      GEN M; 
    //  GEN T = ffinit(stoi(p), k, -1);
    //  GEN x;
    //  x = ffgen(T, 1);
    //  M = gmul(p1, gpowgs(T, 3));

      GEN T = ffinit(stoi(p), k, -1);
      T = ffgen(T, -1);
      //pari_printf("primitive element = %Ps\n ", T);
      M = gmul(p1, gpowgs(T, 0));

    //  // print Matrix 
    //  for(int i=1; i<=n_rows; i++){
    //	  printf("\n");
    //	  for(int j=1; j<=n_cols; j++)
    //		  pari_printf("%Ps, \t ", gcoeff(M,i,j));
    //  }; printf("\n");

      //return rank(M);
        return 0;
    }
    
    class FiniteLinearEntropy : public SF{
	private:
		Eigen::MatrixXd LinearMapping; // LinearMapping matrix over the ring Z/mZ, i.e., entries from {0, ..., m-1}

	public:
		int pField;
		int kField;
		/*
		Construct a submodular function as the entropy function of a hypergraphical source specified by the incidence matrix.
		@param M LinearMapping is a |V| x |S| matrix s.t. Z_V = M * X_S, X_i i\in
S are independent and uniformly distributed over Z/qZ
		*/
		FiniteLinearEntropy(Eigen::MatrixXd M, int p, int k) {
			LinearMapping = M;
			pField = p;
			kField = k;
		}
		/*
		Calculate the entropy of a finite linear source with a given linear mapping matrix.
		@param M LinearMapping matrix.
		@return The entropy of the finite linear source with linear mapping matrix M.
		*/
		static double h(const Eigen::MatrixXd &M, int p, int k) {
			double ent = k*log2(p) * grank(M,p,k); 
			return ent;
		}

		/*
		Calculate the entropy of a hypergraphical source.
		@param B subvector of elements from the ground set.
		@return Entropy of the component sources indexed by elements in B.
		*/
		double operator() (const std::vector<size_t> &B) const {
			size_t n = B.size();
			Eigen::MatrixXd M(n, LinearMapping.cols());
			for (size_t i = 0; i < n; i++) {
				M.row(i) = LinearMapping.row(B[i]);
			}
			return h(M,pField,kField);
		}

		size_t size() const {
			return LinearMapping.rows();
		}
	};

}

#endif