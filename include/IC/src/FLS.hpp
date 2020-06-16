#ifndef _FLS_h_included_
#define _FLS_h_included_

#include <pari/pari.h>

namespace IC { 
    long grank(Eigen::MatrixXd A, int p, int k)
    {
      GEN p1;	  /* vec */
      size_t n_cols = A.cols();
      size_t n_rows = A.rows();
      p1 = cgetg(n_cols+1, t_MAT);
      for(size_t i = 0; i < n_cols; i++){ 
         gel(p1, i+1) = cgetg(n_rows+1, t_COL);
      }


      for(size_t i = 0; i < n_rows; i++){ 
          for(size_t j = 0; j < n_cols; j++){ 
              gcoeff(p1, i+1, j+1) = stoi(A(i,j));
          }
      }


      GEN M; 
      GEN T = ffinit(stoi(p), k, -1);
      T = ffgen(T, -1);
      M = gmul(p1, gpowgs(T, 0));

//		// check if x^{p^k-1} = 1 and print the FF nonzero elements
//		pari_printf("sanity check: x^{p^k-1} = %Ps\n", gpowgs(T, pow(p,k)-1));
//		printf("Field's nonzero elements: ");
//		for(int i=0; i<pow(p,k)-1; i++)
//			pari_printf("%Ps,\t ", gpowgs(T, i)); printf("\n");

//		// print Matrix 
//		for(int i=1; i<=n_rows; i++){
//		printf("\n");
//		for(int j=1; j<=n_cols; j++)
//		  pari_printf("%Ps, \t ", gcoeff(M,i,j));
//		}; printf("\n");

//		std::cout << "r = " << rank(M) << std::endl;

      return rank(M);
    }
    
    class FiniteLinearEntropy : public SF{
	private:
		std::vector<Eigen::MatrixXd> LinearMappings; // LinearMapping matrix over the ring Z/mZ, i.e., entries from {0, ..., m-1}

	public:
		int pField;
		int kField;
		/**
		Construct a submodular function as the entropy function of a finite linear source specified by a vector of linear mappings.
		@param vecM a vector of LinearMapping M0, M1, ... where Mi is a matrix s.t. the features of i are given by z_i = M_i * x for x independent and uniformly distributed over GF(p,k).
		*/
		FiniteLinearEntropy(std::vector<Eigen::MatrixXd> M, int p, int k) {
			LinearMappings = M;
			pField = p;
			kField = k;
		}
		/**
		Calculate the entropy of a finite linear source with a given linear mapping matrix.
		@param M LinearMapping matrix.
		@return The entropy of the finite linear source with linear mapping matrix M.
		*/
		static double h(const Eigen::MatrixXd &M, int p, int k) {
			double ent = k*log2(p) * grank(M,p,k); 
			return ent;
		}

		/**
		Calculate the entropy of a hypergraphical source.
		@param B subvector of elements from the ground set.
		@return Entropy of the component sources indexed by elements in B.
		*/
		double operator() (const std::vector<size_t> &B) const {
			size_t n_rows = 0;
			for (size_t i=0; i<B.size(); i++) {
				n_rows += LinearMappings[B[i]].rows();
			}
			Eigen::MatrixXd M(n_rows, LinearMappings[0].cols());
			size_t idx_rows = 0;
			for (size_t i=0; i<B.size(); i++)
				for (size_t j = 0; j<LinearMappings[B[i]].rows(); j++){
					M.row(idx_rows) = LinearMappings[B[i]].row(j);
					idx_rows++;
			}
			//std::cout << "Concatenated Mappings = \n" << M << std::endl;
			return h(M,pField,kField);
		}

		size_t size() const {
			return LinearMappings.size();
		}
	};

}

#endif
