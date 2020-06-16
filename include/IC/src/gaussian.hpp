#ifndef _gaussian_h_included_
#define _gaussian_h_included_

#define _USE_MATH_DEFINES

namespace IC { 

    	class GaussianEntropy : public SF {
	private:
		std::vector<Eigen::MatrixXd> M; // vector of linear mappings
        size_t m, n;

	public:
		/**
		Construct a submodular function as the entropy function of a gaussian source by a vector of linear mappings.
        
		@param LinearMappings vector of linear mappings in terms of matrices with the same number of columns.
		*/
		GaussianEntropy(std::vector<Eigen::MatrixXd> LinearMappings) {
			M = LinearMappings;
            m = M.size();
            n = (m==0)? 0: M[0].cols();             
            for (size_t i=0; i<m; i++) {
                if (M[i].cols() != n)
                    throw std::runtime_error("Linear mappings must have the same number of columns.");
            }
        }
		/**
		Calculate the entropy (in bits) of a gaussian source with a given linear mapping.
		@param LinearMapping Covariance matrix.
		@return differential entropy of the gaussian vector generated using the linear mapping LinearMapping.
		*/
		static double h(const Eigen::MatrixXd &LinearMapping) {
            size_t k = LinearMapping.rows();
            if (k==0 || LinearMapping.cols()==0)
                return 0;
            Eigen::MatrixXd S = LinearMapping * LinearMapping.transpose();
			const double c = log2(2 * M_PI*M_E) / 2;
			Eigen::LLT<Eigen::MatrixXd> lltOfS(S);
            if (lltOfS.info() == Eigen::NumericalIssue)
				throw std::runtime_error("Rows must be linearly independent.");
			double det = c*k;
			Eigen::MatrixXd L = lltOfS.matrixL();
			for (size_t i = 0; i < k; i++) {
				det += std::log2(L(i, i));
			}
			return det;
		}

		/**
		Calculate the entropy (in bits) of a subvector of elements in the ground set.
		@param B subvector of elements from the ground set.
		@return Entropy of elements in B.
		*/
		double operator() (const std::vector<size_t> &B) const {
			size_t k = B.size();
            std::vector<size_t> row_idx(k+1); // row_idx[k] is the number of rows
            for (size_t i=0; i<B.size(); i++) {
                row_idx[i] = row_idx[k];
                row_idx[k] += M[B[i]].rows();
            }
			Eigen::MatrixXd LinearMapping(row_idx[k], n);
			for (size_t i = 0; i < k; i++) {
                LinearMapping.block(row_idx[i],0,row_idx[i+1]-row_idx[i],n) = M[B[i]];
			}
			return h(LinearMapping);
		}

		size_t size() const {
			return m;
		}
	};
    
	class ScalarGaussianEntropy : public SF {
	private:
		Eigen::MatrixXd Sigma; // covariance matrix

	public:
		/**
		Construct a submodular function as the entropy function of a gaussian random vector with specified covariance matrix.
		@param S Covariance matrix.
		*/
		ScalarGaussianEntropy(Eigen::MatrixXd S) {
			Sigma = S;
			if (Sigma.rows() != Sigma.cols())
				throw std::runtime_error("S must be an SPD matrix.");
		}
		/**
		Calculate the entropy of a gaussian vector with a given covariance matrix.
		@param S Covariance matrix.
		@return differential entropy of the gaussian vector with covariance matrix S.
		*/
		static double h(const Eigen::MatrixXd &S) {
			const double c = log(2 * M_PI*M_E) / 2;
			Eigen::LLT<Eigen::MatrixXd> lltOfS(S);
			if (lltOfS.info() == Eigen::NumericalIssue)
			{
				std::cerr << "S submatrix " << std::endl << S << std::endl;
				throw std::runtime_error("S must be an SPD matrix.");
			}
			size_t k = S.rows();
			double det = c*k;
			Eigen::MatrixXd L = lltOfS.matrixL();
			for (size_t i = 0; i < k; i++) {
				det += std::log(L(i, i));
			}
			return det;
		}

		/**
		Calculate the entropy of a gaussian subvector.
		@param B subvector of elements from the ground set.
		@return Entropy of the gaussian subvector indexed by elements in B.
		*/
		double operator() (const std::vector<size_t> &B) const {
			size_t n = B.size();
			Eigen::MatrixXd S(n, n);
			for (size_t i = 0; i < n; i++) {
				for (size_t j = 0; j < n; j++) {
					S(i, j) = Sigma(B[i], B[j]);
				}
			}
			return h(S);
		}

		size_t size() const {
			return Sigma.rows();
		}
	};
}
#endif