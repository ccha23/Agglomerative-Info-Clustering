#ifndef _gaussian_h_included_
#define _gaussian_h_included_

#define _USE_MATH_DEFINES

namespace IC { 

	class GaussianEntropy : public SF {
	private:
		Eigen::MatrixXd Sigma; // covariance matrix

	public:
		/*
		Construct a submodular function as the entropy function of a gaussian random vector with specified covariance matrix.
		@param S Covariance matrix.
		*/
		GaussianEntropy(Eigen::MatrixXd S) {
			Sigma = S;
			if (Sigma.rows() != Sigma.cols())
				throw std::runtime_error("S must be an SPD matrix.");
		}
		/*
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

		/*
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