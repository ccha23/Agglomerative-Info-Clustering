#ifndef _hypergraph_h_included_
#define _hypergraph_h_included_

namespace IC { 

	class HypergraphEntropy : public SF {
	private:
		Eigen::MatrixXd Incidence; // Incidence matrix

	public:
		/*
		Construct a submodular function as the entropy function of a hypergraphical source specified by the incidence matrix.
		@param M Incidence matrix.
		*/
		HypergraphEntropy(Eigen::MatrixXd M) {
			Incidence = M;
		}
		/*
		Calculate the entropy of a hypergraphical source with a given incidence matrix.
		@param M Incidence matrix.
		@return The entropy of the hypergraphical source with incidence matrix M.
		*/
		static double h(const Eigen::MatrixXd &M) {
			size_t k = M.cols();
			double ent = 0;
			for (size_t i = 0; i < k; i++) {
				ent += M.col(i).maxCoeff();
			}
			return ent;
		}

		/*
		Calculate the entropy of a hypergraphical source.
		@param B subvector of elements from the ground set.
		@return Entropy of the component sources indexed by elements in B.
		*/
		double operator() (const std::vector<size_t> &B) const {
			size_t n = B.size();
			Eigen::MatrixXd M(n, Incidence.cols());
			for (size_t i = 0; i < n; i++) {
				M.row(i) = Incidence.row(B[i]);
			}
			return h(M);
		}

		size_t size() const {
			return Incidence.rows();
		}
	};

}
#endif