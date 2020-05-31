#ifndef _SF_h_included_
#define _SF_h_included_

#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

// using namespace std;
// using namespace Eigen;

namespace IC { 

	/*
	Base functor for submodular function. Ground set V is {0 ... n-1} for some non-negative number n.
	*/
	class SF {
	public:
		/*
		The value of the submodular function evaluated on a set.
		@param B A subvector of elements of V.
		@return The value of the submodular function at B.
		*/
		virtual double operator() (const std::vector<size_t> &B) const = 0;
		/*
		@return The size of the ground set.
		*/
		virtual size_t size() const=0;
	};

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
				std::cout << "S submatrix " << std::endl << S << std::endl;
				throw std::runtime_error("S must be an SPD matrix.");
			}
			size_t k = S.rows();
			double det = c*k;
			Eigen::MatrixXd L = lltOfS.matrixL();
			for (size_t i = 0; i < k; i++) {
				det += log(L(i, i));
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

	/*
	Find the point in the base polytope of the normalized version of a submodular function
	that minimizes the weighted sum of coordinates.
	@param f The submodular function
	@param w The weight vector with length equal to that of the ground set of f.
	@return The vector x in B(f-f(emptySet)) such w'x is minimized.
	*/
	Eigen::VectorXd edmonds_greedy(const SF &f, const Eigen::VectorXd &w) {
		size_t n = f.size();
		if (n != w.size())
			throw std::runtime_error("w must have the same size as the ground set of f.");
		std::vector<size_t> idx(n);
		Eigen::VectorXd x(n);
		iota(idx.begin(), idx.end(), 0);
		sort(idx.begin(), idx.end(), [&](size_t i1, size_t i2) {return w[i1] < w[i2]; });
		std::vector<size_t> B;
		double F = f(B);
		for (auto i : idx) {
			x[i] = -F;
			B.push_back(i);
			F = f(B);
			x[i] += F;
		}
		return x;
	}

	/*
	Compute the minimum norm base in the base polytope of the normalized version of a submodular function
	with specified tolerance and precision. N.b., set function tolerance to 1E-6 or above due to the numerical issue with Eigen.
	@param f The submodular function
	@param fn_tol Function Tolerance
	@param eps Precision
	@return The minimum norm point of the base B(f-f(emptySet)).
	*/
	Eigen::VectorXd min_norm_base(const SF &f, double fn_tol, double eps) {
		size_t n = f.size();
		Eigen::VectorXd x;
		if (n < 1)
			return x;
		x = Eigen::VectorXd::Ones(n);
		if (n == 1) {
			std::vector<size_t> V(n);
			iota(V.begin(), V.end(), 0);
			return x*f(V);
		}
		// Step 0 (Initialize a trivial corral Q and x=Nr Q)
		x = edmonds_greedy(f, x);
		//std::cout << f(std::vector<size_t> {0, 1}) - f(std::vector<size_t> {0}) << std::endl;
		//std::cout << x.transpose() << std::endl;
		Eigen::MatrixXd Q(n, n + 1);
		Q.col(0) = x;
		size_t numCorral = 1; // Corral: First numCorral columns of Q

		Eigen::VectorXd w = Eigen::VectorXd::Zero(n + 1);
		w(0) = 1; // maintain x == Q * w

				// initialization for adaptive computation of the projection to affine hull of the corral
		Eigen::VectorXd e = Eigen::VectorXd::Ones(n + 1);
		Eigen::MatrixXd L = Eigen::MatrixXd::Zero(n + 1, n + 1);
		L(0, 0) = sqrt(1 + x.dot(x));	// maintain L*L.transpose() == e*e.transpose() + Q.transpose() * Q
										// for the first numCorral rows and columns.
		double old_first_order_opt,first_order_opt = INFINITY;
		size_t stuck_count;
		while (1) {
			// Step 1 (Major cycle: move towards origin along the direction of x)
			Eigen::VectorXd q = edmonds_greedy(f, x);
			old_first_order_opt = first_order_opt;
			first_order_opt = abs(x.dot(q) - x.dot(x));
			//log<LOG_INFO>(L"[min_norm_base] First order optimality : %1% ") % first_order_opt;
			// stopping criteria
			if (first_order_opt < fn_tol) { 
				return x; // x separates the base polytope from the origin and is therefore the minimum norm point.
			}
			if (old_first_order_opt <= first_order_opt) { // stuck
				stuck_count++;
				std::cout << "stuck : " << stuck_count << " " << first_order_opt << std::endl;
				if (stuck_count > 100) { // avoid getting stuck for too long
					std::cout << "Unstuck at gap to optimality: " << first_order_opt << std::endl;
					return x;
				}
			}
			else {
				stuck_count = 0;
			}
			{
				// Step 2 (Major cycle: project to affine hull)
				auto L_ = L.topLeftCorner(numCorral, numCorral).triangularView<Eigen::Lower>();
				auto Q_ = Q.topLeftCorner(n, numCorral);
				auto e_ = e.topRows(numCorral);
				Eigen::VectorXd r = L_.solve(e_ + Q_.transpose() * q);
				L.block(numCorral, 0, 1, numCorral) = r.transpose();
				L(numCorral, numCorral) = 1 + q.dot(q) - r.dot(r);
				//std::cout << "Optimality test:" << std::endl;
				//std::cout << "x : " << x.transpose() << std::endl;
				//std::cout << "q : " << q.transpose() << std::endl;
				//std::cout << "Q : " << std::endl << Q_ << std::endl;
				if (1 + q.dot(q) - r.dot(r) < 0) {
					//std::cout << "Update L to include q:" << std::endl;
					//std::cout << "L: " << std::endl << (Eigen::MatrixXd) L_ << std::endl;
					//std::cout << "r: " << r.transpose() << std::endl;
					throw std::runtime_error("Try to set a higher functional tolerance.");
				}
				L(numCorral, numCorral) = sqrt(1 + q.dot(q) - r.dot(r));
				Q.col(numCorral++) = q;
			}
			int minor_cycle_count = 0; // count number of iterations of minor cycles in iteration of the major cycle 
			while (1) {
				minor_cycle_count = minor_cycle_count + 1;
				auto L_ = L.topLeftCorner(numCorral, numCorral).triangularView<Eigen::Lower>();
				auto Q_ = Q.topLeftCorner(n, numCorral);
				auto e_ = e.topRows(numCorral);
				Eigen::VectorXd v = L_.transpose().solve(L_.solve(e_));
				v = (e_ * e_.transpose() + Q_.transpose() * Q_).llt().solve(e_); // cc
				v = v / e_.dot(v);
				//std::cout << "Affine projection: origin onto aff(Q) " << std::endl;
				//std::cout << "L : " << std::endl << (Eigen::MatrixXd)L_ << std::endl;
				//std::cout << "Q : " << std::endl << Q_ << std::endl;
				//std::cout << "v : " << v.transpose() << std::endl;
				bool minor_cycle = false;
				{
					size_t i = 0;
					while (!minor_cycle && i < numCorral)
						minor_cycle = (v(i++) <= eps);
				}
				auto w_ = w.topRows(numCorral);
				if (!minor_cycle) {
					// continue with Main Cycle
					w_ = v;
					x = Q_*w_;
					//log<LOG_INFO>(L"[min_norm_base] # minor cycles : %1% ") % minor_cycle_count;
					break;
				}
				// Step 3: Minor cycle
				double theta = 1;
				for (size_t i = 0; i < numCorral; i++) {
					if (w_(i) - v(i)>= eps) {
						theta = std::min(theta, w_(i) / (w_(i) - v(i)));
					}
				}
				w_ = (1 - theta) * w_ + theta *v;
				size_t i_ = 0;
				std::vector<bool> toDelete(numCorral);
				for (size_t i = 0; i < numCorral; i++) {
					if ((toDelete[i] = (w_(i) <= eps))) {
						for (size_t j = i; j < numCorral - 1; j++) {
							// Maintain L triangular after deletion
							double a = L(j + 1, j), b = L(j + 1, j + 1);
							double c = sqrt(a*a + b*b);
							Eigen::Matrix2d T;
							T << a / c, -b / c, b / c, a / c;
							L.block(j + 1, j, numCorral - j, 2) *= T;
						}
					}
				}
				// perform deletion
				for (size_t i = 0; i < numCorral; i++) {
					if (!toDelete[i]) {
						if (i_<i) {
							w_(i_) = w_(i);
							Q.col(i_) = Q.col(i);
							L.row(i_) = L.row(i);
						} // else, i.e., i_=i : no need to move
						i_++;
					} // else: delete ith corral
				}
				numCorral = i_;
				x = Q.topLeftCorner(n, numCorral)*w.topRows(numCorral);
			}
		}
		return x;
	}
		
	/*
	Compute the minimum norm base in the base polytope of the normalized version of a submodular function with
	the function tolerance 1E-10, optimality tolerance 1E-10 and precision epsilon 1E-15.
	@param f The submodular function
	@return The minimum norm point of the base B(f-f(emptySet)).
	*/
	Eigen::VectorXd min_norm_base(const SF &f) {
		return min_norm_base(f, 1E-10, 1E-15);
	}

}
#endif