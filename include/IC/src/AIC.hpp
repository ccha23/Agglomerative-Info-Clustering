#ifndef _AIC_h_included_
#define _AIC_h_included_

namespace IC { 
    
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
		std::iota(idx.begin(), idx.end(), 0);
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
			std::iota(V.begin(), V.end(), 0);
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
		L(0, 0) = std::sqrt(1 + x.dot(x));	// maintain L*L.transpose() == e*e.transpose() + Q.transpose() * Q
										// for the first numCorral rows and columns.
		double old_first_order_opt,first_order_opt = INFINITY;
		size_t stuck_count;
		while (1) {
			// Step 1 (Major cycle: move towards origin along the direction of x)
			Eigen::VectorXd q = edmonds_greedy(f, x);
			old_first_order_opt = first_order_opt;
			first_order_opt = std::abs(x.dot(q) - x.dot(x));
			//log<LOG_INFO>(L"[min_norm_base] First order optimality : %1% ") % first_order_opt;
			// stopping criteria
			if (first_order_opt < fn_tol) { 
				return x; // x separates the base polytope from the origin and is therefore the minimum norm point.
			}
			if (old_first_order_opt <= first_order_opt) { // stuck
				stuck_count++;
// 				std::cout << "stuck : " << stuck_count << " " << first_order_opt << std::endl;
				if (stuck_count > 100) { // avoid getting stuck for too long
					std::cerr << "Unstuck at gap to optimality: " << first_order_opt << std::endl;
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
    
    /**
	Construct the exact info-clustering solution using the minimum normal base algorithm.
	*/
	class AIC : public HC {
		// for the submodular function f on {0,...,n-1}, and index j in {0,...,n-1}, define
		// f_(B) = f(B\cup {j})
		// for B subset of {0,...,j-1}
		struct SF__ : SF {
			SF &f;
			size_t j;
			SF__(SF &f, size_t j) : f(f), j(j) {
			}
			double operator() (const std::vector<size_t> &B) const {
				size_t n = f.size();
				size_t m = B.size();
				std::vector<size_t> B_(m + 1);
				for (size_t i = 0; i < m; i++) {
					if (B[i] < j) {
						B_[i] = B[i];
					}
				}
				B_[m] = j;
				return f(B_);
			}
			size_t size() const {
				return j;
			}
		};
		// for partition P={C_1,C_2,...,C_n} of {0,...,n'-1} and submodular function f on {0,..,n'-1}, define the submodular function 
		// f_(B):= f(\cup_{i\in B} C_i) - sum_{i\in B} f(C_i)
		// for B subset {0,...,n-1}
		struct SF_ : SF {
			SF &f;
			std::vector<std::vector<size_t> > P;
			std::vector<double> fi;

			SF_(SF &f, std::vector<std::vector<size_t> > P) : f(f), P(P) {
				size_t n = P.size();
				fi.resize(n);
				for (size_t i = 0; i < n; i++) {
					fi[i] = f(P[i]);
				}
			}
			SF_(SF &f) : f(f) {
				size_t n = f.size();
				P.resize(n);
				fi.resize(n);
				for (size_t i = 0; i<n; i++) {
					P[i] = std::vector<size_t>{ i };
					fi[i] = f(P[i]);
				}
			}
			double operator() (const std::vector<size_t> &B) const {
				size_t n = P.size();
				size_t m = B.size();
				std::vector<size_t> B_;
				double ss = 0;
				for (size_t i = 0; i < m; i++) {
					size_t i_ = B[i];
					//if (i_<0 || i_ >= n) throw std::runtime_error("B must be a subset of the ground set.");
					B_.insert(B_.end(), P[i_].begin(), P[i_].end());
					ss += fi[i_];
				}
				return f(B_) - ss;
			}
			size_t size() const {
				return P.size();
			}
		};
	private:
		SF& f;
	public:
		/**
		Construct PSP for the submodular function f in an agglomerative fashion.
		@param f the submodular function.
		*/
		AIC(SF &f_) : f(f_), HC(f_.size()) {}

		/**
		Agglomerate the cluster. Return true if clusters are agglomerated, and false otherwise.
		@param fn_tol Function Tolerance
		@param eps Precision
		*/
		bool agglomerate(double fn_tol, double eps) {
			double gamma = -INFINITY;
			std::vector<std::vector<size_t>> P=this->getPartition(gamma);
			if (P.size() <= 1) return false;
			SF_ f_(f,P);
			size_t k = P.size();
			//std::cout << "P : " << P << endl;
			std::vector<Eigen::VectorXd> x(k - 1);
			for (size_t j = 1; j < k; j++) {
				SF__ f__(f_, j);
				//std::cout << "j : " << j << endl;
				x[j-1] = min_norm_base(f__,fn_tol,eps);
				//std::cout << x[j - 1];
				gamma = std::max(gamma, -(x[j-1].minCoeff()));
			}
			//std::cout << "gamma : " << gamma << endl;
			for (size_t j = 1; j < k; j++) {
				for (size_t i = 0; i < j; i++) {
					//std::cout << x[j - 1](i) << endl;
					if (x[j - 1](i) <= -gamma + fn_tol) // instead of eps, = is required to handle -INFINITY
						//std::cout << "fusing " << P[j][0] << " " << P[i][0] << endl;
						merge(P[j][0], P[i][0], gamma);
				}
			}
			return true;
		}

		bool agglomerate() {
			return agglomerate(1E-6, 1E-15);
		}
	};

}

#endif