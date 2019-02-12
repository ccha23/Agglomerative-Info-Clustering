#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "SF.h"
//#include "Log.h"

using namespace std;
using namespace Eigen;

size_t GaussianEntropy::size() const {
	return Sigma.rows();
}

GaussianEntropy::GaussianEntropy(MatrixXd S) {
	Sigma = S;
	if (Sigma.rows() != Sigma.cols())
		throw std::runtime_error("S must be an SPD matrix.");
}

double GaussianEntropy::h(const MatrixXd &S) {
	//const double c = log(2 * M_PI*M_E) / 2;
	//size_t k = S.rows();
	//if (k < 1) return 0; // determinant of 0x0 matrix is 1. LDLT fails to handle this.
	//LDLT<MatrixXd> ldltOfS(S);
	//double det = c*k;
	//if (!ldltOfS.isPositive())
	//{
	//	VectorXd D = ldltOfS.vectorD();
	//	for (size_t i = 0; i < k; i++) {
	//		if (D(i) < -1E-10) { // handle numerical issues
	//			throw std::runtime_error("S must be an SPD matrix.");
	//		}
	//		else if (D(i) < 0) {
	//			D(i) = 0;
	//		}
	//		det += 2 * log(D(i));
	//	}
	//}
	//return det;
	const double c = log(2 * M_PI*M_E) / 2;
	LLT<MatrixXd> lltOfS(S);
	if (lltOfS.info() == Eigen::NumericalIssue)
	{
		cout << "S submatrix " << endl << S << endl;
		throw std::runtime_error("S must be an SPD matrix.");
	}
	size_t k = S.rows();
	double det = c*k;
	MatrixXd L = lltOfS.matrixL();
	for (size_t i = 0; i < k; i++) {
		det += log(L(i, i));
	}
	return det;
}

double GaussianEntropy::operator() (const vector<size_t> &B) const {
	size_t n = B.size();
	MatrixXd S(n, n);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			S(i, j) = Sigma(B[i], B[j]);
		}
	}
	return h(S);
}

VectorXd edmonds_greedy(const SF &f, const VectorXd &w) {
	size_t n = f.size();
	if (n != w.size())
		throw runtime_error("w must have the same size as the ground set of f.");
	vector<size_t> idx(n);
	VectorXd x(n);
	iota(idx.begin(), idx.end(), 0);
	sort(idx.begin(), idx.end(), [&](size_t i1, size_t i2) {return w[i1] < w[i2]; });
	vector<size_t> B;
	double F = f(B);
	for (auto i : idx) {
		x[i] = -F;
		B.push_back(i);
		F = f(B);
		x[i] += F;
	}
	return x;
}

VectorXd min_norm_base(const SF &f) {
	return min_norm_base(f, 1E-10, 1E-15);
}

VectorXd min_norm_base(const SF &f, double fn_tol, double eps) {
	size_t n = f.size();
	VectorXd x;
	if (n < 1)
		return x;
	x = VectorXd::Ones(n);
	if (n == 1) {
		vector<size_t> V(n);
		iota(V.begin(), V.end(), 0);
		return x*f(V);
	}
	// Step 0 (Initialize a trivial corral Q and x=Nr Q)
	x = edmonds_greedy(f, x);
	//cout << f(vector<size_t> {0, 1}) - f(vector<size_t> {0}) << endl;
	//cout << x.transpose() << endl;
	MatrixXd Q(n, n + 1);
	Q.col(0) = x;
	size_t numCorral = 1; // Corral: First numCorral columns of Q

	VectorXd w = VectorXd::Zero(n + 1);
	w(0) = 1; // maintain x == Q * w

			  // initialization for adaptive computation of the projection to affine hull of the corral
	VectorXd e = VectorXd::Ones(n + 1);
	MatrixXd L = MatrixXd::Zero(n + 1, n + 1);
	L(0, 0) = sqrt(1 + x.dot(x));	// maintain L*L.transpose() == e*e.transpose() + Q.transpose() * Q
									// for the first numCorral rows and columns.
	double old_first_order_opt,first_order_opt = INFINITY;
	size_t stuck_count;
	while (1) {
		// Step 1 (Major cycle: move towards origin along the direction of x)
		VectorXd q = edmonds_greedy(f, x);
		old_first_order_opt = first_order_opt;
		first_order_opt = abs(x.dot(q) - x.dot(x));
		//log<LOG_INFO>(L"[min_norm_base] First order optimality : %1% ") % first_order_opt;
		// stopping criteria
		if (first_order_opt < fn_tol) { 
			return x; // x separates the base polytope from the origin and is therefore the minimum norm point.
		}
		if (old_first_order_opt <= first_order_opt) { // stuck
			stuck_count++;
			cout << "stuck : " << stuck_count << " " << first_order_opt << endl;
			if (stuck_count > 100) { // avoid getting stuck for too long
				cout << "Unstuck at gap to optimality: " << first_order_opt << endl;
				return x;
			}
		}
		else {
			stuck_count = 0;
		}
		{
			// Step 2 (Major cycle: project to affine hull)
			auto L_ = L.topLeftCorner(numCorral, numCorral).triangularView<Lower>();
			auto Q_ = Q.topLeftCorner(n, numCorral);
			auto e_ = e.topRows(numCorral);
			VectorXd r = L_.solve(e_ + Q_.transpose() * q);
			L.block(numCorral, 0, 1, numCorral) = r.transpose();
			L(numCorral, numCorral) = 1 + q.dot(q) - r.dot(r);
			//cout << "Optimality test:" << endl;
			//cout << "x : " << x.transpose() << endl;
			//cout << "q : " << q.transpose() << endl;
			//cout << "Q : " << endl << Q_ << endl;
			if (1 + q.dot(q) - r.dot(r) < 0) {
				//cout << "Update L to include q:" << endl;
				//cout << "L: " << endl << (MatrixXd) L_ << endl;
				//cout << "r: " << r.transpose() << endl;
				throw runtime_error("Try to set a higher functional tolerance.");
			}
			L(numCorral, numCorral) = sqrt(1 + q.dot(q) - r.dot(r));
			Q.col(numCorral++) = q;
		}
		int minor_cycle_count = 0; // count number of iterations of minor cycles in iteration of the major cycle 
		while (1) {
			minor_cycle_count = minor_cycle_count + 1;
			auto L_ = L.topLeftCorner(numCorral, numCorral).triangularView<Lower>();
			auto Q_ = Q.topLeftCorner(n, numCorral);
			auto e_ = e.topRows(numCorral);
			VectorXd v = L_.transpose().solve(L_.solve(e_));
			v = (e_ * e_.transpose() + Q_.transpose() * Q_).llt().solve(e_); // cc
			v = v / e_.dot(v);
			//cout << "Affine projection: origin onto aff(Q) " << endl;
			//cout << "L : " << endl << (MatrixXd)L_ << endl;
			//cout << "Q : " << endl << Q_ << endl;
			//cout << "v : " << v.transpose() << endl;
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
					theta = min(theta, w_(i) / (w_(i) - v(i)));
				}
			}
			w_ = (1 - theta) * w_ + theta *v;
			size_t i_ = 0;
			vector<bool> toDelete(numCorral);
			for (size_t i = 0; i < numCorral; i++) {
				if ((toDelete[i] = (w_(i) <= eps))) {
					for (size_t j = i; j < numCorral - 1; j++) {
						// Maintain L triangular after deletion
						double a = L(j + 1, j), b = L(j + 1, j + 1);
						double c = sqrt(a*a + b*b);
						Matrix2d T;
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

