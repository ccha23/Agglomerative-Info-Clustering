#ifndef _SF_h_included_
#define _SF_h_included_

#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

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
	virtual double operator() (const vector<size_t> &B) const = 0;
	/*
	@return The size of the ground set.
	*/
	virtual size_t size() const=0;
};

class GaussianEntropy : public SF {
private:
	MatrixXd Sigma; // covariance matrix

public:
	/*
	Construct a submodular function as the entropy function of a gaussian random vector with specified covariance matrix.
	@param S Covariance matrix.
	*/
	GaussianEntropy(MatrixXd S);
	/*
	Calculate the entropy of a gaussian vector with a given covariance matrix.
	@param S Covariance matrix.
	@return differential entropy of the gaussian vector with covariance matrix S.
	*/
	static double h(const MatrixXd &S);

	/*
	Calculate the entropy of a gaussian subvector.
	@param B subvector of elements from the ground set.
	@return Entropy of the gaussian subvector indexed by elements in B.
	*/
	double operator() (const vector<size_t> &B) const;

	size_t size() const;
};

/*
Find the point in the base polytope of the normalized version of a submodular function
that minimizes the weighted sum of coordinates.
@param f The submodular function
@param w The weight vector with length equal to that of the ground set of f.
@return The vector x in B(f-f(emptySet)) such w'x is minimized.
*/
VectorXd edmonds_greedy(const SF &f, const VectorXd &w);

/*
Compute the minimum norm base in the base polytope of the normalized version of a submodular function with
the function tolerance 1E-10, optimality tolerance 1E-10 and precision epsilon 1E-15.
@param f The submodular function
@return The minimum norm point of the base B(f-f(emptySet)).
*/
VectorXd min_norm_base(const SF &f);

/*
Compute the minimum norm base in the base polytope of the normalized version of a submodular function
with specified tolerance and precision. N.b., set function tolerance to 1E-6 or above due to the numerical issue with Eigen.
@param f The submodular function
@param fn_tol Function Tolerance
@param eps Precision
@return The minimum norm point of the base B(f-f(emptySet)).
*/
VectorXd min_norm_base(const SF &f, double fn_tol, double eps);

#endif