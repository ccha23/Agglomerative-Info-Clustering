#ifndef _AIC_h_included_
#define _AIC_h_included_

#include <iostream>
#include <vector>
#include <list>
#include <math.h>
#include <queue>
#include <numeric>
#include "SF.h"

using namespace std;

/**
HC is a data structure that records the solution of an agglomerative clustering algorithm. It provides the function
merge to construct the clustering solution, and the functions getPartition, getCriticalValues and similarity to retrieve
information about the clusters.
*/
class HC {
	// Use custom doubly linked list to have constant time deletion and search.
	struct LinkList {
		struct LinkNode {
			size_t key, prev, next;
			LinkNode(size_t k, size_t p, size_t n) {
				key = k;
				prev = p;
				next = n;
			}
		};
		size_t head = 0;
		vector<LinkNode*> lnode;
		LinkList(size_t size, size_t mode) {
			// assume size>0
			lnode.resize(size);
			switch (mode)
			{
			case 0: // isolated nodes
				for (size_t i = 0; i < size; i++) {
					lnode[i] = new LinkNode(i, i, i);
				}
				break;
			case 1: // linked nodes
				if (size > 1) {
					lnode[0] = new LinkNode(0, 0, 1);
					for (size_t i = 1; i < size-1; i++) {
						lnode[i] = new LinkNode(i, i - 1, i + 1);
					}
					lnode[size-1] = new LinkNode(size - 1, size - 2, size - 1);
				} else {
					lnode[0] = new LinkNode(0, 0, 0);
				}
				break;
			}

		}
		~LinkList() {
			for (size_t i = 1, ilen = lnode.size(); i<ilen; i++) {
				delete lnode[i];
			}
		}
		void remove(size_t i) {
			if (head != i) { // ignore if attempt to remove head
				size_t prev = lnode[i]->prev;
				size_t next = lnode[i]->next;
				if (prev != i) {
					lnode[prev]->next = (next != i? next : prev);
					lnode[i]->prev = i;
				}
				if (next != i) {
					lnode[next]->prev = (prev !=i ? prev : next);
					lnode[i]->next = i;
				}
			}
		}
		bool next() {
			if (lnode[head]->next == head) return false;
			head = lnode[head]->next;
			return true;
		}
	};
protected:
	vector<double> gamma;
	vector<size_t> parent;
	vector<vector<size_t> > children;
	vector<double> weight;
	vector<size_t> rank;

	friend std::ostream& operator<<(std::ostream&, const HC&);

	/**
	Find the current cluster label of a node. O(log(|V|)) time.
	@param i node label.
	@return The current cluster label of node i.
	*/
	size_t find(size_t i) const;

	/**
	Merge nodes i and j into the same cluster for all threshold values strictly smaller than gamma. 
	Complexity: O(log n).
	@param i node label.
	@param j node label.
	@param gamma similarity level.
	@return true if a merge is done, i.e., i and j were not connected before.
	*/
	bool merge(size_t i, size_t j, double gamma);

public:
	HC(size_t n);

	/**
	Return the partition P of V, the non-singleton elements of which are the clusters at threshold gamma.
	Complexity: O(n).
	@param gamma similarity threshold.
	@return The partition P that gives the clusters at threshold gamma.
	*/
	vector<vector<size_t> > getPartition(double gamma) const;

	/**
	Return the similarity of node i and j. 
	Complexity: O(log n).
	@param i node label.
	@param j node label.
	@return The similarity of node i and j.
	*/
	double similarity(size_t i, size_t j) const;

	/**
	Return the set of critical values. 
	Complexity: O(n).
	@param i node label.
	@param j node label.
	@return The vector of critical values in descending order.
	*/
	vector<double> getCriticalValues();
};

class AIC : public HC {
	// for the submodular function f on {0,...,n-1}, and index j in {0,...,n-1}, define
	// f_(B) = f(B\cup {j})
	// for B subset of {0,...,j-1}
	struct SF__ : SF {
		SF &f;
		size_t j;
		SF__(SF &f, size_t j) : f(f), j(j) {
		}
		double operator() (const vector<size_t> &B) const {
			size_t n = f.size();
			size_t m = B.size();
			vector<size_t> B_(m + 1);
			for (size_t i = 0; i < m; i++) {
				if (B[i] < j && B[i] >= 0) {
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
		vector<vector<size_t> > P;
		vector<double> fi;

		SF_(SF &f, vector<vector<size_t> > P) : f(f), P(P) {
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
				P[i] = vector<size_t>{ i };
				fi[i] = f(P[i]);
			}
		}
		double operator() (const vector<size_t> &B) const {
			size_t n = P.size();
			size_t m = B.size();
			vector<size_t> B_;
			double ss = 0;
			for (size_t i = 0; i < m; i++) {
				size_t i_ = B[i];
				//if (i_<0 || i_ >= n) throw runtime_error("B must be a subset of the ground set.");
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
	AIC(SF &f);

	/**
	Agglomerate the cluster. Return true if clusters are agglomerated, and false otherwise.
	@param fn_tol Function Tolerance
	@param eps Precision
	*/
	bool agglomerate(double fn_tol, double eps);
	bool agglomerate();
};

/**
Construct a info-clustering solution by Chow-Liu tree approximation.
*/
class CL : public HC {
public:
	/**
	Construct the clusters for a graph given a list of non-negative weighted edges.
	@param size The size |V| of the graph.
	@param first_node Vector of first incident nodes of the edges.
	@param first_node Vector of second incident nodes of the edges.
	@param gamma Vector of the weights of the edges.
	@return The PSP of the graph.
	*/
	CL(size_t size, vector<size_t> first_node, vector<size_t> second_node, vector<double> gamma);
};

/**
Represent vectors by strings.
*/
template <typename T>

std::ostream& operator<<(std::ostream &strm, const vector<T> &v) {
	strm << "[ ";
	for (size_t i = 0, ilen = v.size(); i < ilen; i++)
		strm << v[i] << " ";
	return strm << "]";
}

/**
Represent pairs by strings.
*/
template <typename T1,typename T2>
std::ostream& operator<<(std::ostream &strm, const pair<T1,T2> &p) {
	return strm << "(" << p.first << ", " << p.second << ")";
}

/**
Sort the indices of a vector in descending order of the corresponding elements.
@param v vector
@return The vector of sorted indices.
*/
template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

	// initialize original index locations
	vector<size_t> idx(v.size());
	iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });

	return idx;
}

#endif