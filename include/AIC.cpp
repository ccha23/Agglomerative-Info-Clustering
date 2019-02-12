#include "AIC.h"
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <math.h>
#include <queue>
#include <algorithm>
#include "SF.h"

using namespace std;

// Constructors

AIC::AIC(SF &f_) : f(f_), HC::HC(f_.size()) {
}

bool AIC::agglomerate() {
	return agglomerate(1E-6, 1E-15);
}

bool AIC::agglomerate(double fn_tol, double eps) {
	double gamma = -1;
	vector<vector<size_t>> P=this->getPartition(gamma);
	if (P.size() <= 1) return false;
	SF_ f_(f,P);
	size_t n = P.size();
	//cout << "P : " << P << endl;
	vector<VectorXd> x(n - 1);
	for (size_t j = 1; j < n; j++) {
		SF__ f__(f_, j);
		//cout << "j : " << j << endl;
		x[j-1] = min_norm_base(f__,fn_tol,eps);
		//cout << x[j - 1];
		gamma = max(gamma, -(x[j-1].minCoeff()));
	}
	//cout << "gamma : " << gamma << endl;
	for (size_t j = 1; j < n; j++) {
		for (size_t i = 0; i < j; i++) {
			//cout << x[j - 1](i) << endl;
			if (x[j - 1](i) <= -gamma + fn_tol) // instead of eps, = is required to handle -INFINITY
				//cout << "fusing " << P[j][0] << " " << P[i][0] << endl;
				merge(P[j][0], P[i][0], gamma);
		}
	}
	return true;
}

size_t HC::find(size_t i) const {
	return (parent[i]==i) ? i : find(parent[i]);
}

bool HC::merge(size_t i, size_t j, double gamma) {
	size_t iroot = find(i);
	size_t jroot = find(j);
	if (iroot == jroot) return false; // i and j already merged
	if (this->gamma.empty() || this->gamma.back() > gamma) {
		this->gamma.push_back(gamma); // update the list of critical values of similarity thresholds
	}
	if (rank[iroot] <= rank[jroot]) {
		parent[iroot] = jroot;
		children[jroot].push_back(iroot);
		weight[iroot] = gamma;
		if (rank[iroot] == rank[jroot])
			rank[jroot] += 1;
	} else {
		parent[jroot] = iroot;
		children[iroot].push_back(jroot);
		weight[jroot] = gamma;
	}
	return true;
}

vector<vector<size_t>> HC::getPartition(double gamma) const {
	LinkList* to_cluster = new LinkList(weight.size(),1);
	vector<vector<size_t>> clusters;
	do {
		queue<pair<size_t, size_t>> to_search;
		vector<size_t> cluster;
		to_search.push(make_pair(to_cluster->head, to_cluster->head));
		while (!to_search.empty()) {
			pair<size_t, size_t> p = to_search.front();
			size_t i = p.second;
			size_t j = parent[i];
			if (j != i && j != p.first && weight[i] > gamma) {
				to_search.push(make_pair(i,j));
			}
			for (auto j : children[i]) {
				if (j != p.first && weight[j] > gamma) {
					to_search.push(make_pair(i, j));
				}
			}
			cluster.push_back(i);
			to_cluster->remove(i);
			to_search.pop();
		}
		clusters.push_back(cluster);
	} while (to_cluster->next());
	return clusters;
}

double HC::similarity(size_t i, size_t j) const {
	if (i == j) {
		return INFINITY;
	} 
	else if (parent[i] == i && parent[j] == j) { // not necessary if the disjoint-set forest is connected, i.e., is a tree.
		return 0;
	}
	else if (parent[i] == j) {
		return weight[i];
	} 
	else if (parent[j] == i) {
		return weight[j];
	}
	else if (weight[i] >= weight[j]) {
		return similarity(parent[i], j);
	}
	else {
		return similarity(parent[j], i);
	}
}

vector<double> HC::getCriticalValues() {
	return gamma;
}

HC::HC(size_t n) {
	if (n < 2)
		runtime_error("size must be at least 2.");
	parent.resize(n);
	children.resize(n);
	weight.resize(n);
	rank.resize(n);
	for (size_t i = 0; i < n; i++) {
		parent[i] = i;
		weight[i] = -INFINITY;
		rank[i] = 0;
	}
}

CL::CL(size_t n, vector<size_t> first_node, vector<size_t> second_node, vector<double> gamma) : HC(n) {
	size_t i = 0, ilen = gamma.size(), edge_count=0;
	if (first_node.size() != ilen || second_node.size() != ilen)
		runtime_error("Input vectors must have same length.");
	for (auto i : sort_indexes(gamma)) {
		if (first_node[i] >= n || second_node[i] >= n)
			runtime_error("node index must be smaller than size.");
		if (gamma[i]>0) {
			if (merge(first_node[i], second_node[i], gamma[i]))
				edge_count++;
		}
		if (edge_count >= n - 1) break; // graph is connected
	}
	if (edge_count < n - 1) {
		this->gamma.push_back(0); // 0 is a critical value when the graph is disconnected
	}
}