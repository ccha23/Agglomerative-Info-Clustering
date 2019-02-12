#include "PSP.h"
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <math.h>
#include <queue>
#include <algorithm>

using namespace std;

// Constructors

PSP::PSP(size_t size) {
	if (size < 2) // at least two nodes for non-trivial PSP
		runtime_error("size must be at least 2."); 
	this->size = size;
	node.resize(size);
	for (size_t i = 0; i < size; i++) {
		node[i] = new Node(i);
	}
}

PSP::PSP(size_t size, vector<size_t> first_node, vector<size_t> second_node,vector<double> weight) {
	if (size < 2)
		runtime_error("size must be at least 2.");
	this->size = size;
	node.resize(size);
	for (size_t i = 0; i < size; i++) {
		node[i] = new Node(i);
	}
	size_t i = 0, ilen = weight.size(); 
	if (first_node.size() != ilen || second_node.size() != ilen)
		runtime_error("Input vectors must have same length.");
	for (auto i : sort_indexes(weight)) {
		if (first_node[i] >= size || second_node[i] >= size)
			runtime_error("node index must be smaller than size.");
		if (weight[i]>0) {
			this->fuse(first_node[i], second_node[i], weight[i]);
		}
		if (edge.size() >= size - 1) break; // graph is connected
	}
	if (edge.size() < size-1) { 
		gamma.push_back(0); // 0 is a critical value when graph is disconnected
	}
}

PSP::~PSP() {
	for (size_t i = 0; i < node.size(); i++) {
		delete node[i];
	}
	for (size_t i = 0; i < edge.size(); i++) {
		delete edge[i];
	}
}

PSP::Node::Node(size_t i) {
	parent = i; 
	rank = 0; 
	gamma = 0;
}

PSP::Edge::Edge(size_t i, size_t j, double w) {
	node.first = i;
	node.second = j;
	weight = w;
}

size_t PSP::find(size_t i) const {
	return (node[i]->parent==i) ? i : find(node[i]->parent);
}

void PSP::fuse(size_t i, size_t j, double gamma) {
	size_t iroot = find(i);
	size_t jroot = find(j);
	if (iroot == jroot) return;
	if (this->gamma.empty() || this->gamma.back() > gamma) {
		this->gamma.push_back(gamma);
	}
	if (node[iroot]->rank <= node[jroot]->rank) {
		node[iroot]->parent = jroot;
		node[iroot]->gamma = gamma;
		if (node[iroot]->rank == node[jroot]->rank)
			node[jroot]->rank = node[jroot]->rank + 1;
	} else {
		node[jroot]->parent = iroot;
		node[jroot]->gamma = gamma;
	}
	Edge* e = new Edge(i, j, gamma);
	node[i]->edge.push_back(e);
	node[j]->edge.push_back(e);
	this->edge.push_back(e);
}

vector<size_t> PSP::getCluster(size_t i, double gamma) const {
	// BFS over the maximum weight forrest
	queue<pair<size_t,Edge*>> to_search;
	vector<size_t> cluster;
	to_search.push(make_pair(i,(Edge*) NULL));
	while (!to_search.empty()) {
		pair<size_t,Edge*> p = to_search.front();
		for (Edge* edge : node[p.first]->edge) {
			if (p.second != edge && edge->weight > gamma) {
				to_search.push(make_pair(edge->node.first + edge->node.second - p.first,edge));
			}
		}
		cluster.push_back(p.first);
		to_search.pop();
	}
	return cluster;
}


vector<vector<size_t>> PSP::getPartition(double gamma) const {
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
		LinkList(size_t size) {
			lnode.resize(size);
			lnode[0] = new LinkNode(0, 0, 1);
			for (size_t i = 1; i < size; i++) {
				lnode[i] = new LinkNode(i, i - 1, i + 1);
			}
		}
		~LinkList() {
			for (size_t i = 1, ilen = lnode.size(); i<ilen; i++) {
				delete lnode[i];
			}
		}
		void remove(size_t i) {
			if (head != i) {
				if (lnode[i]->prev>=0) {
					lnode[lnode[i]->prev]->next = lnode[i]->next;
				}
				if (lnode[i]->next < lnode.size()) {
					lnode[lnode[i]->next]->prev = lnode[i]->prev;
				}
			}
		}
		bool next() {
			if (lnode[head]->next >= lnode.size()) return false;
			head = lnode[head]->next;
			return true;
		}
	};
	
	LinkList* to_cluster = new LinkList(size);
	vector<vector<size_t>> clusters;
	do {
		queue<pair<size_t, Edge*>> to_search;
		vector<size_t> cluster;
		to_search.push(make_pair(to_cluster->head, (Edge*)NULL));
		while (!to_search.empty()) {
			pair<size_t, Edge*> p = to_search.front();
			for (Edge* edge : node[p.first]->edge) {
				if (p.second == edge) continue;
				pair<size_t,Edge*> q = make_pair(edge->node.first + edge->node.second - p.first,edge);
				if (edge->weight > gamma) to_search.push(q);
			}
			cluster.push_back(p.first);
			to_cluster->remove(p.first);
			to_search.pop();
		}
		clusters.push_back(cluster);
	} while (to_cluster->next());
	return clusters;
}


double PSP::similarity(size_t i, size_t j) const {
	if (i == j) {
		return INFINITY;
	} 
	else if (node[i]->parent == i && node[j]->parent == j) {
		return 0;
	}
	else if (node[i]->parent == j) {
		return node[i]->gamma;
	} 
	else if (node[j]->parent == i) {
		return node[j]->gamma;
	}
	else if (node[i]->gamma >= node[i]->gamma) {
		return similarity(node[i]->parent, j);
	}
	else {
		return similarity(node[j]->parent, i);
	}
}

vector<double> PSP::getCriticalValues() {
	return gamma;
}

std::ostream& operator<<(std::ostream &strm, const PSP &psp) {
	strm << "critical values: " << psp.gamma << endl;
	strm << "node, parent, rank, gamma, [edge:(first node,second node,weight)]... : " << endl;
	for (size_t i = 0; i < psp.size; i++) {
		strm << i << "," << psp.node[i]->parent << ","
		<< psp.node[i]->rank << "," << psp.node[i]->gamma << ",";
		for (size_t j = 0, jlen = psp.node[i]->edge.size(); j < jlen; j++) {
			strm << "(" << psp.node[i]->edge[j]->node.first << ","
				<< psp.node[i]->edge[j]->node.second << ","
				<< psp.node[i]->edge[j]->weight << ")";
		}
		strm << endl;
	}
	return strm;
}