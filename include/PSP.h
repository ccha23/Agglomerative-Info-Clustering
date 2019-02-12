#ifndef _PSP_h_included_
#define _PSP_h_included_

#include <iostream>
#include <vector>
#include <list>
#include <math.h>
#include <queue>
#include <numeric>

using namespace std;

/**
PSP is a data structure for the principal sequence of partition (PSP) under Chow-Liu tree approximation.
It is based on Kruskel's algorithm and the data structure of the disjoint-set forests with union by rank:
https://en.wikipedia.org/wiki/Kruskal%27s_algorithm

The structure contains:
1) a vector of critical values in descending order;
2) a maximum weighted forest implemented as adjacency list 
   https://en.wikipedia.org/wiki/Adjacency_list 
   that allows partitions at a similarity level to be computed in O(|V|) time;
3) an annotated disjoint-set forest that allows similarity level between two nodes to be computed in O(log(|V|)) time. 
   https://en.wikipedia.org/wiki/Disjoint-set_data_structure
   The structure basically represents a dendrogram with log(|V|) depth.

Complexity:
The structure is constructed in O((|E|+|V|)log(|V|)) time where V is the set of nodes to cluster and |E| 
is the number of pairs of nodes for which mutual information is available for the computation.
*/
class PSP {
	struct Edge {
		double weight;
		pair<size_t, size_t> node;

		/**
		Construct an Edge.
		@param i First incident node.
		@param j Second incident node.
		@param w Edge weight.
		*/
		Edge(size_t i, size_t j, double weight);
	};
	struct Node {
		size_t parent;
		size_t rank;
		double gamma;
		vector<Edge*> edge;

		/**
		Construct a Node to store properties associated with object i.
		@param i Node label.
		*/
		Node(size_t i);
	};
private:
	size_t size;
	vector<Node*> node;
	vector<Edge*> edge;
	vector<double> gamma;
	friend std::ostream& operator<<(std::ostream&, const PSP&);
protected:
	/**
	Find the current cluster label of a node. O(log(|V|)) time.
	@param i node label.
	@return The current cluster label of node i.
	*/
	size_t find(size_t i) const;

	/**
	Fuse two nodes together at a similarity level. If the two nodes are not already fused,
	the node with larger rank will become the parent of the node with smaller rank, to ensure
	maximum rank of O(log(|V|)).
	@param i node label.
	@param j node label.
	@param gamma similarity level.
	*/
	void fuse(size_t i, size_t j, double gamma);
public:
	/**
	Construct a PSP for a graph of isolated nodes.
	@param the size |V| of the graph
	*/
	PSP(size_t size);

	/**
	Construct a PSP for a graph given a list of non-negative weighted edges.
	@param size The size |V| of the graph.
	@param first_node Vector of first incident nodes of the edges.
	@param first_node Vector of second incident nodes of the edges.
	@param weight Vector of the weights of the edges.
	@return The PSP of the graph.
	*/
	PSP(size_t size, vector<size_t> first_node, vector<size_t> second_node, vector<double> weight);

	~PSP();

	/**
	Get the cluster of a node at similarity strictly larger than a threshold.
	O(|C|) time where C is the cluster containing the node.
	@param i Node label.
	@param gamma Threshold value.
	@return The cluster containing node i with similarly strictly larger than gamma.
	*/
	vector<size_t> getCluster(size_t i, double gamma) const;

	/**
	Get all the partition in PSP at similarity strictly larger than a threshold.
	O(|V|) time.
	@param gamma threshold value.
	@return The clusters with similarly strictly larger than gamma.
	*/
	vector<vector<size_t> > getPartition(double gamma) const;

	/**
	Return the similarity of two nodes. O(log(|V|)).
	@param i node label.
	@param j node label.
	@return The similarity of node i and j.
	*/
	double similarity(size_t i, size_t j) const;

	/**
	Return the set of critical values. O(|V|).
	@param i node label.
	@param j node label.
	@return The similarity of node i and j.
	*/
	vector<double> getCriticalValues();
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