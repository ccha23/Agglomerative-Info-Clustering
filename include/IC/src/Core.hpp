#ifndef _Core_h_included_
#define _Core_h_included_

#include <queue>

/**
Represent vectors by strings.
*/
template <typename T>

std::ostream& operator<<(std::ostream &strm, const std::vector<T> &v) {
	strm << "[ ";
	for (size_t i = 0, ilen = v.size(); i < ilen; i++)
		strm << v[i] << " ";
	return strm << "]";
}

/**
Represent pairs by strings.
*/
template <typename T1,typename T2>
std::ostream& operator<<(std::ostream &strm, const std::pair<T1,T2> &p) {
	return strm << "(" << p.first << ", " << p.second << ")";
}

/**
Sort the indices of a vector in descending order of the corresponding elements.
@param v vector
@return The vector of sorted indices.
*/
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

	// initialize original index locations
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });

	return idx;
}

namespace IC { 

	/**
	Base functor for submodular function. Ground set V is {0 ... n-1} for some non-negative number n.
	*/
	class SF {
	public:
		/**
		The value of the submodular function evaluated on a set.
		@param B A subvector of elements of V.
		@return The value of the submodular function at B.
		*/
		virtual double operator() (const std::vector<size_t> &B) const = 0;
		/**
		@return The size of the ground set.
		*/
		virtual size_t size() const = 0;
        
		/**
		Compute the mutual information between two disjoint sets of nodes.
		@param B A subvector of elements of V.
        @param C A subvector of elements of \f$V\setminus B\f$.
		@return \f$I(Z_B\wedge Z_C) = h(B) + h(C) - h(B \cup C)\f$.
		*/
		double mutual_info(const std::vector<size_t> &B, 
                                   const std::vector<size_t> &C) const {
            std::vector<size_t> U(B);
            U.insert( U.end(), C.begin(), C.end() );
            return operator()(B) + operator()(C) 
                - operator()(U);
            
        }
        
		/**
		Compute the mutual information between two distinct nodes.
		@param i Node label.
        @param j Node label.
		@return \f$I(Z_i\wedge Z_j) = h(\{i\}) + h(\{j\}) - h(\{i,j\})\f$.
		*/
		double mutual_info(size_t i, size_t j) const {
            return operator()(std::vector<size_t> {i}) + operator()(std::vector<size_t> {j}) 
                - operator()(std::vector<size_t> {i,j});
        }        
	};    
    
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
			std::vector<LinkNode*> lnode;
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
		std::vector<double> gamma;
		std::vector<size_t> parent;
		std::vector<std::vector<size_t> > children;
		std::vector<double> weight;
		std::vector<size_t> rank;

		friend std::ostream& operator<<(std::ostream&, const HC&);

		/**
		Find the current cluster label of a node. O(log(|V|)) time.
		@param i node label.
		@return The current cluster label of node i.
		*/
		size_t find(size_t i) const {
			return (parent[i]==i) ? i : find(parent[i]);
		}

		/**
		Merge nodes i and j into the same cluster for all threshold values strictly smaller than gamma. 
		Complexity: O(log n).
		@param i node label.
		@param j node label.
		@param gamma similarity level.
		@return true if a merge is done, i.e., i and j were not connected before and gamma is smaller than existing critical values.
		*/
		bool merge(size_t i, size_t j, double gamma)  {
			size_t iroot = find(i);
			size_t jroot = find(j);
			if (iroot == jroot || // i and j already merged
               not (this->gamma.empty() || 
                    this->gamma.back() >= gamma)) // gamma is too large
                return false;
            if (this->gamma.empty() || this->gamma.back() > gamma)
                this->gamma.push_back(gamma); // update the list of distinct critical values
			if (rank[iroot] <= rank[jroot]) {
                // make jroot the parent of iroot
				parent[iroot] = jroot;
				children[jroot].push_back(iroot);
				weight[iroot] = gamma;
				if (rank[iroot] == rank[jroot])
					rank[jroot] += 1; // depth of subtree rooted at jroot increases by 1
			} else { // make iroot the parent of jroot for a more balanced tree
				parent[jroot] = iroot;
				children[iroot].push_back(jroot);
				weight[jroot] = gamma;
			}
			return true;
		}

	public:
		HC(size_t n) {
			if (n < 2)
				std::runtime_error("size must be at least 2.");
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

		/**
		Return the partition P of V, the non-singleton elements of which are the clusters at threshold gamma.
		Complexity: O(n).
		@param gamma similarity threshold.
		@return The partition P that gives the clusters at threshold gamma.
		*/
		std::vector<std::vector<size_t> > getPartition(double gamma) const {
			LinkList* to_cluster = new LinkList(weight.size(),1);
			std::vector<std::vector<size_t>> clusters;
			do {
				std::queue<std::pair<size_t, size_t>> to_search;
				std::vector<size_t> cluster;
				to_search.push(std::make_pair(to_cluster->head, to_cluster->head));
				while (!to_search.empty()) {
					std::pair<size_t, size_t> p = to_search.front();
					size_t i = p.second;
					size_t j = parent[i];
					if (j != i && j != p.first && weight[i] > gamma) {
						to_search.push(std::make_pair(i,j));
					}
					for (auto j : children[i]) {
						if (j != p.first && weight[j] > gamma) {
							to_search.push(std::make_pair(i, j));
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
        
		/**
		Return the coaseast partition of V generated. This is equivalent to getPartition(-INFINITY).
		Complexity: O(n).
		@param gamma similarity threshold.
		@return The partition P that gives the clusters at threshold gamma.
		*/
		std::vector<std::vector<size_t> > getCoarsestPartition() const {
            return getPartition(-INFINITY);
        }
        
        /**
		Return the cluster containing a node at a similarity threshold.
		Complexity: O(n).
        @param i node label.
		@param gamma similarity threshold.
		@return The cluster containing node i at threshold gamma.
		*/
		std::vector<size_t> getCluster(size_t i, double gamma) const {
            std::queue<std::pair<size_t, size_t>> to_search;
            std::vector<size_t> cluster;
            to_search.push(std::make_pair(i, i));
            while (!to_search.empty()) {
                std::pair<size_t, size_t> p = to_search.front();
                size_t i = p.second;
                size_t j = parent[i];
                if (j != i && j != p.first && weight[i] > gamma) {
                    to_search.push(std::make_pair(i,j));
                }
                for (auto j : children[i]) {
                    if (j != p.first && weight[j] > gamma) {
                        to_search.push(std::make_pair(i, j));
                    }
                }
                cluster.push_back(i);
                to_search.pop();
            }
			return cluster;
		}

		/**
		Return the similarity of node i and j. 
		Complexity: O(log n).
		@param i node label.
		@param j node label.
		@return The similarity of node i and j.
		*/
		double similarity(size_t i, size_t j) const {
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

		/**
		Return the set of critical values. 
		Complexity: O(n).
		@return The vector of critical values in descending order.
		*/
		std::vector<double> getCriticalValues() const {
			return gamma;
		}
        
        /**
		Return the set of critical values for clusters containing a node.
		Complexity: O(n).
		@param i node label.
		@return The vector of critical values in descending order for clusters containing i.
		*/
		std::vector<double> getCriticalValues(size_t i) const {
            std::vector<double> gamma;
            for (auto j : children[i]) { // add weights of children
                if (gamma.empty() || gamma.back() > weight[j])
                    gamma.push_back(weight[j]);
            }
            while (parent[i] != i) { // add weights of ancestors
                if (gamma.empty() || gamma.back() > weight[i])
                    gamma.push_back(weight[i]);
                i = parent[i];
            }
			return gamma;
        }
	};
}

#endif