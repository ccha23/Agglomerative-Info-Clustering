#ifndef _ChowLiu_h_included_
#define _ChowLiu_h_included_


namespace IC {
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
		CL(size_t n, std::vector<size_t> first_node, std::vector<size_t> second_node, std::vector<double> gamma) : HC(n) {
			size_t i = 0, ilen = gamma.size(), edge_count=0;
			if (first_node.size() != ilen || second_node.size() != ilen)
				std::runtime_error("Input vectors must have same length.");
			for (auto i : sort_indexes(gamma)) {
				if (first_node[i] >= n || second_node[i] >= n)
					std::runtime_error("node index must be smaller than size.");
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
	};

}

#endif