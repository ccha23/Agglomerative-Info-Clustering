#ifndef _ChowLiu_h_included_
#define _ChowLiu_h_included_


namespace IC {
    /**
	Construct an approximate info-clustering solution by Chow-Liu tree approximation.
	*/
	class CL : public HC {
    private:
        size_t n;
        std::vector<size_t> first_node, second_node;
        std::vector<double> mutual_info;
        void agglomerate() {
            size_t i = 0, ilen = mutual_info.size(), edge_count=0;
//             if (first_node.size() != ilen || second_node.size() != ilen)
//     std::runtime_error("Input vectors must have same length.");
            for (auto i : sort_indexes(mutual_info)) {
// 				if (first_node[i] >= n || second_node[i] >= n)
// 					std::runtime_error("node index must be smaller than size.");
                if (merge(first_node[i], second_node[i], mutual_info[i]))
                    edge_count++;
            }
        }
	public:
		/**
		Construct the clusters for a graph given a list of non-negative weighted edges.
		@param size The size |V| of the graph.
		@param first_node Vector of first incident nodes of the edges.
		@param first_node Vector of second incident nodes of the edges.
		@param mutual_info Vector of the weights of the edges.
		@return The PSP of the graph.
		*/
		CL(size_t n, std::vector<size_t> first_node, std::vector<size_t> second_node, std::vector<double> mutual_info) : 
        n(n), first_node(first_node), second_node(second_node), mutual_info(mutual_info),
        HC(n) {
            agglomerate();
		}
        
        /**
		Construct the clusters for the given entropy function.
		@param h the entropy function
		@return The PSP of the graph.
		*/
        CL(IC::SF &h) : HC(h.size()) {
            for (size_t i = 0; i < h.size(); i++) {
                for (size_t j = 0; j < i; j++) {
                    first_node.push_back(i);
                    second_node.push_back(j);
                    double I = h(std::vector<size_t> {i}) 
                        + h(std::vector<size_t> {j}) 
                        - h(std::vector<size_t> {i, j});
                    mutual_info.push_back(I);
                }
            }
            agglomerate();
        }
	};

}

#endif