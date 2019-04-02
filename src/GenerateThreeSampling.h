/**********************************************************************
 * Copyright (C) 2017 Christopher Morris <christopher.morris@udo.edu>
 *********************************************************************/

#ifndef WLFAST_GENERATETHREESAMPLING_H
#define WLFAST_GENERATETHREESAMPLING_H

#include <cmath>
#include <unordered_map>
#include <queue>

#include "Graph.h"

using ThreeTuple = tuple<Node, Node, Node>;

using namespace GraphLibrary;

namespace GenerateThreeSampling {
    class GenerateThreeSampling {
    public:
        GenerateThreeSampling(const GraphDatabase &graph_database);

        GramMatrix
        compute_gram_matrix(const uint num_iterations, const bool use_labels, const uint num_samples);

        Graph generate_local_graph(const Graph &g, const uint num_iterations, const uint num_samples, bool use_labels);

        ~GenerateThreeSampling();

    private:
        GraphDatabase m_graph_database;

        // Manage indices of of labels in feature vectors.
        ColorCounter m_label_to_index;

        // Counts number of distinct labels over all graphs.
        ulong m_num_labels;

        // Computes labels for vertices of graph.
        ColorCounter
        compute_colors(const Graph &g, const uint num_iterations, const uint num_samples, const bool use_labels);


        // Get neighborhood of a node.
        void explore_neighborhood(const Graph &g, const ThreeTuple &triple, const uint num_iterations,
                                  unordered_map<ThreeTuple, uint> &triple_to_int, Graph &new_graph,
                                  unordered_map<Edge, uint> &edge_type,
                                  unordered_map<Edge, uint> &vertex_id,
                                  unordered_map<Edge, uint> &local, unordered_map<Node,Label> &node_label_map, const bool use_labels);


        Label compute_label(Label c_i, Label c_j, Label c_k, uint e_1, uint e_2, uint e_3);

    };
}

#endif //WLFAST_GENERATETHREESAMPLING_H
