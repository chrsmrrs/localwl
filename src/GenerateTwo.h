/**********************************************************************
 * Copyright (C) 2017 Christopher Morris <christopher.morris@udo.edu>
 *********************************************************************/

#ifndef WLFAST_GENERATETWO_H
#define WLFAST_GENERATETWO_H

#include <cmath>
#include <unordered_map>
#include <queue>

#include "Graph.h"

using TwoTuple = tuple<Node, Node>;

using namespace GraphLibrary;

namespace GenerateTwo {
    class GenerateTwo {
    public:
        GenerateTwo(const GraphDatabase &graph_database);

        GramMatrix compute_gram_matrix(const uint num_iterations, const bool use_labels, const string algorithm);

        Graph generate_local_graph(const Graph &g, const bool use_labels);

        Graph generate_global_graph(const Graph &g, const bool use_labels);

        Graph generate_global_graph_malkin(const Graph &g, const bool use_labels);

        ~GenerateTwo();

    private:
        GraphDatabase m_graph_database;

        // Computes labels for vertices of graph.
        ColorCounter
        compute_colors(const Graph &g, const uint num_iterations, const bool use_labels, const string algorithm);

        // Manage indices of of labels in feature vectors.
        ColorCounter m_label_to_index;

        // Counts number of distinct labels over all graphs.
        uint m_num_labels;
    };
}

#endif //WLFAST_GENERATETWO_H
