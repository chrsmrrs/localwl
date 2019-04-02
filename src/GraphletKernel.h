/**********************************************************************
 * Copyright (C) 2017 Christopher Morris <christopher.morris@udo.edu>
 *
 * This file is part of globalwl.
 *
 * globalwl can not be copied and/or distributed without the express
 * permission of Christopher Morris.
 *********************************************************************/

#ifndef WLFAST_GRAPHLETKERNEL_H
#define WLFAST_GRAPHLETKERNEL_H

#include <algorithm>
#include <unordered_map>

#include "AuxiliaryMethods.h"
#include "Graph.h"

using Graphlet = Label;
using GraphletCounter = map<Graphlet, uint>;

using namespace GraphLibrary;

namespace GraphletKernel {
    class GraphletKernel {
    public:
        GraphletKernel(const GraphDatabase &graph_database);

        // Computes gram matrix for the graphlet kernel.
        GramMatrix compute_gram_matrix(bool use_labels);

        ~GraphletKernel();

    private:
        // Computes number of graphlets in graph.
        GraphletCounter compute_graphlet_count(const Graph &g, bool use_labels);

        // Manages graphs.
        GraphDatabase m_graph_database;

        // Manage indices of labels in feature vectors.
        ColorCounter m_label_to_index;

        // Counts number of distinct labels over all graphs.
        uint m_num_labels;
    };
}
#endif //WLFAST_GRAPHLETKERNEL_H
