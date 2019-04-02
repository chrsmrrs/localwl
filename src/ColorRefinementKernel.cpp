/**********************************************************************
 * Copyright (C) 2019 Christopher Morris <christopher.morris@udo.edu>
 *********************************************************************/

#include "ColorRefinementKernel.h"

namespace ColorRefinement {
    ColorRefinementKernel::ColorRefinementKernel(const GraphDatabase &graph_database) : m_graph_database(
            graph_database),
                                                                                        m_label_to_index(),
                                                                                        m_num_labels(0) {}

    GramMatrix
    ColorRefinementKernel::compute_gram_matrix(const uint num_iterations, bool use_labels, bool use_edge_labels) {
        vector<ColorCounter> color_counters;
        color_counters.reserve(m_graph_database.size());

        // Compute labels for each graph in graph database.
        uint i = 0;
        for (Graph &graph: m_graph_database) {
            cout << i << " " << m_graph_database.size() << endl;

            color_counters.push_back(move(compute_colors(graph, num_iterations, use_labels, use_edge_labels)));


            i++;
        }

        size_t num_graphs = m_graph_database.size();
        vector<S> nonzero_compenents;



        // Compute feature vectors.
        ColorCounter c;
        for (Node i = 0; i < num_graphs; ++i) {
            c = color_counters[i];

            for (const auto &j: c) {
                Label key = j.first;
                uint value = j.second;
                uint index = m_label_to_index.find(key)->second;
                nonzero_compenents.push_back(S(i, index, value));
            }
        }

        // Compute Gram matrix.
        GramMatrix feature_vectors(num_graphs, m_num_labels);
        feature_vectors.setFromTriplets(nonzero_compenents.begin(), nonzero_compenents.end());
        GramMatrix gram_matrix(num_graphs, num_graphs);
        gram_matrix = feature_vectors * feature_vectors.transpose();

        return gram_matrix;
    }

    ColorCounter ColorRefinementKernel::compute_colors(const Graph &g, const uint num_iterations, bool use_labels,
                                                       bool use_edge_labels) {
        size_t num_nodes = g.get_num_nodes();


        Labels coloring;
        Labels coloring_temp;

        // Assign labels to nodes.
        if (use_labels) {
            coloring.reserve(num_nodes);
            coloring_temp.reserve(num_nodes);
            coloring = g.get_labels();
            coloring_temp = coloring;
        } else {
            coloring.resize(num_nodes, 1);
            coloring_temp = coloring;
        }

        EdgeLabels edge_labels;
        if (use_edge_labels) {
            edge_labels = g.get_edge_labels();
        }

        ColorCounter color_map;
        if (use_labels) {
            for (Node v = 0; v < num_nodes; ++v) {
                Label new_color = coloring[v];

                ColorCounter::iterator it(color_map.find(new_color));
                if (it == color_map.end()) {
                    color_map.insert({{new_color, 1}});
                    m_label_to_index.insert({{new_color, m_num_labels}});
                    m_num_labels++;
                } else {
                    it->second++;
                }
            }
        }


        uint h = 1;
        while (h <= num_iterations) {
            // Iterate over all nodes.
            for (Node v = 0; v < num_nodes; ++v) {
                Labels colors;

                Nodes neighbors(g.get_neighbours(v));
                colors.reserve(neighbors.size() + 1);

                // New color of node v.
                Label new_color;
                if (!use_edge_labels) {
                    // Get colors of neighbors.
                    for (const Node &n: neighbors) {
                        colors.push_back(coloring[n]);
                    }
                    sort(colors.begin(), colors.end());
                    colors.push_back(coloring[v]);

                    // Compute new label using composition of pairing function of Matthew Szudzik to map two integers to on integer.
                    new_color = colors.back();
                    colors.pop_back();
                    for (const Label &c: colors) {
                        new_color = AuxiliaryMethods::pairing(new_color, c);
                    }
                    coloring_temp[v] = new_color;
                } else {

                    // Get colors of neighbors.
                    for (const Node &n: neighbors) {
                        const auto t = edge_labels.find(make_tuple(v, n));
                        colors.push_back(AuxiliaryMethods::pairing(coloring[n], t->second));
                        colors.push_back(coloring[n]);
                    }
                    sort(colors.begin(), colors.end());
                    colors.push_back(coloring[v]);

                    // Compute new label using composition of pairing function of Matthew Szudzik to map two integers to on integer.
                    new_color = colors.back();
                    colors.pop_back();
                    for (const Label &c: colors) {
                        new_color = AuxiliaryMethods::pairing(new_color, c);
                    }
                    coloring_temp[v] = new_color;
                }

                // Keep track how often "new_label" occurs.
                auto it(color_map.find(new_color));
                if (it == color_map.end()) {
                    color_map.insert({{new_color, 1}});
                    m_label_to_index.insert({{new_color, m_num_labels}});
                    m_num_labels++;
                } else {
                    it->second++;
                }
            }

            // Assign new colors.
            coloring = coloring_temp;
            h++;
        }


        return color_map;
    }

    ColorRefinementKernel::~ColorRefinementKernel() {}
}