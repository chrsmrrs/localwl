/**********************************************************************
 * Copyright (C) 2019 Christopher Morris <christopher.morris@udo.edu>
 *********************************************************************/

#include "AuxiliaryMethods.h"
#include "GenerateThreeSampling.h"


namespace GenerateThreeSampling {
    GenerateThreeSampling::GenerateThreeSampling(const GraphDatabase &graph_database) : m_graph_database(
            graph_database), m_label_to_index(), m_num_labels(0) {}


    GramMatrix
    GenerateThreeSampling::compute_gram_matrix(const uint num_iterations, const bool use_labels,
                                               const uint num_samples) {
        vector<ColorCounter> color_counters;
        color_counters.reserve(m_graph_database.size());

        uint i = 0;
        // Compute labels for each graph in graph database.
        for (Graph &graph: m_graph_database) {
            //cout << i << endl;
            color_counters.push_back(move(compute_colors(graph, num_iterations, num_samples, use_labels)));
            i++;
        }

        size_t num_graphs = m_graph_database.size();
        vector<S> nonzero_compenents;

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
        feature_vectors = feature_vectors * (1.0 / m_num_labels);
        GramMatrix gram_matrix(num_graphs, num_graphs);
        gram_matrix = feature_vectors * feature_vectors.transpose();

        return gram_matrix;
    }

    ColorCounter
    GenerateThreeSampling::compute_colors(const Graph &g, const uint num_iterations, const uint num_samples,
                                          const bool use_labels) {

        Graph tuple_graph(true);
        tuple_graph = generate_local_graph(g, num_iterations, num_samples, use_labels);

        size_t num_nodes = tuple_graph.get_num_nodes();

        Labels coloring;
        Labels coloring_temp;

        coloring.reserve(num_nodes);
        coloring_temp.reserve(num_nodes);
        coloring = tuple_graph.get_labels();
        coloring_temp = coloring;

        EdgeLabels edge_labels = tuple_graph.get_edge_labels();
        EdgeLabels vertex_id = tuple_graph.get_vertex_id();
        EdgeLabels local = tuple_graph.get_local();

        ColorCounter color_map;
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

        uint h = 1;
        while (h <= num_iterations) {
            // Iterate over all nodes.
            for (Node v = 0; v < num_nodes; ++v) {
                Labels colors_local;
                Labels colors_global;
                Nodes neighbors(tuple_graph.get_neighbours(v));
                colors_local.reserve(neighbors.size() + 1);
                colors_global.reserve(neighbors.size() + 1);

                // New color of node v.
                Label new_color;

                vector<vector<Label>> set_m_local;
                vector<vector<Label>> set_m_global;
                unordered_map<uint, uint> id_to_position_local;
                unordered_map<uint, uint> id_to_position_global;

                uint dl = 0;
                uint dg = 0;
                // Get colors of neighbors.
                for (const Node &n: neighbors) {
                    const auto t = edge_labels.find(make_tuple(v, n));
                    Label l = AuxiliaryMethods::pairing(coloring[n], t->second);

                    const auto type = local.find(make_tuple(v, n));

                    if (type->second == 1) {
                        const auto s = vertex_id.find(make_tuple(v, n));
                        const auto pos(id_to_position_local.find(s->second));
                        if (pos != id_to_position_local.end()) {
                            set_m_local[pos->second].push_back(l);
                        } else {
                            id_to_position_local.insert({{s->second, dl}});
                            set_m_local.push_back(vector<Label>());
                            set_m_local[dl].push_back(l);
                            dl++;
                        }
                    } else {
                        const auto s = vertex_id.find(make_tuple(v, n));
                        const auto pos(id_to_position_global.find(s->second));
                        if (pos != id_to_position_global.end()) {
                            set_m_global[pos->second].push_back(l);
                        } else {
                            id_to_position_global.insert({{s->second, dg}});
                            set_m_global.push_back(vector<Label>());
                            set_m_global[dg].push_back(l);
                            dg++;
                        }
                    }
                }

                for (auto &m: set_m_local) {
                    sort(m.begin(), m.end());
                    new_color = m.back();
                    m.pop_back();
                    for (const Label &c: m) {
                        new_color = AuxiliaryMethods::pairing(new_color, c);
                    }
                    colors_local.push_back(new_color);
                }
                sort(colors_local.begin(), colors_local.end());

                for (auto &m: set_m_global) {
                    sort(m.begin(), m.end());
                    new_color = m.back();
                    m.pop_back();
                    for (const Label &c: m) {
                        new_color = AuxiliaryMethods::pairing(new_color, c);
                    }
                    colors_global.push_back(new_color);
                }
                sort(colors_global.begin(), colors_global.end());

                for (auto &c: colors_global) {
                    colors_local.push_back(c);
                }
                colors_local.push_back(coloring[v]);

                // Compute new label using composition of pairing function of Matthew Szudzik to map two integers to on integer.
                new_color = colors_local.back();
                colors_local.pop_back();
                for (const Label &c: colors_local) {
                    new_color = AuxiliaryMethods::pairing(new_color, c);
                }
                coloring_temp[v] = new_color;

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


    Graph
    GenerateThreeSampling::generate_local_graph(const Graph &g, const uint num_iterations, const uint num_samples,
                                                bool use_labels) {
        // New tuple graph.
        Graph new_graph(false);
        size_t num_nodes = g.get_num_nodes();

        // Random device for sampling node in the graph.
        random_device rand_dev;
        mt19937 mt(rand_dev());
        uniform_int_distribution<Node> uniform_node_sampler(0, num_nodes - 1);

        unordered_map<Edge, uint> edge_type;
        // Manages vertex ids
        unordered_map<Edge, uint> vertex_id;
        // Manage type of neighborhood.
        unordered_map<Edge, uint> local;
        // Map tuples to node ids in new graph.
        unordered_map<ThreeTuple, uint> triple_to_int;
        // Node labels of new graph.
        unordered_map<Node, Label> node_label_map;

        Labels node_labels;
        if (use_labels) {
            node_labels = g.get_labels();
        }

        for (uint c = 0; c < num_samples; ++c) {
            // Sample a triple.
            Node t_0 = uniform_node_sampler(mt);
            Node t_1 = uniform_node_sampler(mt);
            Node t_2 = uniform_node_sampler(mt);

            ThreeTuple triple = make_tuple(t_0, t_1, t_2);

            // Add node to graph representing the triple.
            triple_to_int.insert({{triple, new_graph.get_num_nodes()}});
            Node new_node = new_graph.add_node();

            Label l;
            if (use_labels) {
                l = compute_label(node_labels[t_0], node_labels[t_1], node_labels[t_2], g.has_edge(t_0, t_1),
                                  g.has_edge(t_0, t_2), g.has_edge(t_1, t_2));
            } else {
                l = g.has_edge(t_0, t_1) + g.has_edge(t_0, t_2) + g.has_edge(t_1, t_2) + 1;
            }
            node_label_map.insert({{new_node, l}});


            // Get neighborhood around the triple up to depth.
            explore_neighborhood(g, triple, num_iterations, triple_to_int, new_graph, edge_type, vertex_id, local,
                                 node_label_map, use_labels);

        }


        Labels new_node_labels;
        uint n = new_graph.get_num_nodes();
        for (uint i = 0; i < n; ++i) {
            new_node_labels.push_back(node_label_map.find(i)->second);
        }

        new_graph.set_edge_labels(edge_type);
        new_graph.set_vertex_id(vertex_id);
        new_graph.set_local(local);
        new_graph.set_labels(new_node_labels);


        return new_graph;
    }

    void
    GenerateThreeSampling::explore_neighborhood(const Graph &g, const ThreeTuple &triple, const uint num_iterations,
                                                unordered_map<ThreeTuple, uint> &triple_to_int, Graph &new_graph,
                                                unordered_map<Edge, uint> &edge_type,
                                                unordered_map<Edge, uint> &vertex_id,
                                                unordered_map<Edge, uint> &local,
                                                unordered_map<Node, Label> &node_label_map, const bool use_labels) {

        // Manage depth of node in k-disk.
        unordered_map<ThreeTuple, uint> depth;

        // Node is in here iff it has been visited.
        unordered_set<ThreeTuple> visited;

        // Queue of nodes for DFS.
        queue<ThreeTuple> queue;

//        unordered_map<ThreeTuple, uint> triple_to_int;
//        triple_to_int.insert({{triple, new_graph.get_num_nodes() - 1}});

        // Push center tuple to queue.
        queue.push(triple);

        visited.insert(triple);
        depth.insert({{triple, 0}});

        Labels labels;
        if (use_labels) {
            labels = g.get_labels();
        }

        while (!queue.empty()) {
            ThreeTuple q(queue.front());
            queue.pop();

            Node v = get<0>(q);
            Node w = get<1>(q);
            Node u = get<2>(q);
            Node current_node = triple_to_int.find(make_tuple(v, w, u))->second;
            uint current_depth = depth.find(q)->second;

            if (current_depth <= num_iterations) {
                vector<ThreeTuple> neighbours;

                // Exchange first node.
                Nodes v_neighbors = g.get_neighbours(v);
                for (const auto &v_n: v_neighbors) {
                    auto t = triple_to_int.find(make_tuple(v_n, w, u));

                    Node new_node;
                    if (t == triple_to_int.end()) {
                        triple_to_int.insert({{make_tuple(v_n, w, u), new_graph.get_num_nodes()}});
                        new_node = new_graph.add_node();
                    } else {
                        new_node = t->second;
                    }
                    new_graph.add_edge(current_node, new_node);
//                    new_graph.add_edge(new_node, current_node);


                    Label new_label;
                    if (use_labels) {
                        new_label = compute_label(labels[v_n], labels[w], labels[u], g.has_edge(v_n, w),
                                                  g.has_edge(v_n, u), g.has_edge(w, u));
                    } else {
                        new_label = g.has_edge(v_n, w) + g.has_edge(v_n, u) + g.has_edge(w, u);
                    }

                    node_label_map.insert({{new_node, new_label}});
                    edge_type.insert({{make_tuple(current_node, new_node), 1}});
                    edge_type.insert({{make_tuple(new_node, current_node), 1}});
                    vertex_id.insert({{make_tuple(current_node, new_node), v_n}});
                    vertex_id.insert({{make_tuple(new_node, current_node), v_n}});
                    local.insert({{make_tuple(current_node, new_node), 1}});
                    local.insert({{make_tuple(new_node, current_node), 1}});

                    neighbours.push_back(make_tuple(v_n, w, u));
                }

                // Exchange second node.
                Nodes w_neighbors = g.get_neighbours(w);
                for (const auto &w_n: w_neighbors) {
                    auto t = triple_to_int.find(make_tuple(v, w_n, u));

                    Node new_node;
                    if (t == triple_to_int.end()) {
                        triple_to_int.insert({{make_tuple(v, w_n, u), new_graph.get_num_nodes()}});
                        new_node = new_graph.add_node();
                    } else {
                        new_node = t->second;
                    }
                    new_graph.add_edge(current_node, new_node);
//                    new_graph.add_edge(new_node, current_node);

                    Label new_label;
                    if (use_labels) {
                        new_label = compute_label(labels[v], labels[w_n], labels[u], g.has_edge(v, w_n),
                                                  g.has_edge(v, u), g.has_edge(w_n, u));
                    } else {
                        new_label = g.has_edge(v, w_n) + g.has_edge(v, u) + g.has_edge(w_n, u);
                    }
                    node_label_map.insert({{new_node, new_label}});
                    edge_type.insert({{make_tuple(current_node, new_node), 2}});
                    edge_type.insert({{make_tuple(new_node, current_node), 2}});
                    vertex_id.insert({{make_tuple(current_node, new_node), w_n}});
                    vertex_id.insert({{make_tuple(new_node, current_node), w_n}});
                    local.insert({{make_tuple(current_node, new_node), 1}});
                    local.insert({{make_tuple(new_node, current_node), 1}});

                    neighbours.push_back(make_tuple(v, w_n, u));
                }


                // Exchange third node.
                Nodes u_neighbors = g.get_neighbours(u);
                for (const auto &u_n: u_neighbors) {
                    auto t = triple_to_int.find(make_tuple(v, w, u_n));

                    Node new_node;
                    if (t == triple_to_int.end()) {
                        triple_to_int.insert({{make_tuple(v, w, u_n), new_graph.get_num_nodes()}});
                        new_node = new_graph.add_node();
                    } else {
                        new_node = t->second;
                    }
                    new_graph.add_edge(current_node, new_node);
//                    new_graph.add_edge(new_node, current_node);

                    Label new_label;
                    if (use_labels) {
                        new_label = compute_label(labels[v], labels[w], labels[u_n], g.has_edge(v, w),
                                                  g.has_edge(v, u_n), g.has_edge(w, u_n));
                    } else {
                        new_label = g.has_edge(v, w) + g.has_edge(v, u_n) + g.has_edge(w, u_n);
                    }

                    node_label_map.insert({{new_node, new_label}});
                    edge_type.insert({{make_tuple(current_node, new_node), 3}});
                    edge_type.insert({{make_tuple(new_node, current_node), 3}});
                    vertex_id.insert({{make_tuple(current_node, new_node), u_n}});
                    vertex_id.insert({{make_tuple(new_node, current_node), u_n}});
                    local.insert({{make_tuple(current_node, new_node), 1}});
                    local.insert({{make_tuple(new_node, current_node), 1}});

                    neighbours.push_back(make_tuple(v, w, u_n));
                }

                // Preprocessing.
                for (ThreeTuple &n: neighbours) {
                    if (visited.find(n) == visited.end()) {
                        depth.insert({{n, current_depth + 1}});
                        queue.push(n);
                        visited.insert(n);
                    }
                }
            }
        }
    }

    Label GenerateThreeSampling::compute_label(Label c_i, Label c_j, Label c_k, uint e_1, uint e_2, uint e_3) {


        Labels temp(
                {{AuxiliaryMethods::pairing(e_1 + e_2 + 1, c_i + 1),
                         AuxiliaryMethods::pairing(e_1 + e_3 + 1,
                                                   c_j + 1),
                         AuxiliaryMethods::pairing(e_2 + e_3 + 1, c_k + 1)
                 }
                });

        // sort(temp.begin(), temp.end());

        Label new_color = e_1 + e_2 + e_3;
        for (Label d: temp) {
            new_color = AuxiliaryMethods::pairing(new_color, d);
        }
        return new_color;
    }


    GenerateThreeSampling::~GenerateThreeSampling() {}
}
