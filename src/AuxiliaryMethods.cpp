/**********************************************************************
 * Copyright (C) 2019 Christopher Morris <christopher.morris@udo.edu>
 *********************************************************************/

#include "AuxiliaryMethods.h"

using Eigen::IOFormat;
using Eigen::MatrixXd;
using namespace std;

namespace AuxiliaryMethods {
    vector<int> split_string(string s) {
        vector<int> result;
        stringstream ss(s);

        while (ss.good()) {
            string substr;
            getline(ss, substr, ',');
            result.push_back(stoi(substr));
        }

        return result;
    }

    GraphDatabase read_graph_txt_file(string data_set_name) {
        string line;

        string path = "/home/morris/localwl";

        vector<uint> graph_indicator;
        ifstream myfile(
                path + "/data_sets/" + data_set_name + "/" + data_set_name +
                "_graph_indicator.txt");
        if (myfile.is_open()) {
            while (getline(myfile, line)) {
                graph_indicator.push_back(stoi(line));
            }
            myfile.close();
        } else {
            printf("%s", "!!! Unable to open file 1 !!!\n");
            exit(EXIT_FAILURE);
        }

        uint num_graphs = graph_indicator.back() + 1;

        // Get labels from for each node.
        bool label_data = true;
        string label;
        Labels node_labels;
        ifstream labels(
                path + "/data_sets/" + data_set_name + "/" + data_set_name + "_node_labels.txt");
        if (labels.is_open()) {
            while (getline(labels, label)) {
                node_labels.push_back(stoul(label));
            }
            myfile.close();
        } else {
            label_data = false;
        }

        GraphDatabase graph_database;
        unordered_map<int, int> offset;
        int num_nodes = 0;

        // Add vertices to each graph in graph database and assign labels.
        for (uint i = 0; i < num_graphs; ++i) {
            pair<int, int> p(i, num_nodes);
            offset.insert(p);
            unsigned long s = count(graph_indicator.begin(), graph_indicator.end(), i);

            Labels l;
            if (label_data) {
                for (unsigned long j = num_nodes; j < s + num_nodes; ++j) {
                    l.push_back(node_labels[j]);
                }
            }

            num_nodes += s;
            EdgeList edge_list;

            Graph new_graph(false, s, edge_list, l);
            graph_database.push_back(new_graph);
        }


        // Get labels from for each node.
        bool edge_label_data = true;
        Labels edge_labels;
        ifstream elabels(path + "/data_sets/" + data_set_name + "/" + data_set_name + "_edge_labels.txt");
        if (elabels.is_open()) {
            while (getline(elabels, label)) {
                edge_labels.push_back(stoul(label));
            }
            myfile.close();
        } else {
            edge_label_data = false;
        }

        // Insert edges for each graph.
        vector<EdgeLabels> edge_label_vector;
        for (uint i = 0; i < num_graphs; ++i) {
            edge_label_vector.push_back(EdgeLabels());
        }


        uint c = 0;
        vector<int> edges;
        ifstream edge_file(path + "/data_sets/" + data_set_name + "/" + data_set_name + "_A.txt");
        if (edge_file.is_open()) {
            while (getline(edge_file, line)) {
                vector<int> r = split_string(line);

                uint graph_num = graph_indicator[r[0] - 1];
                uint off = offset[graph_num];
                Node v = r[0] - 1 - off;
                Node w = r[1] - 1 - off;

                if (!graph_database[graph_num].has_edge(v, w)) {
                    graph_database[graph_num].add_edge(v, w);
                }

                if (edge_label_data) {
                    edge_label_vector[graph_num].insert({{make_tuple(v,w), edge_labels[c]}});
                    edge_label_vector[graph_num].insert({{make_tuple(w,v), edge_labels[c]}});
                }

                edges.push_back(stoi(line));
                c++;

            }
            edge_file.close();
        } else {
            printf("%s", "!!! Unable to open file 2!!!\n");
            exit(EXIT_FAILURE);
        }

        if (edge_label_data) {
            for (uint i = 0; i < num_graphs; ++i) {
                graph_database[i].set_edge_labels(edge_label_vector[i]);
            }
        }

        return graph_database;
    }

    vector<int> read_classes(string data_set_name) {
        string line;

        string path = "/home/morris/localwl";
        vector<int> classes;

        ifstream myfile(
                path + "/data_sets/" + data_set_name + "/" + data_set_name +
                "_graph_labels.txt");
        if (myfile.is_open()) {
            while (getline(myfile, line)) {
                classes.push_back(stoi(line));
            }
            myfile.close();
        } else {
            printf("%s", "!!! Unable to open file !!!\n");
            exit(EXIT_FAILURE);
        }

        return classes;
    }

    void write_gram_matrix(const GramMatrix &gram_matrix, string file_name) {
        const IOFormat CSVFormat(10, 1, ", ", "\n");
        string path = "/home/morris/localwl/svm/SVM/src/PRO/";

        ofstream file(path + file_name.c_str());

        // Convert sparse matrix to dense matrix to write it out to a file.
        MatrixXd dense_gram_matrix(gram_matrix);
        file << dense_gram_matrix.format(CSVFormat);

        file.close();
    }

    void write_sparse_gram_matrix(const GramMatrix &gram_matrix, string file_name) {
        saveMarket(gram_matrix, file_name);
    }


    void write_libsvm(const GramMatrix &gram_matrix, const vector<int> classes, std::string filename) {
        MatrixXd dense_gram_matrix(gram_matrix);

        int size = classes.size();
        MatrixXd norm(size, size);


        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                double x = sqrt(dense_gram_matrix(i, i)) * sqrt(dense_gram_matrix(j, j));
                if (x != 0) {
                    norm(i, j) = dense_gram_matrix(i, j) / x;
                } else {
                    norm(i, j) = 0.0;
                }
            }
        }

        string path = "/home/morris/localwl/svm/SVM/src/EXP/";

        ofstream file(path + filename);

        if (file.is_open()) {
            for (int i = 0; i < size; i++) {


                file << classes[i] << " 0:" << (i + 1);
                for (int c = 0; c < size; c++) {
                    file << " " << (c + 1) << ":" << norm(i, c);
                }
                file << std::endl;
            }
            file.close();
        } else {
            printf("%s", "!!! Unable to open file 3!!!\n");
            exit(EXIT_FAILURE);
        }
    }


    Label pairing(const Label a, const Label b) {
        return a >= b ? a * a + a + b : a + b * b;
    }
}
