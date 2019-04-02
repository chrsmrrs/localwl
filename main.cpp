/**********************************************************************
 * Copyright (C) 2017 Christopher Morris <christopher.morris@udo.edu>
 *********************************************************************/

#include <cstdio>
#include "src/AuxiliaryMethods.h"
#include "src/GenerateThree.h"

#include <iostream>
#include <chrono>

using namespace std::chrono;
using namespace std;

using namespace std;

int main() {

    string graph_database_name = "ENZYMES";
    cout << graph_database_name << endl;
    for (uint i = 0; i < 5; ++i) {
        cout << i << endl;
        GraphDatabase gdb = AuxiliaryMethods::read_graph_txt_file(graph_database_name);
        gdb.erase(gdb.begin());
        vector<int> classes = AuxiliaryMethods::read_classes(graph_database_name);
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        GenerateThree::GenerateThree generator_1(gdb);
        GramMatrix gm = generator_1.compute_gram_matrix(i, true, "local");
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(t2 - t1).count();
        cout << duration << endl;
        AuxiliaryMethods::write_libsvm(gm, classes, "MUTAG" + graph_database_name + "_" + to_string(i) + ".gram");
    }

    return 0;
}

