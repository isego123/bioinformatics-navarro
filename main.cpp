#include <bits/stdc++.h>
#include "fastqloader.h"
#include "BigraphToDigraph.h"
using namespace std;

int main(void) {
    auto fastqs = loadFastqFromFile("./input/ref10000_simulatedreads.fastq");
    auto alignmentGraph = DirectedGraph::StreamGFAGraphFromFile("./input/ref10000_linear.gfa");
    return 0;
}