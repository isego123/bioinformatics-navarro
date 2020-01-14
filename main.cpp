#include <bits/stdc++.h>
#include "fastqloader.h"
#include "BigraphToDigraph.h"
using namespace std;

vector<int> in[MAXN];
vector<int> out[MAXN];

void propagate(int u, int v, AlignmentGraph& graph, vector<int> &Cold, vector<int> &Cnew) {
    if(Cnew[v] > 1 + Cnew[u]) {
        Cnew[v] = 1 + Cnew[u];
        for(int i = 0; i < out[v].size(); i++) {
            propagate(v, out[v][i], graph, Cold, Cnew);
        }
    }
}

int main(void) {
    auto fastqs = loadFastqFromFile("./input/ref10000_simulatedreads.fastq");
    auto alignmentGraph = DirectedGraph::StreamGFAGraphFromFile("./input/ref10000_linear.gfa");
    return 0;
}
