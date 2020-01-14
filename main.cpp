#include <bits/stdc++.h>
#include "fastqloader.h"
#include "BigraphToDigraph.h"
using namespace std;

const int MAXN = 10000000;

vector<int> in[MAXN];
vector<int> out[MAXN];

/*
recursive check and update of node values in cyclic graphs if insertion is needed
*/
void propagate(int u, int v, AlignmentGraph& graph, vector<int> &Cold, vector<int> &Cnew) {
    if(Cnew[v] > 1 + Cnew[u]) {
        Cnew[v] = 1 + Cnew[u];
        for(int i = 0; i < out[v].size(); i++) {
            propagate(v, out[v][i], graph, Cold, Cnew);
        }
    }
}
/*
Navarro algorithm implementation - https://www.sciencedirect.com/science/article/pii/S0304397599003333
*/

int func(AlignmentGraph& graph, FastQ& query) {
    for(int i = 0; i < MAXN; i++) in[i].clear(), out[i].clear();
    auto componentOrder = graph.TopologicalOrderOfComponents();
    vector<vector<size_t>> components = componentOrder.first;
    vector<size_t> topo_sort;
    vector<char> all_chars;
    int curr_id = 0;
    map<size_t, int> mapping;
    
    for(auto comp: components) {
        for(int i = 0; i < comp.size(); i++) {
            mapping[comp[i]] = curr_id;
            for(int j = 0; j < graph.nodeToSeq[comp[i]].size(); j++) {
                topo_sort.push_back(curr_id);
                all_chars.push_back(graph.nodeToSeq[comp[i]][j]);
                if(j > 0) in[curr_id].push_back(curr_id - 1);
                if(j + 1 < graph.nodeToSeq[comp[i]].size()) out[curr_id].push_back(curr_id + 1);
                curr_id++;
            }
        }
    }
    
    for(auto comp: components) {
        for(int i = 0; i < comp.size(); i++) {
            
            for(int j = 0; j < graph.inNeighbors[comp[i]].size(); j++) {
                int curr = graph.inNeighbors[comp[i]][j];
                int fake_id = mapping[curr];
                out[fake_id + graph.nodeToSeq[curr].size() - 1].push_back(mapping[comp[i]]);
            }
            
            for(int j = 0; j < graph.outNeighbors[comp[i]].size(); j++) {
                int curr = graph.outNeighbors[comp[i]][j];
                in[mapping[curr]].push_back(mapping[comp[i]] + graph.nodeToSeq[comp[i]].size() - 1);
            }
        }
    }
    
    vector<int> v;
    vector<int> Cnew(topo_sort.size(), 0);
    vector<int> Cold(topo_sort.size(), 0);
    
    for(int i = 0; i < query.sequence.size(); i++) {
        for(int j = 0; j < topo_sort.size(); j++) {
            int g = query.sequence.size();
            size_t v = topo_sort[j];
            char curr_char = all_chars[j];
            if(query.sequence[i] == curr_char) {
                for(int k = 0; k < in[v].size(); k++) {
                    g = min(g, Cnew[in[v][k]]);
                }
                 g = min(g, i);
            } else {
                g = Cnew[v];
                for(int k = 0; k < in[v].size(); k++) {
                    g = min(g, Cnew[in[v][k]]);
                }
                g += 1;
            }
            Cold[v] = g;
        }
        
        for(int j = 0; j < topo_sort.size(); j++) {
            size_t v = topo_sort[j];
            Cnew[v] = Cold[v];
        }
        
        for(int j = 0; j < topo_sort.size(); j++) {
            size_t v = topo_sort[j];
            for(int k = 0; k < in[v].size(); k++) {
                propagate(in[v][k], v, graph, Cold, Cnew);
            }
        }
        
    }

    int sol = 1000000000;
    for(int i = 0; i < Cnew.size(); i++) {
        sol = min(sol, Cnew[i]);
    }
    return sol;
}

/*
main function
*/

int main(void) {
    auto fastqs = loadFastqFromFile("./input/ref10000_simulatedreads.fastq");
    auto alignmentGraph = DirectedGraph::StreamGFAGraphFromFile("./input/ref10000_twopath.gfa");
    
    auto timeStart = std::chrono::steady_clock::now();
    int count = 0;
    vector<int> results;
    
    for (auto f : fastqs) {
        int solution;
        auto tstart = std::chrono::steady_clock::now();
        if (count % 2 == 0) {
	        solution = func(alignmentGraph, f);
        } else {
            auto frev = f.reverseComplement();
    	    solution = func(alignmentGraph, frev);
        }
        auto ttend = std::chrono::steady_clock::now();
        size_t ms = std::chrono::duration_cast<std::chrono::microseconds>(ttend - tstart).count();
        cout << "vrijeme: " << (double) ms / 1000000 << endl;
        results.push_back(solution);
        cout << solution << endl;
        count++;
    }
	
    auto timeEnd = std::chrono::steady_clock::now();
    size_t ms = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeStart).count();
	std::cout << ms << " ms" << endl;
    
    for(int i = 0; i < results.size(); i++) {
        cout << results[i] << endl;
    }
    
    return 0;
} 
