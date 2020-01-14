#!/bin/bash
g++ -std=c++11 ThreadReadAssertion.cpp CommonUtils.cpp fastqloader.cpp AlignmentGraph.cpp BigraphToDigraph.cpp main.cpp -o main
./main
