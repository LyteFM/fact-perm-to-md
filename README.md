# fact-perm-to-md

This Code is originally from Adrian Fritz's Master Thesis: "A Heuristic for Cograph-Editing", Saarland University, 2015. The algorithm itselft is based on the following paper:

Christian Capelle, Michel Habib, Fabien Montgolfier. Graph Decompositions and Factorizing Permutations. Discrete Mathematics and Theoretical Computer Science, DMTCS, 2002, 5, pp.55-70.

I have isolated the parts that compute the modular decomposition from a given factorizing permutation and adapted it for directed graphs.

### License

Copyright (C) 2015 Adrian Fritz, 2018 Fynn Leitow

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.


### Compilation

Ensure that boost::graph v.1.55 or higher is installed:

$ cd PATH/TO/MD/build $ cmake .. $ make
or (for four kernels) $ make -j 4


### Usage

$ ./mod_dec -graph_as_string -fact_perm_as_sting > out.dot

convert to PNG via:

$ dot -Tpng out.dot -o out.png

Because I used this inside a Java-Project, the format of the graph is the same as a .toString() - call on a Graph in the JGraphT - Library. The factorising permutation is separated by Commata. An example can be found at the bottom of the modDecopm.cpp - file. Feel free to use it directly on a boost::graph or implement your own parser.
