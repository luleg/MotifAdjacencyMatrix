# Building Motif Adjacency Matrices

This is an efficient C++ software for building Motif Adjacency Matrices (MAM) of networks, for a range of motifs/graphlets:

* All 3-node motifs (Triangles and Wedges),
* All 4-node motifs containing a 4-node cycle (Quadrangles),
* It can also build an undirected network from a directed one.

The exhaustive list of admissible motifs to be found in *GraphletIdentifiers.pdf*.

## Using the Software

### Requirements

A bash shell and a `g++` compiler are enough to compile and use the software.

Tested in an `Ubuntu 18.04` environment emulated via a `Windows Subsystem for Linux 1`, with `gcc version 7.5.0` as compiler.

*While this implementation uses the [SNAP Software v6.0](https://snap.stanford.edu/snap/download.html), the present release contains the required files from the SNAP, and can be used as a standalone.*

### Installation

On a bash command, at the root of the folder:

```bash
cd src/BuildingMAM
make
cd ..
```

### Usage

The following command in the root folder explains how to run the software.
```bash
src/BuildingMAM/buildMAM -h
```

#### Detailed Usage

<img align="right" src="https://github.com/luleg/MotifAdjacencyMatrix/blob/main/Images/BiFan.png" width="8%" height="8%">

The software comes with examples in the folder *Data*. Assume we aim to build the MAM of the network *Data/Hierachy.nwk* upon the Bifan motif (shown on the right), and store the output MAM in a file *Data/Hierachy_Bifan.txt*. In the root folder:

```bash
src/BuildingMAM/buildMAM -i Data/Hierachy.nwk -m Q204 -otxt Data/Hierachy_Bifan.txt -nthreads 8
```

where:

<img align="right" src="https://github.com/luleg/MotifAdjacencyMatrix/blob/main/Images/toyGraph.png" width="15%" height="15%">

* ```-i Data/Hierachy.nwk``` provides the path to the input network.

:warning: The network must be an edgelist with integer nodes and only two columns (i.e. neither weights nor timestamps), separated by a space charcater, as shown on the right.

* ```-m Q204``` indicates the motif to use is the Bifan, that is identified by *Q204* (see *GraphletIdentifiers.pdf*).
* ```-otxt Data/Hierachy_Bifan.txt``` provides the path to the output MAM. The MAM is an edgelist of three columns, correponding to source, target and weight.

:bulb: An edge of weight ```w``` between two nodes ```a``` and ```b``` appears twice in the output file: as ```a b w``` and ```b a w```.

* ```nthreads 8``` 8 threads are used to build the MAM.

:bulb: To use when the network is large and the motif is a wedge or a quadrangle. Number of threads should be 4 and 8.

The same result is obtained using the command:

```bash
src/BuildingMAM/buildMAM -f Data/inputs.dat
```

:microscope: In *Data/inputs.dat* each line is one of the arguments mentioned above, ended by a character "*?*" to prevent portability issue.



## Definition of Motif Adjacency Matrices

Given a directed graph ***G = (V<sub>G</sub>, E<sub>G</sub>)***, and a graph with *k* nodes ***M = (V<sub>M</sub>, E<sub>M</sub>)*** (the so-called graphlet);
the MAM of ***G*** built upon ***M*** is an undirected weighted graph ***B= (V<sub>G</sub>, E<sub>B</sub>)*** such that the weight of an edge between two vertices ***u*** and ***v*** in ***V<sub>G</sub>*** is

<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=|\{\{w_1,...,w_{k-2}\} \subset V_G \setminus \{u,v\}\ \ \text{ s.t. the subgraph induced by }\ \ \{u,v, w_1,...,w_{k-2}\}\ \ \text{ in }\ \ G\ \ \text{ is isomorphic to }\ \ M \}|.">
</p>

Namely,
> A MAM counts the number of times that two nodes co-occur in an induced subgraph isomorphic to the graphlet.


More details about MAM to be found in

Austin Benson, [Tools for Higher-Order Network Analysis, Chap-II](https://arxiv.org/pdf/1802.06820.pdf) (2017)


## Tutorial

Two networks are provided in this software, namely *Data/Hierarchical.nwk* and *Data/Flows.nwk*. They both highlight interesting aspects of MAM, that we briefly investigate below.

:heavy_plus_sign: In addition to the current software, ```Python3``` with ```NetworkX2.5``` is used to build the networks, and ```Gephi``` to visualize them.

<img align="right" src="https://github.com/luleg/MotifAdjacencyMatrix/blob/main/Images/HierarchicalRaw.png" width="30%" height="30%">

### Hierarchical Network

The *Hierarchical* network has been created to imitate the behaviour of the *Twitter* who-follows-whom network. That is, a large number of anonymous users follows a handful of famous ones. The probability for an anonymous to follow a famous user is higher if the latter is famous in an field of interest for the former (e.g. a painter, a singer, a politician, etc.).

Famous users from a same field have a very high probability to follow each other, and a still high probability to follow famous user from other field (they may have met at a Gala night, in a TV show, etc.). Anonymous users interested in a same field have a slightly higher probability to follow each other than other anonymous users, since they have more chance to have met, in a museum, a concert, an online forum, etc. and may be interested in each other tweets.

A network with two fields was built, using a Stochastic Block Model whose odds and block sizes are shown below.

<img align="left" src="https://github.com/luleg/MotifAdjacencyMatrix/blob/main/Images/SBMHierarchy.png" width="25%" height="25%">

```python
import networkx as nx

BlockSizes = [250,250,25,25];
Pobabilities = [ [12e-3, 10e-3, 0.1,   2e-2 ],
                 [10e-3, 12e-3, 2e-2,  0.1  ],
                 [3e-3,  3e-3,  0.15,  0.015],
                 [3e-3,  3e-3,  0.015, 0.15 ] ]
G = nx.stochastic_block_model(BlockSizes,Pobabilities,directed=True)
nx.write_edgelist(G,'Data/Hierarchical.nwk',delimiter,' ',data=False)
```
The network is displayed at the top right, using the layout ```Force Atlas 2``` in ```Gephi```. Red, respectively dark blue, nodes correspond to famous users from one field, let say painting, respectively politics. Yellow, respectively light blue, nodes are anonymous users interested in painting, respectively politics. Users from different fields are roughly located in different part of the layout. However the network seems to consists in one homogeneous component, and without the colors, one cannot say that there are two fields, neither where to partition the network.

<img align="right" src="https://github.com/luleg/MotifAdjacencyMatrix/blob/main/Images/Q2252.png" width="8%" height="8%"> To highlight the two components, we used the software to build the MAM of this network upon motif *Q2252* (shown on the right):
```bash
src/BuildingMAM/buildMAM Data/Hierarchy.nwk -m Q2252 -otxt Data/Hierarchy_Q2252.txt
```
<img align="left" src="https://github.com/luleg/MotifAdjacencyMatrix/blob/main/Images/Hierarchical2252WithSingle.png" width="37%" height="37%"> This provides the network on the left, again displayed using the layout ```Force Atlas 2``` in ```Gephi```. Some nodes do not belong to any subgraph isomorphic to the *Q2252* graphlets, and are let disconnected. In the largest connected component, the two fields are highly separated, and a partition of the network in two components appears clearly, with only a few nodes being wrongly assigned.

<img align="right" src="https://github.com/luleg/MotifAdjacencyMatrix/blob/main/Images/Hierarchical2252_2Comp.png" width="30%" height="30%"> Furthermore, by removing the edges of weight 3 and below, we end up with two disconnected components that fit the expected fields extremely well, as shown on the right. Once again, the figure is obtained using the layout ```Force Atlas 2``` in ```Gephi``` (previous single nodes are not shown and disconnected nodes here come from thresholding the edges).


### Network with Flow

 <img align="left" src="https://github.com/luleg/MotifAdjacencyMatrix/blob/main/Images/FlowsRaw.png" width="20%" height="20%"> The second network, called *Flows*, consists in two blocks of nodes, internally densely connected, with a huge amount of edges from one block to the other, and only a few edges in the opposite direction. <img align="right" src="https://github.com/luleg/MotifAdjacencyMatrix/blob/main/Images/SBMFlows.png" width="20%" height="20%"> Such a network is interestig, because two different partitionings may be wanted. If one is interested in finding nodes involved in source/sink relationships, the partitioning that put all the nodes in one unique block is the most desirable, while if looking for groups of nodes involved in a same reinforcement process, the partitioning with two partitions, one per block densely connected, is prefered. We will investigate such a network using two graphlets, one for each kind of structures we have enunciated.

 We built the network using the Stochastic Block Model whose probabilities and block sizes are shown on the right, using ```NetworkX```:
```python
import networkx as nx

BlockSizes = [250,150];
Probabilities =[ [0.05, 0.2] , [1e-3, 0.075] ]
G = nx.stochastic_block_model(BlockSizes,Probabilities,directed=True)
nx.write_edgelist(G,'Data/Flows.nwk',delimiter,' ',data=False)
```
The resulting network is displayed on the top left, where red nodes are those from the smallest block (150 nodes), and green nodes belong to the largest one (250). The network visualisation was obtained by applying a ```Fruchterman Reingold``` layout followed by a ```Force Atlas 2``` layout in ```Gephi```.

<img align="right" src="https://github.com/luleg/MotifAdjacencyMatrix/blob/main/Images/Flows38.png" width="20%" height="20%"><img align="left" src="https://github.com/luleg/MotifAdjacencyMatrix/blob/main/Images/FFL.png" width="10%" height="10%">If one aims to find groups of nodes involved in source/sink relationships, a sensible motif upon which building the MAM can be motif *T38*, the Feed Forward loop, illustrated on the left. The corresponding MAM was obtained using the software:
```bash
src/BuildingMAM/buildMAM -i Data/Flows.nwk -m T38 -otxt Data/Flows_T38.txt
```
This results in the network highlighted on the right, whose spatialisation was done via a ```Fruchterman Reingold``` followed by a ```Force Atlas 2``` layouts in ```Gephi```. We observe that this MAM seems homogeneous, as it does not highlight block structure, which is the desired target.

_______

<img align="left" src="https://github.com/luleg/MotifAdjacencyMatrix/blob/main/Images/Flows98.png" width="25%" height="25%"><img align="right" src="https://github.com/luleg/MotifAdjacencyMatrix/blob/main/Images/Loop.png" width="10%" height="10%">On the other hand, a reinforcement process can be expressed by a motif containing a loop. We chose the motif *T98*, displayed on the right, to build the second MAM, which was generated with the command:
```bash
src/BuildingMAM/buildMAM -i Data/Flows.nwk -m T98 -otxt Data/Flows_T98.txt
```
The result is highlighted on the left, with the spatialisation obtained by our two usual layouts in ```Gephi```. This time, the network is not homogeneous. Instead, the two blocks can be distinguished, as desired: nodes from different blocks are located in different areas of the layout, and there is a slight bottleneck at the junction between the two blocks.

## Stay Update

A second version of this software is on going, which will contains (at least) two major improvements:
* Dealing with anchors. The current version makes no difference between the roles that the nodes play within the graphlet. But counting nodes that co-occur in a graphlet only if they play specific roles can be useful. For instance, in the *Hierarchical* network, it should allow to separate  anonymous users from famous ones.
* A even more efficient method for computing MAM built upon 3-node graphlets has been implemented using [GraphBlas](https://graphblas.github.io/), and will be integrated to this software.
