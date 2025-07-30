# RACD

This repository contains C++ codes and datasets for the paper:

> Accelerated Coordinate Descent for Directed Densest Subgraph Discovery

## Environment

The codes of our in-depth study are implemented and tested under the following development environment:

- Hardware : Intel(R) Xeon(R) Gold 6338 CPU @ 2.00GHz and 512GB of memory.
- Operation System : Ubuntu 20.04.5 LTS (GNU/Linux 5.15.0-101-generic x86_64)

## Datasets

We use 15 real-world datasets from different domains, which are mainly downloaded from the [Stanford Network Analysis Platform](http://snap.stanford.edu/data/) and [Networks](http://konect.cc/networks/).



### Table: Datasets used in our experiments.

| Dataset | Category | \|V\| | \|E\| |
|---------|----------|-------|-------|
| Maayan-Lake (ML) | Foodweb | 183 | 2,494 |
| Moreno-Oz (MO) | Social | 217 | 2,672 |
| Traffic-Control (TC) | Infrastructure | 1,126 | 2,615 |
| Maayan-Figeys (MF) | Metabolic | 2,239 | 6,452 |
| OpenFlights (OF) | Infrastructure | 2,939 | 30,501 |
| Advogato (AD) | Social | 6,541 | 51,127 |
| Email-EU (EU) | Comments | 265,214 | 420,045 |
| Amazon (AM) | E-commerce | 403,394 | 3,387,388 |
| WikiTalk (WT) | Communication | 2,394,385 | 5,021,410 |
| Amazon-ratings (AR) | E-commerce | 3,376,962 | 5,838,041 |
| Baidu-Zhishi (BA) | Hyperlink | 2,141,300 | 17,794,839 |
| Wikipedia-Link (WL) | Hyperlink | 759,923 | 20,546,603 |
| EW-2013 (EW) | Social | 4,206,785 | 101,355,853 |
| Wiki-En (WE) | Hyperlink | 13,593,032 | 437,217,424 |
| SK-2005 (SK) | Web | 50,636,154 | 1,949,412,601 |


## How to Run the Codes


### A. Code Compilation


After cloning the codes from GitHub, use the following command to compile the codes in the repository :


```sh

cmake CMakeLists.txt

make ./DensestSubgraph

```


### B. Command Line Parameters


In terms of our work, one could run the code in the following format:

```sh
./DensestSubgraph -t d -a a -red appro-xy-core -alloc rac -ext cp -ver cp -dc t -seq t -map t -res t -width 5 -stats t -it 10 -wshrink t -initwcore t -path ./path/to/your/graph.txt -eps 0.1
```

In these parameters, the meaningful ones in our work are:
- `-t` : type of graph, `d` for directed and `u` for undirected. Here we always set it to `d`.
- `-a` : algorithm type, `a` for approximate and `e` for exact. Here we set it to `a`.
- `-eps` : the tolerance of the algorithm, which is set to `0` in exact algorithms and different values in approximate algorithms. Here we set it to `0.1`.
- `-stats` : whether to output the detailed statistics of the algorithm. Here we set it to `t`. To disable it, please set it to `f`.
- `-it` : the fixed iteration number. Here we set it to `10`.
- `-wshrink` : whether to use the shrinking technique. Here we set it to `t`.
- `-initwcore` : whether to use the initial w-core density technique. Here we set it to `t`.
- `-path` : the path of the input graph. You should set it to the path of your graph file. 
