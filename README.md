# random-mixed-local-search-TSP
A new local search based meta-heuristic algorithm for TSP. Paper: https://doi.org/10.3390/app9193985

We have studied the open-loop variant of TSP to find out which local-search operators work best on them. Along with well-known relocation and 2-opt local search methods, we propose another method: link swap. This code runs these methods on open-loop TSPs iteratively. 

The code is written in C. It takes input data as TSPLIB format (https://www.or.uni-bonn.de/lectures/ws17/co_exercises/programming/tsp/tsp95.pdf). Output is only instance name and corresponding path length. But there are other options in writing several details of the output now commented out. If needed, please play with them :)

data folder has some sample instances. More datasets are available: https://cs.uef.fi/ml/tsp/

How to run? (simple)
Keep your instance in the data folder. Your instance should have locations (x and y co-ordinates) of all points in TSPLIB format. You can use Euclidean distance or Haversine distance. Both functions are written in the code. Compile with a gcc compiler and run. Will get the output in output folder.