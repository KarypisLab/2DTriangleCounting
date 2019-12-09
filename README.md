# TriangleCounting
An MPI-based 2D parallel triangle counting algorithm for distributed-memory architectures. 


## Build requirements with versions currently used
 - CMake (v3.10.2). 
 - GCC (v8.2.0). 
 - OpenMPI (v4.0.0). 


## Configuring and building
Build the program by issuing the following commands:
```
make config cc=mpicc
make
```

## Other make commands
    make clean 
         Removes all object files but retains the configuration options.
   
    make distclean 
         Performs clean and completely removes the build directory.


## Running the program 
Run the program by the issuing the following commands based on the type of the input file: 
```
mpirun -np <#ranks> ./build/Linux-x86_64/mpitc  -iftype=metis <input_file>
mpirun -np <#ranks> ./build/Linux-x86_64/mpitc  -iftype=rmat -scale=<rmat_scale>
```

For more information about running the program, use the help command:
```
mpirun -np 1 ./build/Linux-x86_64/mpitc  -help
./build/Linux-x86_64/mpitc -help 

Usage: mpirun -np <#ranks> ./build/Linux-x86_64/mpitc [options] infile
 
 Options
  -tctype=text
     Specifies the type of triangle counting algorithm to use.
     Possible values are:
        mapjik2d      hash-map jik enumeration 2D [default]
 
  -iftype=text
     Specifies the format of the input file. 
     Possible values are:
        rmat    Generate an RMAT graph for a given scale
        metis   METIS format [default]
 
  -scale=int
     Specifies the scale for generating the RMAT graphs.
     Default value is 10.
 
  -dbglvl=int
     Specifies the level of debugging information to be displayed.
     Default value is 0.
 
  -help
     Prints this message.
```

The program supports two formats for its input files: 
- [Metis](http://www.cs.umn.edu/~metis) format used in the graph partitioning program.
- RMAT format generated using an external library called [Grappa](https://github.com/uwsampa/grappa) to generate these graphs. 

Note that the graph has to be undirected and it needs to include both pairs of
edges (i.e., (u,v) and (v,u)).

Here is the output of a sample run:
```
mpirun -np 4 ./build/Linux-x86_64/mpitc  -iftype=metis  data/small-graph.metis
./build/Linux-x86_64/mpitc -iftype=metis data/small-graph.metis 
./build/Linux-x86_64/mpitc -iftype=metis data/small-graph.metis 
./build/Linux-x86_64/mpitc -iftype=metis data/small-graph.metis 
./build/Linux-x86_64/mpitc -iftype=metis data/small-graph.metis 
Reading graph data/small-graph.metis...
Reading graph data/small-graph.metis...
Reading graph data/small-graph.metis...
Reading graph data/small-graph.metis...

-----------------
infile: data/small-graph.metis
per proc #nvtxs: 3
tctype: mapjik2d, otype: incd

Overall separate timer: 0.000668. 
#triangles: 16; rate: 0.0613 MT/sec
wall clock time:  0.001s
pre-processing:   0.000s
triangle counting:  0.000s

-----------------
RUNRESULT: data/small-graph.metis mapjik2d 10 16 0.0613 0.001 0.000 0.000
```

## Citing 
The parallel algorithm implemented is described in

[__"A 2D Parallel Triangle Counting Algorithm for Distributed-Memory Architectures"__ Ancy Sarah Tom and George
Karypis, Proceedings of the 48th International Conference on Parallel Processing (ICPP),
2019.](https://arxiv.org/abs/1907.09575)


