# hdp_mixture
Hierarchical Dirichlet process mixture model

This is a library that implements a Gibbs sampler for a hierarchical Dirichlet process mixture of normals. It was originally developed as a method for this paper: http://www.nature.com/nmeth/journal/v14/n4/full/nmeth.4189.html . Accordingly, some code in the library is specific to nanopore sequencing, but there is a general purpose API as well.

hdp_mixture uses source from RANLIB for random number generation and sonLib for basic data structures:
- https://people.sc.fsu.edu/~jburkardt/c_src/ranlib/ranlib.html
- https://github.com/benedictpaten/sonLib

To build hdp_mixture, the sonLib repository must be cloned into the same directory as the hdp_mixture repository. The system must also have gcc version 4.9.1 or higher.
