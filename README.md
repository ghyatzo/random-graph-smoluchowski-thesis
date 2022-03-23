# random-graph-smoluchowski-thesis
Code I used in my master thesis, for the efficient generation of d-coloured graphs and subsequent plotting.

# Main files
Quick and dirty run down of the functions in this project.

## data
all the data obtained through the simulations.

## imgs
some pregenerated images that were used in the manuscript.

## src
The main files with the code

#### `utils.jl`
A bunch of helper functions
- `mass_norm(a::Vector)` (also in place version)
	+ normalise a distribution with respect to the total mass of the system, calculated as the first moment.
- `sparsify!(a; eps)`
	+ run through the whole array and puts to 0 all elements that are smaller than `eps`.
- `zerolift(a; eps)`
    + rises all elements that are zero to eps.
- `nanmap(a; eps)`
    + if a value is smaller than eps, make it a `NaN` (plotting helper function when dealing with logarithmic scales.)
- `massmultiply(a)`
    + returns the array with `k*a[k]` entries.
- `partitionN(N, eta)`
    + input: `N` the total amount of nodes, `eta` the distribution of subtypes
    + output: `n` a vector of concrete quantities, `n[i] = N * eta[i]`
- `vertex_type(v, n)`
    + given a vertex number `v`, returns the region in which it falls, and therefore the type.

### Numerical
Self sufficient code to integrate numerically the
- 1d burgers equation
- 1d smoluchowski system of odes
- 2d smoluchowski system of odes
(using DifferentialEquations pacakge)

### Plotting
For each type of plot we generated, there is a self sufficient file. You can run them on their own and images will be saved in the `imgs` folder

### Random Graphs

- `convolution_monocomponent.jl`
    + code to run the exact solution obtained through the implicit algebraic system of functional equations.
- `random_graph_monocomponent.jl`
    + an easy example for the monocomponent case
- `random_graph_multicomp.jl`
    + the main file. explained below
    

#### `random_graph_multicomp.jl`

here we implement the main algorithm that runs in `O((d+1)N)` time.
This algorithm was adapted from https://doi.org/10.1007/978-3-642-21286-4_10.

 `C = At`, are the rates at which we create new edges.
 In the most simple case, where we disallow edges of different types to join together, we have that
 The expected degree distribution for black edges then is (c1, 0), similarly the expected degree distribution for red edges is (0, c4)

- `chung_lu_edges(c1, c2, n1, n2)`
    + generates a set of pseudo edges, between `n1` nodes of one type and `n2` nodes of another type, with expected degree between them given by `c1` and `c2` respectively.

- `generate_multi_graph(t, kernel, n)`
    + input:
        * time `t`
        * kernel: defines the quantities `c_ij`
        * n: the partition in regions of the `N` nodes.
    + The idea is, for each pair `(i,j)` we generate a set of pseudoedges through `chung_lu_edges`
    + we run through all the pseudo-edgesets and wire them correctly wrt the type of nodes involved.
    
- `avg_cc_dist`
    + generates multiple graphs, read the component distribution and return an average of the results
    
- `avg_weak_cc_dist`
    + same as above but for the weaklt connected compoentns.
