# Readme
Quick and dirty run down of the functions in this project

#### `ThesisCode.jl`
The main module. We import all the code from the other places and run test/plots and whatnot.

#### `utils.jl`
A bunch of helper functions
- `mass_norm(a::Vector)` (also in place version)
	+ normalise a distribution with respect to the total mass of the system, calculated as the first moment.
- `sparsify!(a, eps)`
	+ run through the whole array and puts to 0 all elements that are smaller than `eps`
- `u(z, ws)`
	+ builds the u function (a bit useless now)
- `weak_components(z)`

## Numerical

## Random Graphs

 UB = exp(c1*(x-1))exp(c2*(y-1))
 UR = exp(c3*(x-1))exp(c4*(y-1))

 C = At, are the rates at which we create new edges.
 In the most simple case, where we disallow edges of different types to join together, we have that
 The expected degree distribution for black edges then is (c1, 0), similarly the expected degree distribution for red edges is (0, c4)

 In the final graph, we encode the type of the vertex as `[n1..., n2...]`
 We then can define the degree distribution between different blocks.

 the aim now is to define the expected degree between all possible combinations of blocks.
 NUMBER OF COMPONENTS: d = 2
 We will generate d^2 = 4 sets of expected degrees sequences, in this case:
	- black -> black
	- black -> red
	- red 	-> black
 	- red   -> red

 The aim of this, is then to generate d^2 separate Chung-Lu graphs with d^2 different weights sets
 Each time, we generate the graphs in linear time and then we will merge all edges to the same graph.

 We now generate a realisation of this distribution of edges adapting the algorithm from https://doi.org/10.1007/978-3-642-21286-4_10
 that runs in O(N+M) where N is the number of nodes and M is the expected number of edges.

 we assume that the numbering of the vertices is `n = [n1..., n2...]`
 so the actual vertex number for the `k`th vertex in the block `nj` is `sum([ni for i = 1:(j-1)])+k` or `cumsum(n)[j-1]+k`