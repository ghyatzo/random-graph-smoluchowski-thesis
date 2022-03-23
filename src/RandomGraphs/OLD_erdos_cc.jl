#### OLD CODE ----- 

using Distributions
using Random
using LinearAlgebra
using Graphs
using StatsBase
using GLMakie
using GraphMakie


weight(d, dhat, m, i, j) = @views dhat[i]*dhat[j]*(1 - (d[i]*d[j])/2m)
weight(d, dhat, m, e::Union{Tuple, Pair}) = weight(d, dhat, m, e...)

function isgraphical(d)
    iseven(sum(d)) || return false
    cumsum = 0
    minsum = mapreduce(i->min(i, d[i]), +, eachindex(d))
    for i ∈ eachindex(d)
        cumsum += d[i]
        minsum -= min(i, d[i])
        (cumsum <= i*(i-1) + minsum) || return false
    end
    return true
end
function sample_graphical_sequence(dist, n)
    seq = rand(dist, n)
    return isgraphical(seq) ? seq : sample_graphical_sequence(dist, n)
end

edgepairs(n) = [(u, v) for u in 1:n-1 for v in u+1:n]
get_pair_linear_idx(i, j, n) = begin
    @assert (1 <= i < n) && (i < j <= n)
    return ((i-1)*(2n-i) ÷ 2)+j-i
end
get_pair_linear_idx(t::Tuple{Int, Int}, n) = get_pair_linear_idx(t[1], t[2], n)

function procedureA(d)
    @assert isgraphical(d)
    E = Set{Int}()
    dhat = deepcopy(d)
    n = length(d)
    m = sum(d) ÷ 2
    edgeset = edgepairs(n)

    vw = map(pair -> weight(d, dhat, m, pair), edgeset)
    sampleset = filter(i -> !iszero(vw[i]), eachindex(edgeset))

    eid = 0
    while length(sampleset) > 1
        @views eid = sample(sampleset, ProbabilityWeights(vw[sampleset]))
        eid ∈ E && continue

        i, j = edgeset[eid]
        push!(E, eid)
        dhat[i] -= 1
        dhat[j] -= 1

        @inbounds for k ∈ reverse(eachindex(sampleset))
            edge = edgeset[sampleset[k]]
            !(i ∈ edge || j ∈ edge) && continue

            w = weight(d, dhat, m, edge)
            w == 0 ? popat!(sampleset, k) : vw[k] = w
        end
    end

    !isempty(sampleset) && push!(E, sampleset...)

    length(E) < m && return false
    return E
end

function component_size_distribution(c, N)
    p = c*inv(n)
    dseq = sample_graphical_sequence(Poisson(c), N)
    E = procedureA(dseq)
    edges = map(eid -> edgeset[eid], collect(E))

    g = Graphs.SimpleGraph(n)
    for e in edges
        Graphs.add_edge!(g, Graphs.Edge(e...))
    end

    cc = Graphs.connected_components(g)
    cc_lengths = sort(length.(cc))

    dist = map(n -> count(==(n), cc_lengths), 1:maximum(cc_lengths))
    return dist
end


## Step 1, choose c, we want c to go from 0..1
## Step 2, set p = c/n and let n get very big (maybe huge steps) like 100k 1M, 10M, 100M
## Step 3, sample n points from a Poisson distribution to be our degree distribution.
## Step 4, use the procedure A to generate our graphs (pick the edges) such that we have the graph
## Step 5, count the connected components

# Step 1 & 2
c = 1.
n = 10000

p = c*inv(n);
ne = floor(Int, c*(n - 1)/2)

g = component_size_distribution(c, n)
g3 = erdos_renyi(n, ne)
g4 = expected_degree_graph(fill(c, n))

ntests = 1000
nedges = zeros(ntests)
CC_L = zeros(ntests)
C1 = zeros(ntests)
C2 = zeros(ntests)
for i = 1:ntests
    g = erdos_renyi(n, p)
    nedges[i] = g.ne
    cc = connected_components(g)
    ccl = sort(map(length, cc), rev=true)
    CC_L[i] = length(cc)
    C1[i] = ccl[1]
    C2[i] = ccl[2]
end
mean(nedges)
mean(CC_L)
mean(C1)
mean(C2)
std(nedges)
hist(nedges)

## Count connected comps
cc = Graphs.connected_components(g)
cc_lengths = sort(map(length, cc))

cc3 = Graphs.connected_components(g3)
cc_lengths3 = sort(map(length, cc3))

cc4 = Graphs.connected_components(g4)
cc_lengths4 = sort(map(length, cc4))






