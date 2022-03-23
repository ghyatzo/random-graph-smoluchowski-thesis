using Graphs
using LinearAlgebra: I, eigvals
using SparseArrays
using Base: Threads


includet("../utils.jl")

function component_composition(component::Vector{Int}, n::Vector{Int})
	# returns the object type of a component in terms of its constituents.
	# [k_1, k_2, ..., k_d]
	d = length(n)
	comp = zeros(Int, d)
	for i = 1:d
		comp[i] = count(v -> vertex_type(v, n) == i, component)
	end
	return comp
end
function count_components(components::Vector{Vector{Int}}, n::Vector{Int}, kmax::Int)
	# labels all components in a graph with the respective object type and counts how many there are.
	cc_comp = @. component_composition(components, (n,))
	cc_dist = zeros(Int, kmax, kmax)
	for comp in cc_comp
		cc_dist[1+comp[1], 1+comp[2]] += one(Int)
	end

	return cc_dist
end
count_components(g::SimpleGraph, n::Vector{Int}, kmax::Int) = count_components(connected_components(g), n, kmax)


function chung_lu_edges(c1::T, c2::T, n1::Int, n2::Int; seed::Int=-1) where T <: Real

	# We can't have more expected edges than there are nodes available.
    @assert all(zero(T) <= c1 <= n2) "c1 needs to be at least 0 and at most n2"
    @assert all(zero(T) <= c2 <= n1) "c2 needs to be at least 0 and at most n1"

    rng = Graphs.getRNG(seed)

    # pseudo edgeset initialisation
    E = Set{Tuple{Int, Int}}()

    S1 = n1*c1
    S2 = n2*c2

    for u = 1:(n1-1)
    	v = u + 1
    	p = min(c1 * c2 / S1, one(T))
    	while v <= n2 && p > zero(p)
    		if p != one(T)
    			v += floor(Int, log(rand(rng)) / log(one(T) - p))
    		end
    		if v <= n2
    			q = min(c1 * c2 / S2, one(T))
    			if rand(rng) < q / p
    				push!(E, (u, v))
    			end
    			p = q
    			v += 1
    		end
    	end
    end
    return E
end

function generate_multi_graph(t, A::Matrix, n::Vector{Int}; seed::Int=-1)
	# n is a vector that lists how many nodes of each type there are.
	d = length(n)
	@assert size(A, 1) == size(A, 2) == d

	# total amount of nodes.
	N = sum(n)

	# computes all the c_ij, n_i combos.
	Ω = [(t*A[i,j], n[i]) for i = 1:d, j = 1:d]

	EE = [Set{Tuple{Int, Int}}() for i = 1:d, j = 1:d]
	for i = 1:d
		for j = 1:d
			c1, n1 = Ω[i,j]
			c2, n2 = Ω[j,i]
			# for each (i,j), generate the corresponding pseudo-edgeset.
			EE[i, j] = chung_lu_edges(c1, c2, n1, n2; seed)
		end
	end

	# compute the node-regions intervals
	On = [0, cumsum(n)[1:end-1]...]
	g = SimpleGraph(N)
	for j = 1:d
		for i = 1:d
			for edge in EE[i, j]
				# for each pseudo edge, place its ends in the corresponding regions
				add_edge!(g, edge[1]+On[i], edge[2]+On[j])
			end
		end
	end

	return g
end

function avg_cc_dist(t, kernel::AbstractMatrix, n::Vector{Int}; nruns = 50, kmax=100)
	# computes the complete average d-dimensional connected components distribution up to kmax.
	lk = ReentrantLock()

	dz = inv(nruns)
	avg = sparse(0.0*I, kmax, kmax)
	Threads.@threads for _ = 1:nruns
		g = generate_multi_graph(t, kernel, n)
		ccs = connected_components(g)
		lock(lk) do
        for comp in ccs
			bcount = count(v -> vertex_type(v, n) == 1, comp)
			rcount = count(v -> vertex_type(v, n) == 2, comp)
			((bcount >= kmax) || (rcount >= kmax)) && continue
			
			avg[1+bcount, 1+rcount] += dz
		end
		end
	end
	return avg
end
function avg_weak_cc_dist(t, kernel::AbstractMatrix, n::Vector{Int}; nruns=50, kmax=100)
	# computes the complete average weakly-connected components distribution up to kmax.
	lk = ReentrantLock()

	avg = zeros(kmax)
	Threads.@threads for _ = 1:nruns
		g = generate_multi_graph(t, kernel, n)
		ccs = connected_components(g)
		cc_lengths = map(length, ccs)

		lock(lk) do
			for i = 1:kmax
				avg[i] += count(==(i), cc_lengths)
			end
		end
	end
	return avg ./ nruns
end

avg_cc_dist(t, kernel::AbstractMatrix, N::Int, η::Vector; nruns = 50, kmax=100) = avg_cc_dist(t, kernel, partitionN(N, η); nruns, kmax)
avg_weak_cc_dist(t, kernel::AbstractMatrix, N::Int, η::Vector; nruns = 50, kmax=100) = avg_weak_cc_dist(t, kernel, partitionN(N, η); nruns, kmax)




