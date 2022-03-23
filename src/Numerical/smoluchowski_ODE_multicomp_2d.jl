includet("../utils.jl")

using DifferentialEquations
using LinearAlgebra
using OffsetArrays
using SparseArrays
using Base: Threads

# \sum_{x\in\Omega_k} w_x w_{k-x} x^TA(k-x)
function source_multiconv(A::Matrix, w::AbstractMatrix, k::Vector)
	# Only two dimensional, first term in the smoluchowski equations
	x = zeros(Int, length(k))
	y = zeros(Int, length(k))

	k1, k2 = k
	res = 0
	@inbounds for i in 0:k1
		x[1] = i
		y[1] = k1 - i
		@inbounds for j in 0:k2
			(i == j == 0) || (i == k1 && j == k2) && continue
			x[2] = j
			y[2] = k2 - j
			res += @views w[x[1], x[2]] * w[y[1], y[2]] * dot(x, A, y)
		end
	end

	return res
end

# \sum_{x\in\N^d} w_xw_kx^TAk
function sink_conv_multi(A, w, k::Vector{T}) where T <: Integer
	# only two dimensional, second term in the smoluchowski equations
	res = 0.0
	s = zeros(T, length(k))
	wk = w[k...]
	@inbounds for j in axes(w, 2)
		s[2] = j
		@inbounds for i in axes(w, 1)
			s[1] = i
			res += wk * w[i,j] * dot(s, A, k)
		end
	end
	return res
end

function smoluchowski2d!(du, u, p ,t)
	A, = p
	for j in axes(du, 2)
		for i in axes(du, 1)
			du[i,j] = 0.5*source_multiconv(A, u, [i,j]) - sink_conv_multi(A, u, [i,j])
		end
	end
end

function solve_smoluchowski_2d(A, w0)
	# what does it mean to mass normalise in this case?
	normA = maximum(eigvals(A))
	tc = inv(normA)
	problem = ODEProblem(smoluchowski2d!, w0, (0., tc), (A,))
	return solve(problem, Tsit5(); progress = true, progress_steps = 1)
end

