using OffsetArrays

# return true, every k iterations. utility to filter iteration messages.
every(k, n) = n % k == 0

# swaps values below eps to a zero.
sparsify!(a; eps=1e-12) = a[abs.(a) .< eps] .= zero(eltype(a));

# rise all zero to eps
zerolift(a; eps=1e-12) = map(a) do e; return e < eps ? eps : e end

# all values below eps become NaN (useful for logarithms)
nanmap(a; eps=1e-12) = map(a) do e; return e < eps ? NaN : e end

# returns k*w[k]
massmultiply(a) = map(eachindex(a)) do k; k*a[k] end


function partitionN(N, η::Vector)
	@assert sum(η) == 1 "the sum of the parts must sum to 1"
	n = @. round(Int, η * N)
	return n
end


vertex_type(v::Int, n::Vector{Int}) = findfirst(>=(v),cumsum(n))


# 1 Dimension Utilities

# normalise the distribution wrt to the mass PDF
mass_norm!(a::Vector) = a ./= mass_norm(a)
mass_norm(a::Vector) = mapreduce(k->k*a[k], +, eachindex(a))

# exact solutions
_w_exact_factorial(k, t) = k^(k-2)*t^(k-1)*exp(-k*t)/factorial(k) # k must be lower than 20.
_w_stirling(k, t) = inv(2sqrt(π))*t^(k-1)*exp(k*(1-t))/(k^(5/2))

w_exact_mono(k, t) = k < 20 ? _w_exact_factorial(k,t) : _w_stirling(k,t)
u_exact_mono(t, z, nterms) = mapreduce(k -> k*w_exact_mono(k,t)*exp(-k*z), +, 1:nterms)
