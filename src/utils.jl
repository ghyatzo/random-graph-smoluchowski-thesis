using OffsetArrays


every(k, n) = n % k == 0

sparsify!(a, eps) = a[abs.(a) .< eps] .= zero(eltype(a));
zerolift(a, eps) = map(a) do e; return e < eps ? eps : e end
nanmap(a; eps=1e-12) = map(a) do e; return e < eps ? NaN : e end
massmultiply(a) = map(eachindex(a)) do k; k*a[k] end


function partitionN(N, η::Vector)
	@assert sum(η) == 1 "the sum of the parts must sum to 1"
	n = @. round(Int, η * N)
	return n
end
vertex_type(v::Int, n::Vector{Int}) = findfirst(>=(v),cumsum(n))


# 1 Dimension Utilities
# builds the first moment generating function
u(z, Ws) = mapreduce(k -> k*Ws[k]*exp(-k*z), +, eachindex(Ws))
# calculates the second moment
m2(Ws) = mapreduce(k->k^2*Ws[k], +, eachindex(Ws))
# critical time monocomponent
t_c(w0) = inv(m2(w0))
# normalise the distribution wrt to the mass PDF
mass_norm!(a::Vector) = a ./= mass_norm(a)
mass_norm(a::Vector) = mapreduce(k->k*a[k], +, eachindex(a))

# 2 dimensions utilities
function weak_components(z::AbstractMatrix)
	@assert size(z, 1) == size(z, 2)
	weak_z = [mapreduce(+, 0:k) do b
		z[1+b, 1+k-b]
	end for k = 1:minimum(size(z))-1]
end
weak_components(z::T) where T <: OffsetMatrix = weak_components(z.parent)

# exact solutions
_w_exact_factorial(k, t) = k^(k-2)*t^(k-1)*exp(-k*t)/factorial(k) # k must be lower than 20.
_w_stirling(k, t) = inv(2sqrt(π))*t^(k-1)*exp(k*(1-t))/(k^(5/2))
w_exact_mono(k, t) = k < 20 ? _w_exact_factorial(k,t) : _w_stirling(k,t)
u_exact_mono(t, z, nterms) = mapreduce(k -> k*w_exact_mono(k,t)*exp(-k*z), +, 1:nterms)