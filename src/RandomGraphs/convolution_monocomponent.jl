using OffsetArrays
using PaddedViews
using FastTransforms

# poisson testing, first 100 terms
nthconv(n,c) = [exp(-n*c)*(n*c)^k*inv(factorial(big(k))) for k in 0:100]

# the degree distribution at each t.
function component_size_distribution(u::Vector)
	N = length(u)

	Ou = OffsetArray(u, -1) 					 #0 indexing
	Ou1 = OffsetArray(zeros(eltype(u), N-1), -1) #0 indexing

	# excess distribution
	mu = mapreduce(k -> k*Ou[k], +, eachindex(Ou))
	@inbounds for k in eachindex(Ou1)
		Ou1[k] = inv(mu)*(k+1)Ou[k+1]
	end

	# W is defined only for n >= 1, makes sense to start at 1 the array
	w = similar(u)
	w[1] = Ou[0]

	# pile on the convolution at each iteration.
	Ou1conv = deepcopy(Ou1)
	for n in 2:length(w)
		Ou1conv = OffsetArray(slowerFastLinearConvolution(parent(Ou1conv), parent(Ou1)), -1)
		w[n] = mu*inv(n*(n-1))*Ou1conv[n-2]
	end
	return w
end

function fastLinearConvolution(f, g)
	N = length(f)
	M = length(g)

	f_pad = PaddedView(0, f, (N+M-1,))
	g_pad = PaddedView(0, g, (N+M-1,))

	F = plan_fft(f_pad)

	real.(ifft( (F*f_pad).*(F*g_pad)))
end