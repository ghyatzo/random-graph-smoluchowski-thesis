include("../utils.jl")

using DifferentialEquations

# implement the monocomponent smoluchowski System of equations as a system of ODEs
ker(x::Integer, y::Integer, w::Vector) =  begin
	@views x*y*w[x]*w[y]
end

source_conv(ker::Function, w::Vector, k) = mapreduce(s->ker(s, k-s, w), +, 1:(k-1))
sink_conv(ker::Function, w::Vector, k, kmax) = mapreduce(s->ker(s, k, w), +, 1:(kmax))

function smoluchowski!(du, u, p, t)
	ker, kmax = p
	du[1] = -sink_conv(ker, u, 1, kmax)

	@inbounds for k in 2:length(du)
		du[k] = 0.5*source_conv(ker, u, k) - sink_conv(ker, u, k, kmax)
	end
end

function solve_smoluchowski_mono(ker, w0, kmax)
	p = (ker, kmax)

	problem = ODEProblem(smoluchowski!, mass_normalize!(w0), (0., 1.), p)
	return sol = solve(problem, Tsit5())
end
