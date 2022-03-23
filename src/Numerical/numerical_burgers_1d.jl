include("../utils.jl")

using DifferentialEquations

# Burgers inviscid equation:
# uₜ = -uₓu - uₓ
# uₜ = -uₓ(u - M)
# u = u(x, t)

# linear extrapolation for ghost nodes.
function bc(u, i, N, step)
	if i == 0
		m = inv(step)*(u[2]-u[1])
		return u[1] - m*step
	elseif i == N+1
		m = inv(step)*(u[end]-u[end-1])
		return u[end] + m*step
	else
		return u[i]
	end
end

function burgers_1d!(du, u, p, t)
	M, dx, N = p
	@inbounds for i in 1:N
		um1, up1 = bc(u, i-1, N, dx), bc(u, i+1, N, dx)
		du[i] = -inv(2dx)*(up1 - um1)*(u[i] - M)
	end
end

function solve_burgers_1d(N, grid, u0, tc)
	p = (1.0, step(grid), N)

	problem = ODEProblem(burgers_1d!, mass_normalize!(u0), (0.0, tc), p)
	# Solvers:
	# Tsit5(), Rosenbrock23(), radau(), RadauIIA3(), Rodas5()
	sol = solve(problem, RadauIIA3())

	return sol
end
