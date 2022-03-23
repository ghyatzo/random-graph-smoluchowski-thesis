# using OrdinaryDiffEq, DiffEqOperators
using LinearAlgebra
using BenchmarkTools
using NLsolve
using GLMakie
using LazyGrids
using ImageFiltering
using StaticArrays
using JLD2
using Base.Threads

# Burgers inviscid equation:
# uₜ = -uₓu - uₓ
# uₜ = -uₓ(u - M)
# u = u(x, t)

## discretise the function in the domain: z ∈ [0,1]
# FINITE DIFFERENCE METHOD
Δz = 0.01
Δt = 0.01
z = 0:Δz:5
tspan = 0.0:Δt:1.0

∂z(u, i, Δz) = (u[i+1]-u[i-1])*inv(2Δz)
stepfunc(Δz, Δt) = u -> begin
    Du = ∂z.((BorderArray(u, Pad(1)),), eachindex(u), Δz)
    return u .- (Du .* (u .- 1)).*Δt
end
residuals(x, b) = stepfunc(Δz, Δt)(x) .- b
solstep(u0) = (U, x) -> begin
    U .= residuals(x, u0)
end

sol = zeros(length(z))
u = exp.(-z)
sol .= u
jldopen("1dburgers.jld2", "a+") do file
    file["sol/1"] = sol
end
jldopen("1dburgers.jld2", "r+") do file
    for i in 2:length(tspan)
        @time sol .= nlsolve(solstep(u), u).zero
        u .= u1
        file["sol/$i"] = sol
    end
end

surface(z, tspan, sol, axis=(type=Axis3,))
wireframe!(z[1:2:end], collect(tspan)[1:5:end], sol[1:2:end, 1:5:end], axis=(type=Axis3,), color=:black)


# n = length(u0)
# W = Array{Float64, 2}(undef, n, size(sol)[2])
# A = [k*exp(-k*z) for z in z, k in 1:n]
# Adag = adag(A, 0.001)
# Adag*A
# for t in axes(sol, 2)
#     W[:, t] .= Adag*sol[:, t]
# end

# W
# surface(z, tspan, W, axis=(type=Axis3,))

# lines(tspan, W[1, :])
# lines!(tspan, W[2, :])
# lines!(tspan, W[3, :])
# lines!(tspan, W[4, :])
# lines!(tspan, W[5, :])
# lines!(tspan, W[6, :])

∂x(u, comp, I::CartesianIndex, Δx) = begin
    i, j = Tuple(I)
    return (u[comp, i-1, j] - u[comp, i+1, j])/(2Δx)
end
∂y(u, comp, I::CartesianIndex, Δy) = begin
    i, j = Tuple(I)
    return (u[comp, i, j-1] - u[comp, i, j+1])/(2Δy)
end


Δx = Δy = 0.05
Δt = 0.001
tspan = 0.0:Δt:1.0
x = Δx:Δx:2
y = Δy:Δy:2

xg, yg = ndgrid(x, y)
Ju = @MArray zeros(2, 2)
u = zeros(2, length(x), length(y))
u[1, :, :] .= exp.(-xg)
u[2, :, :] .= exp.(-yg)
sol = zeros(2, length(x), length(y))
sol[1, :, :] .= u[1, :, :]
sol[2, :, :] .= u[2, :, :]
# jldopen("2dburgers.jld2", "w") do file
#     file["sol/1"] = sol
# end

stepfunc2d!(Ju, xg, Δx, Δy, Δt) = (R, u) -> begin
    au = BorderArray(u, Pad(0, 1, 1))
    @inbounds @threads for Id in eachindex(xg)
        Ju[1, 1] = ∂x(au, 1, Id, Δx)
        Ju[1, 2] = ∂y(au, 1, Id, Δy)
        Ju[2, 1] = ∂x(au, 2, Id, Δx)
        Ju[2, 2] = ∂y(au, 2, Id, Δy)
        @views R[:, Id] .= u[:, Id] .+ (Ju*(u[:, Id] .- 1)).*Δt
    end
end

F! = stepfunc2d!(Ju, xg, Δx, Δy, Δt)
solstep2d(u) = (R, x) -> begin
    F!(R, x)
    R .-= u
end

nlsolve(solstep2d(u), u; show_trace=true, xtol=1e-8)

for i in 2:length(tspan)
    @show tspan[i]
    @time sol .= nlsolve(solstep2d(u), u).zero
    # file["sol/$i"] = sol
    u = sol
end


jldopen("2dburgers.jld2", "r+") do file
    for i in 2:length(tspan)
        @show tspan[i]
        @time sol .= nlsolve(solstep2d(u), u).zero
        file["sol/$i"] = sol
        u = sol
    end
end


sol = zeros(2, length(x), length(y), length(tspan))
jldopen("2dburgers.jld2", "r") do file
    for i in 1:length(tspan)
        sol[:, :, :, i] = file["sol/$i"]
    end
end

surface(x, y, sol[1, :, :, end])
surface!(x, y, sol[2, :, :, end])

u1 = sol[1, :, :, :]
u2 = sol[2, :, :, :]
norm = map(x -> sqrt(x[1]^2 + x[2]^2), zip(u1, u2))

z = Node(norm[:, :, 1])
u1 = Node(sol[1, :, :, 1])
u2 = Node(sol[2, :, :, 1])
fig, ax, surf = wireframe(x, y, z)
surface!(x, y, u2)

n = length(tspan)
i = 1
direction = 1
while true
    if i < 1
        direction = 1
        i = 2
    elseif i > n
        direction = -1
        i = n-1
    end
    # u1[] = sol[1, :, :, i]
    # u2[] = sol[2, :, :, i]
    z[] = norm[:, :, i]
    sleep(0.1)
    i = i + direction
end
