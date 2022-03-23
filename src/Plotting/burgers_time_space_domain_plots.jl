using JLD2
using CairoMakie, GraphMakie, NetworkLayout
CairoMakie.activate!(type = "pdf")

includet("../utils.jl")
includet("../RandomGraphs/random_graph_multicomp.jl")

## --------- 1D - Burgers' Equation Space/Time Domain of the derivative

## spatial/time domain plots look at burgers equation usual rapresentation in 1d.
# we use the weak component data.
# numerical weak components distribution through time - reference
nWCE_data = load("../data/weak_comp_evo_multicomp_kmax=99_numerical.jld2")
nWCE = nWCE_data["evo"]
T = nWCE_data["tspan"]

# random weak components distribution through time
rWCE = load("../data/weak_comp_evo_multicomp_kmax=99_N=1e6_rgraph.jld2")["evo"]

# compute Burgers's derivative
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
function burgers_deriv(u)
    N = length(u)
    dx = 0.01
    du = similar(u)
    @inbounds for i in 1:N
        um1, up1 = bc(u, i-1, N, dx), bc(u, i+1, N, dx)
        du[i] = -inv(2dx)*(up1 - um1)*(u[i] - 1)
    end
    return du
end

#generate burgers' vector
Z = 0:0.01:5
rU(Z, t) = @. u(Z, (rWCE[:, t],))
nU(Z, t) = @. u(Z, (nWCE[:, t],))
nDU = zeros(length(Z), length(T))
rDU = zeros(length(Z), length(T))
for t in eachindex(T)
    nDU[:, t] .= burgers_deriv(nU(Z, t))
    rDU[:, t] .= burgers_deriv(rU(Z, t))
end
F_burger = Figure(resolution=(1800, 1000))
Ax_num = Axis(F_burger[1,1]; title = "Derivative space/time domain - ODE system", xlabel = L"T", ylabel=L"z", xlabelsize=30, ylabelsize=30)
Ax_ran = Axis(F_burger[1,2]; title = "Derivative space/time domain - Random Graphs",xlabel = L"T", xlabelsize=30)
Ax_num_line = Axis(F_burger[2,1]; xlabel = L"z", ylabel =L"-u_t(z)", xlabelsize=30, ylabelsize=30)
Ax_ran_line = Axis(F_burger[2,2]; xlabel = L"z", xlabelsize=30)

nhm = heatmap!(Ax_num, T, Z, abs.(nDU'), colormap=:roma)
rhm = heatmap!(Ax_ran, T, Z, abs.(rDU'), colormap=:roma)
Colorbar(F_burger[1,3], rhm, label = L"-u_t(z)", labelsize=30)
for t in [5, 15, 30, 40, 50]
    lines!(Ax_num_line, Z, abs.(nDU[:, t]), label="t = $(T[t])")
    lines!(Ax_ran_line, Z, abs.(nDU[:, t]), label="t = $(T[t])")
end

axislegend(Ax_num_line, valign=:top, framevisible=false)
axislegend(Ax_ran_line, valign=:top, framevisible=false)

save("../imgs/burgers/space_time_domain_ODEandRG_burgers.pdf", F_burger, pt_per_unit=2)
