using JLD2
using CairoMakie, GraphMakie, NetworkLayout
CairoMakie.activate!(type = "pdf")

includet("../utils.jl")
includet("../RandomGraphs/random_graph_multicomp.jl")

## ------- Weak Components Study --------------

#Defaults
kernel = [1.0 0.9;0.9 1.0]
N = Int(1e6)
η = [0.5, 0.5]
n = partitionN(N, η)    

kmax = size(w0, 1)
tc = (inv ∘ maximum ∘ eigvals)(kernel)
T = 0:0.01:tc

# numerical weak components distribution through time - reference
 nWCE = load("../../data/weak_comp_evo_multicomp_kmax=99_numerical.jld2")["evo"]
# nWCE = zeros(kmax-1, length(T))
# for (it, t) = enumerate(T)
# 	tmp = mass_norm!(weak_components(num_sol(t)))
# 	sparsify!(tmp, ceps)
# 	nWCE[:, it] = tmp
# end

# random weak components distribution through time
rWCE = load("../../data/weak_comp_evo_multicomp_kmax=99_N=1e6_rgraph.jld2")["evo"]
# rWCE = zeros(kmax-1, length(T))
# for (it, t) = enumerate(T)
# 	rWCE[:, it] = avg_weak_cc_dist(t, kernel, N, η; kmax=kmax-1) ./ N,
# end


## FIGURE 1 - compare ODEsol and RGsol weak component density distribution for 3 times 
nWCDD_start = nWCE[:, 5]
rWCDD_start = rWCE[:, 5]
nWCDD_mid = nWCE[:, 29]
rWCDD_mid = rWCE[:, 29]
nWCDD_end = nWCE[:, end]
rWCDD_end = rWCE[:, end]

ax_kwargs = (yscale=log10,
 yminorticksvisible = true,
 yminorgridvisible = true,
 yminorticks = IntervalsBetween(10),
 limits = (0,100, 8e-9, 2),
 xlabel="k", ylabel=L"\omega_k", aspect=1)

num_lines_kw = (linewidth = 1.5, color=:darkred)
rg_scatter_kw = (marker=:utriangle, color=:darkorange, markerksize=30)
F_WCDD_times = Figure(resolution=(1500, 500))
A_WCDD_start = Axis(F_WCDD_times[1, 1]; title=L"t=0.04", ax_kwargs...)
A_WCDD_mid = Axis(F_WCDD_times[1, 2]; title=L"t=0.28", ax_kwargs...)
A_WCDD_end = Axis(F_WCDD_times[1, 3]; title=L"t=t_c", ax_kwargs...)

lines!(A_WCDD_start, 1..99, nanmap(nWCDD_start); num_lines_kw...)
scatter!(A_WCDD_start, 1..99, nanmap(rWCDD_start); rg_scatter_kw...)

lines!(A_WCDD_mid, 1..99, nanmap(nWCDD_mid); num_lines_kw...)
scatter!(A_WCDD_mid, 1..99, nanmap(rWCDD_mid); rg_scatter_kw...)

lines!(A_WCDD_end, 1..99, nanmap(nWCDD_end); num_lines_kw...)
scatter!(A_WCDD_end, 1..99, nanmap(rWCDD_end); rg_scatter_kw...)

save("../../imgs/weak_cc/ODEvRG_density_3timeevo.pdf", F_WCDD_times)

## FIGURE 2 - compare ODEsol and RGsol weak component density distribution for 3 values of N
# t = T[35]
# Nspan = [10^4, 10^5, 10^6]
# NN = zeros(kmax-1, length(Nspan))
# for (ni, nn) = enumerate(Nspan)
# 	@info "running n: $nn"
# 	test_sol = avg_weak_cc_dist(t, kernel, nn, η; kmax=kmax-1) ./ nn
# 	NN[:, ni] = test_sol
# end
# save("../data/weak_comp_N=1e4-5-6_evo_multicomp_kmax=99_rgraphs.jld2", Dict("evo" => NN, "t" => t, "kernel" => kernel, "Nspan" => Nspan))
NN = load("../data/weak_comp_N=1e4-5-6_evo_multicomp_kmax=99_rgraphs.jld2")["evo"]

num_ref = nWCE[:, 35]
F_WCDD_Ns = Figure(resolution=(1000, 500))
AxWCDD1e4 = Axis(F_WCDD_Ns[1, 1]; title="N = 1e4", ax_kwargs...)
AxWCDD1e5 = Axis(F_WCDD_Ns[1, 2]; title="N = 1e5", ax_kwargs...)
AxWCDD1e6 = Axis(F_WCDD_Ns[1, 2]; title="N = 1e6", ax_kwargs...)

lines!(AxWCDD1e4, 1..99, nanmap(num_ref); num_lines_kw...)
scatter!(AxWCDD1e4, 1..99, nanmap(NN[:, 1]); rg_scatter_kw...)

lines!(AxWCDD1e5, 1..99, nanmap(num_ref); num_lines_kw...)
scatter!(AxWCDD1e5, 1..99, nanmap(NN[:, 2]); rg_scatter_kw...)

lines!(AxWCDD1e6, 1..99, nanmap(num_ref); num_lines_kw...)
scatter!(AxWCDD1e6, 1..99, nanmap(NN[:, 3]); rg_scatter_kw...)

save("../../imgs/weak_cc/ODEvRG_density_2Nevo.pdf", F_WCDD_Ns)