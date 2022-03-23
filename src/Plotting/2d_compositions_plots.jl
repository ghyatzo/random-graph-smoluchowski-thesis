using JLD2
using CairoMakie, GraphMakie, NetworkLayout
CairoMakie.activate!(type = "pdf")

includet("../utils.jl")
includet("../RandomGraphs/random_graph_multicomp.jl")


## ------- Components Composition Distribution Study

# numerical component composition distribution through time
nCCE_data = load("../data/comp_densities_evo_multicomp_kmax=99_numerical.jld2")
kernel = nCCE_data["kernel"]
NCCE = nCCE_data["evo"]
T = nCCE_data["tspan"]

# random component composition distribution through time
rCCE_data = load("../data/comp_densities_evo_multicomp_kmax=99_N=1e6_rgraph.jld2")
RCCE = rCCE_data["evo"]
kernel = rCCE_data["kernel"]
N = rCCE_data["N"]
η = rCCE_data["η"]
T = rCCE_data["tspan"]

F_compcomp = Figure(resolution=(1500, 600))
Axncomp = Axis(F_compcomp[1,1], title="ODEs", xlabel=L"k_1", ylabel=L"k_2", ylabelsize=20, xlabelsize=20)
Axrcomp = Axis(F_compcomp[1,2], title="Random Graphs", xlabel=L"k_1", xlabelsize=20)

nhm = heatmap!(Axncomp, log10.(NCCE[:,:,53]), colormap=:roma)
heatmap!(Axrcomp, log10.(RCCE[:,:,53]), colormap=:roma)
Colorbar(F_compcomp[2,1:2], nhm, vertical = false, flipaxis = false,label=L"\log(w_\mathbf{k})", labelsize=20)

save("../imgs/comp_comp/composition_2d_distribution.pdf", F_compcomp)

# Generate skewed distributions, skewed kernels and/or type distribution
unbal_rCCE_data = load("../data/comp_densities_evo_skewed=70-30_kmax=99_N=1e6_rgraph.jld2")
unbal_rCCE = unbal_rCCE_data["evo"]
kernel = unbal_rCCE_data["kernel"]
N = unbal_rCCE_data["N"]
unbal_η = unbal_rCCE_data["η"]
T = unbal_rCCE_data["tspan"]

# very little cross interaction
lowcross_rCCE_data = load("../data/comp_densities_evo_lowcross_kmax=99_N=1e6_rgraph.jld2")
lowcross_rCCE = lowcross_rCCE_data["evo"]
kernel = lowcross_rCCE_data["kernel"]
N = lowcross_rCCE_data["N"]
lowcross_η = lowcross_rCCE_data["η"]
T = lowcross_rCCE_data["tspan"]

# not uniform diagonal entries
offdiag_rCCE_data = load("../data/comp_densities_evo_offdiag_kmax=99_N=1e6_rgraph.jld2")
offdiag_rCCE = offdiag_rCCE_data["evo"]
kernel = offdiag_rCCE_data["kernel"]
N = offdiag_rCCE_data["N"]
offdiag_η = offdiag_rCCE_data["η"]
T = offdiag_rCCE_data["tspan"]

F_compcomp_exp = Figure(resolution=(1800, 600))
Axunbal = Axis(F_compcomp_exp[1,1], title="skewed distribution", xlabel=L"k_1", ylabel=L"k_2", ylabelsize=20, xlabelsize=20)
Axlowx = Axis(F_compcomp_exp[1,2], title="low off diagonal entries", xlabel=L"k_1", xlabelsize=20)
Axdiag = Axis(F_compcomp_exp[1,3], title="unbalanced diagonals", xlabel=L"k_1", xlabelsize=20)

nhm = heatmap!(Axunbal, log10.(unbal_rCCE[:,:,53]), colormap=:roma)
heatmap!(Axlowx, log10.(lowcross_rCCE[:,:,53]), colormap=:roma)
heatmap!(Axdiag, log10.(offdiag_rCCE[:,:,53]), colormap=:roma)
Colorbar(F_compcomp_exp[2,1:3], nhm, vertical = false, flipaxis = false,label=L"\log(w_\mathbf{k})", labelsize=20)

save("../imgs/comp_comp/composition_experiements.pdf", F_compcomp_exp)