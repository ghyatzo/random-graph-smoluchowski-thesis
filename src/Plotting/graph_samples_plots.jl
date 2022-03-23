using JLD2
using CairoMakie, GraphMakie, NetworkLayout
CairoMakie.activate!(type = "pdf")

includet("../utils.jl")
includet("../RandomGraphs/random_graph_multicomp.jl")

# Standard
N = 200
n = partitionN(N, η)
A1 = [1 0; 0 1]
A2 = [0 1; 1 0]
A3 = [0.5 0.5; 0.5 0.5]
tcAn = 1 #normalised critical time.
node_colors = [fill(:darkblue, n[1])..., fill(:darkorange, n[2])...]

begin
    G = generate_multi_graph(tcAn, A1, n)
    f, ax, p = graphplot(G; layout=Spring(C=2.8), node_size=13, edge_width=1, node_color=node_colors)
    hidedecorations!(ax); hidespines!(ax)
    save("../../imgs/graph_demo/2d/separate.pdf", f, pt_per_unit=2)
end

# Inbalanced population
η_inb = [0.7, 0.3]
n_inb = partitionN(N, η_inb)
node_colors_inb = [fill(:black, n_inb[1])..., fill(:red, n_inb[2])...]
begin
    G_inb = generate_multi_graph(tcAn, A2, n_inb)
    f, ax, p = graphplot(G_inb; layout=Spring(C=1.5), node_size=6, edge_width=0.7, node_color=node_colors_inb)
    hidedecorations!(ax); hidespines!(ax)
    save("../../imgs/graph_demo/inbalanced2d.pdf", f, pt_per_unit=2)
end

# tricolored graphs, standard.
Atri = [
    0 0 1;
    0 1 0;
    1 0 0]
tctri = (inv ∘ maximum ∘ eigvals)(Atri)

η_tri = [inv(3), inv(3), inv(3)]
n_tri = partitionN(N, η_tri)
node_colors_tri = [fill(:darkblue, n_tri[1])..., fill(:darkorange, n_tri[2])..., fill(:darkcyan, n_tri[3])...]
begin
        G_tri = generate_multi_graph(1.0, Atri, n_tri)
    f, ax, p = graphplot(G_tri; layout=Spring(C=2.5), node_size=13, edge_width=0.9, node_color=node_colors_tri)
    hidedecorations!(ax); hidespines!(ax)
    save("../../imgs/graph_demo/3d/red-red-alt-blue-cyan.pdf", f, pt_per_unit=2)
end