using Graphs

function random_component_size_distribution(c, N)
    p = c*inv(N)
    ne::Integer = c*(N-1)÷2

    # Generate an Erdos-Renyi Graph
    g = Graphs.erdos_renyi(N, ne)

    # Count the connected Components
    cc = Graphs.connected_components(g)

    # sort them by length
    cc_lengths = sort(length.(cc))

    # count how many connected components there are of size k
    dist = map(n -> count(==(n), cc_lengths), 1:maximum(cc_lengths))

    # normalise by the total mass of the system
    return dist ./ N
end