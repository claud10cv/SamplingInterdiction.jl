using MatrixNetworks

function bip_matching(tails, heads, weights, demands, adjlist)
    optval, edges = bip_matching(tails, heads, -weights)
    -optval, edges, []
end

function bip_matching(tails, heads, weights)
    nnodes = max(maximum(tails), maximum(heads))
    nedges = length(tails)
    W = sparse(tails, heads, weights, nnodes, nnodes)
    match = bipartite_matching(W)
    optweight = match.weight
    if abs(optweight - round(Int64, optweight)) < 1e-5
        optweight = round(Int64, optweight)
    end
    incmx = sparse(tails, heads, collect(1 : nedges))
    (restails, resheads) = edge_list(match)
    optedges = [incmx[restails[i], resheads[i]] for i in eachindex(restails)]
    sort!(optedges)
    optweight, optedges
end
