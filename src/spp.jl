using LightGraphs
using DataStructures
using LinearAlgebra

function spp(tails, heads, weights, demands, adjlist, s, t)
    nodes = vcat(tails, heads)
    weights = convert(Array{Float64, 1}, weights)
    sort!(nodes)
    unique!(nodes)
    nnodes = length(nodes)
    nedges = length(tails)
    full_nnodes = max(maximum(tails), maximum(heads))
    index = nnodes > 0.75 * full_nnodes ? zeros(Int64, full_nnodes) : spzeros(Int64, full_nnodes)
    for (i, v) in enumerate(nodes)
        index[v] = i
    end
    costmx = []
    incmx = []
    let
        I = [index[tails[i]] for i in 1 : nedges]
        J = [index[heads[i]] for i in 1 : nedges]
        costmx = sparse(I, J, weights, nnodes, nnodes)#[]
        incmx = sparse(I, J, collect(1 : nedges), nnodes, nnodes)
    end
    g = DiGraph(nnodes)
    for e in 1 : nedges
        i, j = index[tails[e]], index[heads[e]]
        add_edge!(g, i, j)
    end
    dstate = dijkstra_shortest_paths(g, index[s], costmx)
    optedges = []
    r = index[t]
    while r != index[s]
        pr = dstate.parents[r]
        push!(optedges, incmx[pr, r])
        r = pr
    end
    optweight = dstate.dists[index[t]]
    if abs(optweight - round(Int64, optweight)) < 1e-5
        optweight = round(Int64, optweight)
    end
    sort!(optedges)
    optweight, optedges, []#parents
end

function spp_st(s, t)::Function
    fixed_spp(tails, heads, weights, demands, adjlist, att) = spp(tails, heads, weights, demands, adjlist, s, t)
    return fixed_spp
end
