using LightGraphsMatching
using LightGraphs

function matching(tails, heads, weights, demands, adjlist)
    optval, edges = matching(tails, heads, weights)
    optval, edges, []
end

function matching(tails, heads, weights, incmx)
    nnodes = max(maximum(tails), maximum(heads))
    nedges = length(tails)
    incmx = sparse(tails, heads, collect(1 : nedges), nnodes, nnodes)
    g = Graph(nnodes)
    w = Dict{Edge, Int64}()
    scale = typeof(weights) == Array{Float64, 1} ? 1e+7 : 1
    for e in 1 : nedges
        t, h, u = tails[e], heads[e], weights[e]
        if t > h println("error, $t > $h")
            exit()
        end
        add_edge!(g, t, h)
        w[Edge(t, h)] = round(Int64, u * scale)
    end
    match = minimum_weight_perfect_matching(g, w)
    optweight = match.weight / scale
    if abs(optweight - round(Int64, optweight)) < 1e-5
        optweight = round(Int64, optweight)
    end
#    println(match.mate)
    optedges = [incmx[i, match.mate[i]] for i in 1 : nnodes if i < match.mate[i]]
    sort!(optedges)
    optweight, optedges
end
