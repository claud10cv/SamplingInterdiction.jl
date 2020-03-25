function mst(tails, heads, weights, demands, adjlist)
    return mst(tails, heads, weights)
end

function mst(tails, heads, weights)
    nnodes = max(maximum(tails), maximum(heads))
    nedges = length(tails)
    component = collect(1 : nnodes)
    edges = [(tails[i], heads[i], weights[i], i) for i in 1 : nedges]
    sort!(edges, lt = (x, y) -> x[3] < y[3])

    totalWeight = 0
    optEdges = []
    for e in 1 : nedges
        t, h, w, i = edges[e]
        ct = component[t]
        ch = component[h]
        if ct == ch
            continue
        end
        for u in 1 : nnodes
            if component[u] == ct || component[u] == ch
                component[u] = min(ct, ch)
            end
        end
        totalWeight += w
        push!(optEdges, i)
    end
    sort!(optEdges)
    parents = zeros(Int64, nnodes)
    parents[1] = 1
    numparents = 1
    while numparents < nnodes
        for e in optEdges
            t, h = tails[e], heads[e]
            if parents[t] > 0 && parents[h] <= 0
                parents[h] = t
                numparents += 1
            elseif parents[h] > 0 && parents[t] <= 0
                parents[t] = h
                numparents += 1
            end
        end
    end

    totalWeight, optEdges, parents
end
