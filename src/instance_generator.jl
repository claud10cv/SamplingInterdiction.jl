using LightGraphs
# using LightGraphsMatching
using Distributions

function knapsack_random_generate(n, pL, pF, Ef = 0.2)
    d = Normal(100, 30)
    D = 100
    C = 1000
    E = round(Int64, Ef * C)
    vals = rand(1 : C, n)
    fcons = max.(ceil.(Int64, rand(d, n)), 1)
    lcons = max.(ceil.(Int64, rand(d, n)), 1)
    attvals = rand(1 : E, n)
    for i in 1 : n
	attvals[i] = min(vals[i], attvals[i])
    end
    avgdem = 100
    cap = floor(Int64, pF * avgdem)
    budget = floor(Int64, pL * avgdem)
    vals, attvals, lcons, fcons, cap, budget
end

function knapsack_random_generate_martello(n, id)
    vals = rand(1: 100, n)
    followercons = rand(1 : 100, n)
    cap = ceil(Int64, (id / 11) * sum(followercons))
    leadercons = rand(1 : 100, n)
    budget = rand(cap - 10 : cap + 10)
    vals, copy(vals), leadercons, followercons, cap, budget
end

function knapsack_random_generate_fullinterdict(n, d, C)
    vals = rand(1 : C, n)
    demands = rand(1 : d, n)
    att_vals = copy(vals)
    vals, att_vals, demands
end

function bipartite_matching_generate(n, density = 1)
    allarcs = [(i, j) for i in 1 : n for j in 1 : n]
    shuffle!(allarcs)
    target = max(1, floor(Int64, length(allarcs) * density))
    tails = [allarcs[i][1] for i in 1 : target]
    heads = [allarcs[i][2] for i in 1 : target]
    weights = rand(1 : 1000, target)
    attweights = [rand(1 : weights[i]) for i in 1 : target]
    n, tails, heads, weights, attweights
end

function cflp_random_generate_contardo(nfacs, ncusts, targetfacs = 2)
    dems = rand(1 : 100, ncusts)
    sumdems = sum(dems)
    usedfacs = targetfacs
    demperfac = ceil(Int64, sumdems / usedfacs)

    caps = ceil.(Int64, [demperfac * (0.9 + 0.2 * rand()) for f in 1 : nfacs])
    fcosts = rand(8000 : 10000, nfacs)
    attfcosts = rand(1000 : 2000, nfacs)
    pos = rand(0 : 1000, nfacs + ncusts, 2)
    asscosts = zeros(Int64, nfacs, ncusts)
    for i in 1 : nfacs, j in 1 : ncusts
        dij = sqrt((pos[i, 1] - pos[j + nfacs, 1]) * (pos[i, 1] - pos[j + nfacs, 1]) + (pos[i, 2] - pos[j + nfacs, 2]) * (pos[i, 2] - pos[j + nfacs, 2]))
        dij = ceil(Int64, dij * (0.9 + 0.2 * rand()))
        asscosts[i, j] = dij
    end
    nfacs, ncusts, dems, caps, fcosts, asscosts, attfcosts
end

function perfect_matching_random_generate(n, p, C, D)
    mindeg = 5
    deg = zeros(Int64, n)
    alledges = [(i, j) for i in 1 : n - 1 for j in i + 1 : n]
    maxedges = length(alledges)
    nedges = ceil(Int64, p * maxedges)
    shuffle!(alledges)
    tails = []
    heads = []
    numedges = 0
    while numedges < nedges || minimum(deg) < mindeg
        numedges += 1
        t = alledges[numedges][1]
        h = alledges[numedges][2]
        push!(tails, t)
        push!(heads, h)
        deg[t] += 1
        deg[h] += 1
    end
    weights = rand(1 : C, numedges)
    weightsp = rand(1 : D, numedges)
    tails, heads, weights, weightsp
end


function sp_grid_random_generate(n, C, D)
    s = 1
    t = n * n + 2
    tails = []
    heads = []
    function nodeid(r, c)
        (r - 1) * n + c + 1
    end
    for r in 1 : n
        push!(tails, s)
        push!(heads, nodeid(r, 1))
        push!(tails, nodeid(r, n))
        push!(heads, t)
    end
    for r in 1 : n
        for c in 1 : n
            u = nodeid(r, c)
            vvec = []
            if r < n
                push!(vvec, nodeid(r + 1, c))
                if c < n
                    push!(vvec, nodeid(r + 1, c + 1))
                end
            end
            if r > 1
                push!(vvec, nodeid(r - 1, c))
                if c < n
                    push!(vvec, nodeid(r - 1, c + 1))
                end
            end
            if c < n
                push!(vvec, nodeid(r, c + 1))
            end
            for v in vvec
                push!(tails, u)
                push!(heads, v)
            end
        end
    end
    nedges = length(tails)
    weights = rand(1 : C, nedges)
    weightsp = rand(1 : D, nedges)
    tails, heads, weights, weightsp, s, t
end

function mst_grid_random_generate(n, C, D)
    tails = []
    heads = []
    function nodeid(r, c)
        (r - 1) * n + c
    end
    for r in 1 : n
        for c in 1 : n
            u = nodeid(r, c)
            vvec = []
            if r < n
                push!(vvec, nodeid(r + 1, c))
                if c < n
                    push!(vvec, nodeid(r + 1, c + 1))
                end
            end
            if r > 1
                push!(vvec, nodeid(r - 1, c))
                if c < n
                    push!(vvec, nodeid(r - 1, c + 1))
                end
            end
            if c < n
                push!(vvec, nodeid(r, c + 1))
            end
            for v in vvec
                if v > u
                    push!(tails, u)
                    push!(heads, v)
                end
            end
        end
    end
    nedges = length(tails)
    weights = rand(1 : C, nedges)
    weightsp = rand(1 : D, nedges)
    tails, heads, weights, weightsp
end
function mst_random_generate(n, p, mindeg)
    while true
        xpos = rand(1:1000, n)
        ypos = rand(1:1000, n)

        dists = Array{Tuple{Int64, Int64, Int64}, 1}()
        includes = trues(n, n)
        dmat = zeros(Int64, n, n)
        deg = (n - 1) * ones(Int64, n)
        for i in 1 : n, j in i + 1 : n
            d = eucd(xpos[i], ypos[i], xpos[j], ypos[j])
            dmat[i, j] = ceil(Int64, d * rand(0.9 : 1.1))
            push!(dists, (i, j, d))
        end

        sort!(dists, lt = (x, y) -> x[3] > y[3])
        todel = round(Int64, (1 - p) * length(dists))
        currdel = 0
        for (i, j, d) in dists
            if currdel >= todel
                break
            elseif min(deg[i], deg[j]) <= mindeg
                continue
            else
                includes[i, j] = false
                deg[i] -= 1
                deg[j] -= 1
                currdel += 1
            end
        end

        empty!(dists)
        for i in 1 : n, j in i + 1 : n
            if includes[i, j]
                push!(dists, (i, j, dmat[i, j]))
            end
        end

        tails = Int64[]
        heads = Int64[]
        weights = Int64[]
        for (i, j, c) in dists
            push!(tails, i)
            push!(heads, j)
            push!(weights, c)
        end
        # nedges = length(tails)
        # weights = rand(1 : C, nedges)
        # weightsp = rand(1 : D, nedges)
        mst, edges = MSTInterdict.mst(tails, heads, ones(length(tails)))
        if length(edges) >= n - 1
            maxw = round(Int64, maximum(weights) / 10)
            nedges = length(tails)
            weightsp = rand(1 : maxw, nedges)
            return tails, heads, weights, weightsp
        end
    end
end

function eucd(x1, y1, x2, y2)
    return max(1, ceil(Int64, sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2))))
end
