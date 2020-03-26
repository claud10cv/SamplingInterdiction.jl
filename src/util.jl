#using LightGraphsFlows

function build_adj_lists(tails, heads, is_symmetric::Bool)
    nnodes = max(maximum(tails), maximum(heads))
    nedges = length(tails)
    adjlist = [[] for i in 1 : nnodes]
    for e in 1 : nedges
        t, h = tails[e], heads[e]
        push!(adjlist[t], (e, h))
        if is_symmetric
            push!(adjlist[h], (e, t))
        end
    end
    adjlist
end

function build_cost_matrix(tails, heads, weights, is_symmetric::Bool)
    nnodes = max(maximum(tails), maximum(heads))
    nedges = length(tails)
    costmx = spzeros(Int64, nnodes, nnodes)
    incmx = spzeros(Int64, nnodes, nnodes)
    for e in 1 : nedges
        t, h, w = tails[e], heads[e], weights[e]
        if is_symmetric
            costmx[t, h] = costmx[h, t] = w
            incmx[t, h] = incmx[h, t] = e
        else
            costmx[t, h] = w
            incmx[t, h] = e
        end
    end
    costmx, incmx
end

function print_sp_instance(tails, heads, weights, attacked, s, t, filename)
    nnodes = max(maximum(tails), maximum(heads))
    nedges = length(tails)
    open(filename, "w") do f
        write(f, "p sp $nnodes $nedges\n")
        for e in 1 : nedges
            write(f, "a $(tails[e] - 1) $(heads[e] - 1) $(weights[e]) $(attacked[e])\n")
        end
    end
end
function print_mst_instance(tails, heads, weights, attacked, filename)
    nnodes = max(maximum(tails), maximum(heads))
    nedges = length(tails)
    open(filename, "w") do f
        write(f, "p mst $nnodes $nedges\n")
        for e in 1 : nedges
            write(f, "e $(tails[e] - 1) $(heads[e] - 1) $(weights[e]) $(attacked[e])\n")
        end
    end
end

function print_knapsack_instance(vals, attacked, leadercons, followercons, cap, budget, filename)
    nnodes = length(vals)
    open(filename, "w") do f
        write(f, "p kp $nnodes\n")
        write(f, "c f $cap\n")
        write(f, "c l $budget\n")
        for n in 1 : nnodes
            write(f, "n $(vals[n]) $(attacked[n]) $(leadercons[n]) $(followercons[n])\n")
        end
    end
end

function print_knapsack_instance(vals, attacked, demands, k, filename)
    nnodes = length(vals)
    open(filename, "w") do f
        write(f, "p kp $nnodes\n")
        write(f, "k $k\n")
        for n in 1 : nnodes
            write(f, "n $(vals[n]) $(attacked[n]) $(demands[n])\n")
        end
    end
end

function print_perfect_matching_instance(tails, heads, weights, attacked, filename)
    nnodes = max(maximum(tails), maximum(heads))
    nedges = length(tails)
    open(filename, "w") do f
        write(f, "p pm $nnodes $nedges\n")
        for e in 1 : nedges
            write(f, "e $(tails[e] - 1) $(heads[e] - 1) $(weights[e]) $(attacked[e])\n")
        end
    end
end

function print_cflp_instance(nfacs, ncusts, dems, caps, fcosts, asscosts, attfcosts, filename)
    open(filename, "w") do f
        println(f, "p cflp $nfacs $ncusts")
        for i in 1 : nfacs
            println(f, "f $(caps[i]) $(fcosts[i]) $(attfcosts[i])")
        end
        for i in 1 : ncusts
            println(f, "n $(dems[i])")
        end
        for i in 1 : nfacs, j in 1 : ncusts
            println(f, "a $i $j $(asscosts[i, j])")
        end
    end
end


function batch_cflp_random_generator()
#    N = [500, 1000]
#    F = [50, 100]
    F = [50]
    N = [50]
    TF = [5]
    S = collect(1 : 9)
    for n in N, f in F, tf in TF, s in S
        filepath = "./instances/flp/cflp/contardo-sefair"
        filename = "$filepath/Rand$(n)_$(f)_$(tf)_$(s).txt"
        nfacs, ncusts, dems, caps, fcosts, asscosts, attfcosts = cflp_random_generate_contardo(f, n, tf)
        print_cflp_instance(nfacs, ncusts, dems, caps, fcosts, asscosts, attfcosts, filename)
    end
end

function batch_perfect_matching_random_generator()
    N = [50, 100, 200, 500, 1000, 2000]
    CD = [(10, 5), (10, 10), (10, 20), (100, 50), (100, 100), (100, 200)]
    P = [0.1, 0.25, 0.5, 0.75]
    for n in N, (c, d) in CD, p in P, i in 0 : 9
        print("printing instance for n = $n, (c, d) = ($c, $d), p = $p and i = $i...")
        tails, heads, weights, weightsp = perfect_matching_random_generate(n, p, c, d)
        filename = "Rand$(n)_$(c)-$(d)_$(round(Int64, 100 * p))_$(i).txt"
        filepath = "./instances/perfectmatching/"
        print_perfect_matching_instance(tails, heads, weights, weightsp, "$filepath/$filename")
        println("done")
    end
end

function batch_spp_grid_generator()
    N = [100, 200, 300, 400, 500]
    CD = [(10, 5), (10, 10), (10, 20), (100, 50), (100, 100), (100, 200)]

    for n in N, (c, d) in CD
        for i in 0 : 9
            print("printing instance for n = $n, (c, d) = ($c, $d) and i = $i...")
            tails, heads, weights, weightsp, s, t = sp_grid_random_generate(n, c, d)
            filename = "Grid$(n)x$(n)_$(c)-$(d)_$(i).txt"
            filepath = "./instances/shortestpath/grotescas"
            print_sp_instance(tails, heads, weights, weightsp, s, t, "$filepath/$filename")
            println("done")
        end
    end
end

function batch_knapsack_random_generator()
  #  N = [10000, 100000, 1000000]
  #  N = [1000000]
  #  N = [100, 1000]
    N = [50]
    PL = [5, 10, 20]
    PF = [20, 40]
    S = collect(1 : 9)
    EF = [0.125, 0.25]
    for n in N, pl in PL, pf in PF, s in S, ef in EF
        print("printing instance for n = $n, pl = $pl, pf = $pf, s = $s...")
	intef = round(Int64, 1000 * ef)
        vals, attacked, leadercons, followercons, cap, budget = knapsack_random_generate(n, pl, pf, ef)
        filename = "Rand$(n)_$(pl)_$(pf)_$(intef)_$(s).txt"
        filepath = "./instances/knapsack/contardo-sefair"
        print_knapsack_instance(vals, attacked, leadercons, followercons, cap, budget, "$filepath/$filename")
        println("done")
    end
end

function batch_mst_random_generator()
    N = [100, 200, 500, 1000, 5000]
    P = [0.1, 0.25, 0.5, 0.75]
    mindeg = [10]
    for n in N, p in P, m in mindeg
        for i in 0 : 9
            print("printing instance for n = $n, p = $p, m = $m and i = $i...")
            tails, heads, weights, weightsp = mst_random_generate(n, p, m)
            filename = "Rand$(n)_$(round(Int64, 100 * p))_$(m)_$(i).txt"
            filepath = "./instances/mst"
            print_mst_instance(tails, heads, weights, weightsp, "$filepath/$filename")
            println("done")
        end
    end
end

function batch_mst_grid_generator()
    N = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    CD = [(10, 5), (10, 10), (10, 20), (100, 50), (100, 100), (100, 200)]

    for n in N, (c, d) in CD
        for i in 0 : 9
            print("printing instance for n = $n, (c, d) = ($c, $d) and i = $i...")
            tails, heads, weights, weightsp = mst_grid_random_generate(n, c, d)
            filename = "Grid$(n)x$(n)_$(c)-$(d)_$(i).txt"
            filepath = "./instances/mst"
            print_mst_instance(tails, heads, weights, weightsp, "$filepath/$filename")
            println("done")
        end
    end
end

function compute_edge_lower_bounds_spp(tails, heads, weights)
    nnodes = max(maximum(tails), maximum(heads))
    nedges = length(tails)
    s = 1
    t = nnodes
    adjlists_fw = build_adj_lists(tails, heads, false)
    costmx_fw, incmx_fw = build_cost_matrix(tails, heads, weights, false)
    adjlists_bw = build_adj_lists(heads, tails, false)
    costmx_bw, incmx_bw = build_cost_matrix(heads, tails, weights, false)
    fwg = DiGraph(nnodes)
    bwg = DiGraph(nnodes)
    for e in 1 : nedges
        i, j = tails[e], heads[e]
        add_edge!(fwg, i, j)
        add_edge!(bwg, j, i)
    end
    dstate_fw = dijkstra_shortest_paths(fwg, s, costmx_fw)
    dstate_bw = dijkstra_shortest_paths(bwg, t, costmx_bw)
    edge_lb = zeros(Int64, nedges)
    for e in 1 : nedges
        spp_fw = dstate_fw.dists[tails[e]]
        spp_bw = dstate_bw.dists[heads[e]]
        cij = weights[e]
        lb = spp_fw + cij + spp_bw
        edge_lb[e] = lb
    end
    edge_lb
end

function find_attack(att, others)
    att1 = Set(att.attacked)
    for a in others
        att0 = Set(a.attacked)
        if att0 == att1
            return true
        end
    end
    return false
end

function find_tree(t, others)
    t1 = Set(t.edges)
    for s in others
        s1 = Set(s.edges)
        if s1 == t1
            return true
        end
    end
    return false
end

function draw_grid(tails, heads, weights, attacked_weights)
    nnodes = max(maximum(tails), maximum(heads))
    gsize = round(Int64, sqrt(nnodes - 2))
    half = round(Int64, gsize / 2)
    nedges = length(tails)
    println("$nnodes $nedges $gsize $half")
    f = open("grid_example.tex", "w")
    println(f, "\\documentclass{minimal}")
    println(f, "\\usepackage{tikz}")
    println(f, "\\usepackage{graphicx}")
    println(f, "\\begin{document}")
    println(f, "%\\begin{figure}")
    println(f, "\\begin{tikzpicture}[darkstyle/.style={circle,draw,fill=gray!40,minimum size=2},")
    println(f, "\t\tendpoint/.style={circle,draw,fill=red!40,minimum size=5}]")
    println(f, "\t\\node [endpoint] (s) at (-2, $half) {s};")
    println(f, "\t\\node [endpoint] (t) at ($(gsize + 1), $(half)) {t};")
    for n in 0 : nnodes - 3
        i = floor(Int64, n % gsize)
        j = floor(Int64, n / gsize)
        println(f, "\t\\node [darkstyle] ($n) at ($j, $i) {};")
    end
    for e in 1 : nedges
        t, h, w, aw = tails[e], heads[e], weights[e], attacked_weights[e]
        if t == 1
            t = "s"
        elseif t == nnodes
            t = "t"
        else
            t -= 2
        end
        if h == 1
            h = "s"
        elseif h == nnodes
            h = "t"
        else
            h -= 2
        end
        println(f, "\\draw [->] ($t) -- ($h);")
    end
    println(f, "\\end{tikzpicture}")
    println(f, "%\\caption{Original Grid network, \$10\\times 10\$}")
    println(f, "%\\end{figure}")
    println(f, "\\end{document}")
    close(f)
end

function draw_algorithm(tails, heads, paths, attacks)
    nnodes = max(maximum(tails), maximum(heads))
    gsize = round(Int64, sqrt(nnodes - 2))
    half = round(Int64, gsize / 2)
    nedges = length(tails)
    println("$nnodes $nedges $gsize $half")
    f = open("grid_example_alg.tex", "w")
    println(f, "\\documentclass{minimal}")
    println(f, "\\usepackage{tikz}")
    println(f, "\\usepackage{graphicx}")
    println(f, "\\begin{document}")
    println(f, "%\\begin{figure}")
    println(f, "\\begin{tikzpicture}[darkstyle/.style={circle,draw,fill=gray!40,minimum size=2},")
    println(f, "\t\tendpoint/.style={circle,draw,fill=red!40,minimum size=5},")
    println(f, "\t\tspp/.style={red},")
    println(f, "\t\treg/.style={black},")
    println(f, "\t\tatt/.style={blue}]")
    println(f, "\t\\node [draw, endpoint] (s) at (-2, $half) {s};")
    println(f, "\t\\node [draw, endpoint] (t) at ($(gsize + 1), $(half)) {t};")
    for n in 0 : nnodes - 3
        i = floor(Int64, n % gsize)
        j = floor(Int64, n / gsize)
        println(f, "\t\\node [draw, darkstyle] ($n) at ($j, $i) {};")
    end
    maxits = length(attacks)
    slideno = 1
    for nit in 1 : maxits
        for p in 1 : nit
            if p < nit
                style = "reg"
            else
                style = "spp"
            end
            for e in paths[p]
                t, h = tails[e], heads[e]
                if t == 1
                    t = "s"
                elseif t == nnodes
                    t = "t"
                else
                    t -= 2
                end
                if h == 1
                    h = "s"
                elseif h == nnodes
                    h = "t"
                else
                    h -= 2
                end
                println(f, "\\draw<$slideno-$(slideno + 1)> [$style, ->] ($t) -- ($h);")
            end
        end
        for e in attacks[nit]
            t, h = tails[e], heads[e]
            if t == 1
                t = "s"
            elseif t == nnodes
                t = "t"
            else
                t -= 2
            end
            if h == 1
                h = "s"
            elseif h == nnodes
                h = "t"
            else
                h -= 2
            end
            println(f, "\\draw<$(slideno + 1)> [att, ->] ($t) -- ($h);")
        end
        slideno += 2
    end
    println(f, "\\end{tikzpicture}")
    println(f, "%\\caption{Original Grid network, \$10\\times 10\$}")
    println(f, "%\\end{figure}")
    println(f, "\\end{document}")
    close(f)
end

function knapsack_load(nodes, demands)
    if isempty(nodes) return 0
    else return sum(demands[u] for u in nodes)
    end
end

function build_ucinet(filename)
    nnodes = 0
    inTSP = falses(0)
    arcs = []
    open(filename) do f
        readline(f)
        let
            line = readline(f)
            tok = split(line, "=")
            nnodes = parse(Int64, tok[2])
            inTSP = falses(nnodes)
        end
        readline(f)
        readline(f)
        while !eof(f)
            line = readline(f)
            tok = split(line, " ")
            t = parse(Int64, tok[1])
            h = parse(Int64, tok[2])
            w = parse(Int64, tok[3])
            if h - t == 1 || (t == nnodes && h == 1)
        		inTSP[t] = true
            end
            push!(arcs, (t, h, w))
        end
    end
    tails = [arcs[i][1] for i in eachindex(arcs)]
    heads = [arcs[i][2] for i in eachindex(arcs)]
    weights = [arcs[i][3] for i in eachindex(arcs)]
    inTSP = falses(nnodes)
    for i in 1 : length(tails)
        t = tails[i]
        h = heads[i]
        if h - t == 1 || (t == nnodes && h == 1)
            inTSP[t] = true
        end
    end
    inf = 2 * maximum(weights)
    for t in 1 : nnodes
        if !inTSP[t]
            h = (t % nnodes) + 1
            push!(tails, t)
            push!(heads, h)
            push!(weights, inf)
        end
    end
    nedges = length(tails)
    nnodes, nedges, tails, heads, weights
end

function build_gr_road_network(distance_fn, time_fn)
    dfile = open(distance_fn)
    tfile = open(time_fn)
    nnodes = nedges = 0
    tails = Int64[]
    heads = Int64[]
    weights = Int64[]
    attacks = Int64[]
    inTSP = []
    while !eof(dfile)
        line = readline(dfile)
        tok = split(line, " ")
        if tok[1] == "c" continue
        elseif tok[1] == "p"
            nnodes = parse(Int64, tok[3])
            nedges = parse(Int64, tok[4])
            inTSP = falses(nnodes)
        elseif tok[1] == "a"
            t = parse(Int64, tok[2])
            h = parse(Int64, tok[3])
            d = parse(Int64, tok[4])
            push!(tails, t)
            push!(heads, h)
            push!(weights, d)
            # push!(tails, h)
            # push!(heads, t)
            # push!(weights, d)
            if h - t == 1 || (t == nnodes && h == 1)
        		inTSP[t] = true
        	# elseif t - h == 1 || (h == nnodes && t == 1)
        	# 	inTSP[h] = true
        	end
        end
    end
    while !eof(tfile)
        line = readline(tfile)
        tok = split(line, " ")
        if tok[1] == "c" continue
        elseif tok[1] == "p" continue
        elseif tok[1] == "a"
            t = parse(Int64, tok[4])
            push!(attacks, t)
            # push!(attacks, t)
        end
    end
    close(tfile)
    close(dfile)
    inf = 10000#sum(weights + atts)
    for i in 1 : nnodes
    	t, h = i, (i % nnodes) + 1
    	if !inTSP[t]
            push!(tails, t)
            push!(heads, h)
            push!(weights, inf)
            push!(attacks, inf)
        end
    end

    arcs = []
    for i in eachindex(tails)
        push!(arcs, (tails[i], heads[i], weights[i], attacks[i]))
    end
    sort!(arcs)
    unique!(x -> (x[1], x[2]), arcs)
    tails = [arcs[i][1] for i in eachindex(arcs)]
    heads = [arcs[i][2] for i in eachindex(arcs)]
    weights = [arcs[i][3] for i in eachindex(arcs)]
    attacks = [arcs[i][4] for i in eachindex(arcs)]
    nedges = length(tails)
    nnodes, nedges, tails, heads, weights, attacks
end
function build_road_network(filename)
    f = open(filename)
    nnodes = parse(Int64, readline(f))
    for l in 1 : nnodes
        readline(f)
    end
    arcs = []
    inTSP = falses(nnodes)
    nedges = parse(Int64, readline(f))
	for e in 1 : nedges
        line = readline(f)
        tok = split(line, " ")
        t, h = parse(Int64, tok[1]) + 1,  parse(Int64, tok[2]) + 1
        if t == h
            readline(f)
            continue
        end
        if h - t == 1 || (t == nnodes && h == 1)
    		inTSP[t] = true
    	elseif t - h == 1 || (h == nnodes && t == 1)
    		inTSP[h] = true
    	end
        line = readline(f)
        tok = split(line, " ")
        w, aw = round(Int64, parse(Float64, tok[1])),  round(Int64, parse(Float64, tok[2]))
        push!(arcs, (t, h, w, aw))
        push!(arcs, (h, t, w, aw))
    end
    inf = 10000#sum(weights + atts)
    for i in 1 : nnodes
    	t, h = i, (i % nnodes) + 1
    	if !inTSP[t]
            push!(arcs, (t, h, inf, inf))
        end
    end
    sort!(arcs)
    unique!(x -> (x[1], x[2]), arcs)
    tails = [arcs[i][1] for i in eachindex(arcs)]
    heads = [arcs[i][2] for i in eachindex(arcs)]
    weights = [arcs[i][3] for i in eachindex(arcs)]
    atts = [arcs[i][4] for i in eachindex(arcs)]
    nedges = length(tails)
    nnodes, nedges, tails, heads, weights, atts
end

function build_random_road_nets(nnodes, nedges, tails, heads, weights, atts, N, state)
    for i in 1 : N
        filename = "./instances/shortestpath/road/$state$i.txt"
        f = open(filename, "w")
        println(f, "sp min $nnodes $nedges")
        s = t = 0
        while s == t
            s = rand(1 : nnodes)
            t = rand(1 : nnodes)
        end
        println(f, "$s 1")
        println(f, "$t -1")
        for e in 1 : nedges
            println(f, "$(tails[e]) $(heads[e]) $(weights[e]) $(atts[e])")
        end
        close(f)
        println("instance written for $state$i")
    end
end

function objectsFitInKnapsack(already, demands, budget)
    nedges = length(demands)
    kl = knapsack_load(already, demands)
    [e for e in 1 : nedges if !in(e, already) && demands[e] <= budget - kl]
end
