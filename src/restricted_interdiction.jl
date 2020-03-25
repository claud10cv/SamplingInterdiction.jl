using JuMP
using Gurobi, CPLEX
using MathOptInterface

const MOI = MathOptInterface

function net_restricted_interdict(tails,
                                heads,
                                weights,
                                leadercons,
                                followercons,
                                attacked_weights,
                                firstnew,
                                trees,
                                budget,
                                lb,
                                ub,
                                attacks,
                                def::Function,
                                is_symmetric, tilim)

    newattacks = []
    if gParams.heuristic_restricted_interdiction
        newub, newattacks, newopt, newtrees = net_restricted_interdict_heu(tails,
                                                                        heads,
                                                                        weights,
                                                                        leadercons,
                                                                        followercons,
                                                                        attacked_weights,
                                                                        trees,
                                                                        budget,
                                                                        lb + 1,
                                                                        ub,
                                                                        attacks,
                                                                        def,
                                                                       is_symmetric)
    end
    if isempty(newattacks)
        gParams.heuristic_restricted_interdiction = false
        newub, newattacks, newopt, newtrees, status = net_restricted_interdict_binsearch(tails,
                                                                                heads,
                                                                                weights,
                                                                                leadercons,
                                                                                followercons,
                                                                                attacked_weights,
                                                                                firstnew,
                                                                                trees,
                                                                                budget,
                                                                                lb + 1,
                                                                                ub,
                                                                                attacks,
                                                                                def,
                                                                                is_symmetric, tilim)
    	if status == :timelimit
    		return newub, newattacks, newopt, newtrees, :timelimit
    	elseif isempty(newattacks)
    		return lb, newattacks, newopt, newtrees, :optimal
    	else
                    return newub, newattacks, newopt, newtrees, :feasible
    	end
    else
        return newub, newattacks, newopt, newtrees, :feasible
    end
end

function net_restricted_interdict_heu(tails,
                                    heads,
                                    weights,
                                    leadercons,
                                    followercons,
                                    attacked_weights,
                                    trees,
                                    budget,
                                    lb,
                                    ub,
                                    attacks,
                                    def::Function,
                                    is_symmetric)
    nnodes = max(maximum(tails), maximum(heads))
    nedges = length(tails)
    ntrees = length(trees)
    candsEdges = collect(1 : nedges)
    randTreesToEval = copy(trees)

    edgesPerLoad = [[] for q in 0 : budget]
    for e in 1 : nedges
        for q in leadercons[e] : budget
            push!(edgesPerLoad[q], e)
        end
    end


    edgeInTree = Dict{Tree, BitArray{1}}()
    for t in randTreesToEval
    #    println(t.edges)
        edgeInTree[t] = falses(nedges)
        for e in t.edges
            edgeInTree[t][e] = true
        end
    end

    function eval_attack(t::Tree, attacked::Array{Int64, 1})
        w = t.cost
        for e in attacked
            if edgeInTree[t][e]
                w += attacked_weights[e]
            end
        end
        return w
    end


    function eval_move(att::Attack, del::Int64, add::Int64; abort_at = typemin(Int64))
        edel = att.attacked[del]
        minw = typemax(Int64)
        for t in randTreesToEval
            w = eval_attack(t, att.attacked)
            if edgeInTree[t][edel]
                w -= attacked_weights[edel]
            end
            if edgeInTree[t][add]
                w += attacked_weights[add]
            end
            minw = min(minw, w)
            if minw <= abort_at return minw
            end
        end
        return minw
    end

    function eval_move(att::Attack, del::Int64, add::Tuple{Int64, Int64}; abort_at = typemin(Int64))
        edel = att.attacked[del]
        minw = typemax(Int64)
        for t in randTreesToEval
            w = eval_attack(t, att.attacked)
            if edgeInTree[t][edel]
                w -= attacked_weights[edel]
            end
            if edgeInTree[t][add[1]]
                w += attacked_weights[add[1]]
            end
            if edgeInTree[t][add[2]]
                w += attacked_weights[add[2]]
            end
            minw = min(minw, w)
            if minw <= abort_at return minw
            end
        end
        return minw
    end
    function eval_move(att::Attack, del::Tuple{Int64, Int64}, add::Int64; abort_at = typemin(Int64))
        edel = [att.attacked[del[1]], att.attacked[del[2]]]
        minw = typemax(Int64)
        for t in randTreesToEval
            w = eval_attack(t, att.attacked)
            if edgeInTree[t][edel[1]]
                w -= attacked_weights[edel[1]]
            end
            if edgeInTree[t][edel[2]]
                w -= attacked_weights[edel[2]]
            end
            if edgeInTree[t][add]
                w += attacked_weights[add]
            end
            minw = min(minw, w)
            if minw <= abort_at return minw
            end
        end
        return minw
    end

    function eval_attack(attacked::Array{Int64, 1})
        minw = typemax(Int64)
        for t in randTreesToEval
            w = eval_attack(t, attacked)
            minw = min(minw, w)
        end
        minw
    end


    attackSeeds = []
    let
        bestAttackedVec = []
        for r in 1 : 500
            attacked = Int64[]
            alledges = shuffle(gParams.rng, collect(1 : nedges))
            for e in alledges
                if knapsack_load(attacked, leadercons) + leadercons[e] <= budget
                    push!(attacked, e)
                end
            end
            w = eval_attack(attacked)
            push!(bestAttackedVec, (attacked, w))
        end
        sort!(bestAttackedVec, lt = (u, v) -> u[2] > v[2])
        resize!(bestAttackedVec, 50)

        for (attacked, w) in bestAttackedVec
            new_weights = copy(weights)
	    newatt = zeros(nedges)
            for e in attacked
                new_weights[e] += attacked_weights[e]
		newatt[e] += 1
            end
            newlb, mst_edges, nothing = def(tails, heads, new_weights, followercons, [], newatt)
            att = Attack(attacked, mst_edges, newlb, newlb)
            push!(attackSeeds, att)
        end
        sort!(attackSeeds, lt = (u, v) -> u.cost > v.cost)
        resize!(attackSeeds, 5)
    end

    currAttacks = []
    newTrees = []

    maxb = maximum([length(as.attacked) for as in attackSeeds])
    randEdges = [(i, j) for i in 1 : maxb + 2 for j in candsEdges]
    numshakes = min(maxb, 1)
    maxshakes = 1

    maxTrees = 100
    sort!(randTreesToEval, lt = (u, v) -> u.cost < v.cost)
    sort!(newTrees, lt = (u, v) -> u.cost < v.cost)
    for seed in attackSeeds
        # if !check_attacks([seed], leadercons, budget)
        #     println("wrong seed!")
        # end
        bestAttack = deepcopy(seed)
        if bestAttack.cost >= lb
            if !find_attack(bestAttack, currAttacks)
                push!(currAttacks, deepcopy(bestAttack))
            end
            newtree = Tree(deepcopy(bestAttack.mst), bestAttack.cost)
            if !find_tree(newtree, newTrees)
                push!(newTrees, newtree)
                edgeInTree[newtree] = falses(nedges)
                for k in newtree.edges
                    edgeInTree[newtree][k] = true
                end
            end
        end
        for shake in 1 : maxshakes
            # println("beginning shake")
            pivotAttack = performShake(bestAttack,
                                        candsEdges,
                                        tails,
                                        heads,
                                        weights,
                                        leadercons,
                                        followercons,
                                        budget,
                                        attacked_weights,
                                        def,
                                        is_symmetric,
                                        numshakes)
            if pivotAttack.cost > bestAttack.cost
                bestAttack = deepcopy(pivotAttack)
            end
            if pivotAttack.cost >= lb
                if !find_attack(pivotAttack, currAttacks)
                #    println("attack in heu (pivotAttack) = $(pivotAttack.attacked)")
                    push!(currAttacks, deepcopy(pivotAttack))
                end
                newtree = Tree(deepcopy(pivotAttack.mst), pivotAttack.cost)
                if !find_tree(newtree, newTrees)
                    push!(newTrees, newtree)
                    edgeInTree[newtree] = falses(nedges)
                    for k in newtree.edges
                        edgeInTree[newtree][k] = true
                    end
                end
            end
            improved = true
            while improved
                improved = false
                shuffle!(gParams.rng, randEdges)
                nonimprovits = 0
                for (g, f) in randEdges
                    if nonimprovits > 50 break
                    end
                    nonimprovits += 1
                    if length(currAttacks) >= 5
                         return ub, currAttacks, [], newTrees
                    end
                    if g > length(pivotAttack.attacked) continue
                    end
                    if length(randTreesToEval) > maxTrees
                        sort!(randTreesToEval, lt = (u, v) -> eval_attack(u, bestAttack.attacked) < eval_attack(v, bestAttack.attacked))
                        resize!(randTreesToEval, maxTrees)
                    end
                    if length(newTrees) > maxTrees
                        sort!(newTrees, lt = (u, v) -> eval_attack(u, bestAttack.attacked) < eval_attack(v, bestAttack.attacked))
                        resize!(newTrees, maxTrees)
                    end
                    e = pivotAttack.attacked[g]
                    if !in(f, pivotAttack.attacked)
                        attload = sum(leadercons[e] for e in pivotAttack.attacked)
                        attload += leadercons[f] - leadercons[e]
                        deltaload = budget - attload
                        others_add = deltaload <= 0 ? [] : [k for k in edgesPerLoad[deltaload] if !in(k, pivotAttack.attacked) && k != e && k != f]
                        others_rm = deltaload >= 0 ? [] : [k for k in eachindex(pivotAttack.attacked) if k != g && deltaload + leadercons[pivotAttack.attacked[k]] >= 0]
                        if deltaload < 0 && isempty(others_rm) continue
                        end
                        if isempty(others_add) && isempty(others_rm)
                            eval = eval_move(pivotAttack, g, f; abort_at = pivotAttack.cost)
                            delta_attack = pivotAttack.cost - eval
                            if delta_attack >= 0 continue
                            end
                            #TEST THE SWAPPING OF THE TWO EDGES
                            new_weights = copy(weights)
			    newatt = zeros(nedges)
                            for i in 1 : length(pivotAttack.attacked)
                                a = i == g ? f : pivotAttack.attacked[i]
                                new_weights[a] += attacked_weights[a]
				newatt[a] += 1
                            end

                            newlb, mst_edges, nothing = def(tails, heads, new_weights, followercons, [], newatt)
                            if newlb > pivotAttack.cost
                                newatt = copy(pivotAttack.attacked)
                                newatt[g] = f
                                sort!(newatt)
                                newPivot = deepcopy(Attack(newatt, mst_edges, newlb, newlb))
                                pivotAttack = deepcopy(newPivot)
                                if newlb >= lb
                                    if !find_attack(newPivot, currAttacks)
                                        push!(currAttacks, deepcopy(newPivot))
                                    end
                                    newtree = Tree(deepcopy(mst_edges), newlb)
                                    if length(newTrees) < maxTrees || (newlb < newTrees[end].cost && !find_tree(newtree, newTrees))
                                        push!(newTrees, newtree)
                                        sort!(newTrees, lt = (u, v) -> u.cost < v.cost)
                                        edgeInTree[newtree] = falses(nedges)
                                        for k in newtree.edges
                                            edgeInTree[newtree][k] = true
                                        end
                                    end
                                end
                                if newlb > bestAttack.cost
                                    bestAttack = deepcopy(newPivot)
                                end
                                improved = true
                                nonimprovits = 0
                                continue
                            else
                                nt = Tree(deepcopy(mst_edges), newlb)
                                if length(randTreesToEval) < maxTrees || (newlb < randTreesToEval[end].cost && !find_tree(nt, randTreesToEval))
                                    push!(randTreesToEval, nt)
                                    sort!(randTreesToEval, lt = (u, v) -> u.cost < v.cost)
                                    edgeInTree[nt] = falses(nedges)
                                    for k in nt.edges
                                        edgeInTree[nt][k] = true
                                    end
                                end
                                #println("number of new trees = $(length(newTrees))")
                                if length(newTrees) < maxTrees || (newlb < newTrees[end].cost && !find_tree(nt, newTrees))
                                    push!(newTrees, nt)
                                    sort!(newTrees, lt = (u, v) -> u.cost < v.cost)
                                    edgeInTree[nt] = falses(nedges)
                                    for k in nt.edges
                                        edgeInTree[nt][k] = true
                                    end
                                end
                            end
                        elseif !isempty(others_add)
                            shuffle!(gParams.rng, others_add)
                            for l in others_add
                                if in(f, pivotAttack.attacked) || in(l, pivotAttack.attacked) continue
                                end
                                eval = eval_move(pivotAttack, g, (f, l); abort_at = pivotAttack.cost)
                                delta_attack = pivotAttack.cost - eval
                                if delta_attack >= 0 continue
                                end
                                #TEST THE SWAPPING OF THE TWO EDGES
                                new_weights = copy(weights)
				newatt = zeros(nedges)
                                for i in 1 : length(pivotAttack.attacked)
                                    a = i == g ? f : pivotAttack.attacked[i]
                                    new_weights[a] += attacked_weights[a]
				    newatt[a] += 1
                                end
                                if l != f 
					new_weights[l] += attacked_weights[l]
					newatt[l] += 1
                                end

                                newlb, mst_edges, nothing = def(tails, heads, new_weights, followercons, [], newatt)
                                if newlb > pivotAttack.cost
                                    newatt = copy(pivotAttack.attacked)
                                    newatt[g] = f
                                    if l != f push!(newatt, l)
                                    end
                                    sort!(newatt)
                                    newPivot = deepcopy(Attack(newatt, mst_edges, newlb, newlb))
                                    pivotAttack = deepcopy(newPivot)
                                    if newlb >= lb
                                        if !find_attack(newPivot, currAttacks)
                                            push!(currAttacks, deepcopy(newPivot))
                                        end
                                        newtree = Tree(deepcopy(mst_edges), newlb)
                                        if length(newTrees) < maxTrees || (newlb < newTrees[end].cost && !find_tree(newtree, newTrees))
                                            push!(newTrees, newtree)
                                            sort!(newTrees, lt = (u, v) -> u.cost < v.cost)
                                            edgeInTree[newtree] = falses(nedges)
                                            for k in newtree.edges
                                                edgeInTree[newtree][k] = true
                                            end
                                        end
                                    end
                                    if newlb > bestAttack.cost
                                        bestAttack = deepcopy(newPivot)
                                    end
                                    improved = true
                                    nonimprovits = 0
                                    break
                                else
                                    nt = Tree(deepcopy(mst_edges), newlb)
                                    if length(randTreesToEval) < maxTrees || (newlb < randTreesToEval[end].cost && !find_tree(nt, randTreesToEval))
                                        push!(randTreesToEval, nt)
                                        sort!(randTreesToEval, lt = (u, v) -> u.cost < v.cost)
                                        edgeInTree[nt] = falses(nedges)
                                        for k in nt.edges
                                            edgeInTree[nt][k] = true
                                        end
                                    end
                                    #println("number of new trees = $(length(newTrees))")
                                    if length(newTrees) < maxTrees || (newlb < newTrees[end].cost && !find_tree(nt, newTrees))
                                        push!(newTrees, nt)
                                        sort!(newTrees, lt = (u, v) -> u.cost < v.cost)
                                        edgeInTree[nt] = falses(nedges)
                                        for k in nt.edges
                                            edgeInTree[nt][k] = true
                                        end
                                    end
                                end
                            end
                            if improved continue
                            end
                        else
                            shuffle!(gParams.rng, others_rm)
                            for k in others_rm
                                eval = eval_move(pivotAttack, (g, k), f; abort_at = pivotAttack.cost)
                                delta_attack = pivotAttack.cost - eval
                                if delta_attack >= 0 continue
                                end
                                #TEST THE SWAPPING OF THE TWO EDGES
                                new_weights = copy(weights)
				newatt = zeros(nedges)
                                for i in 1 : length(pivotAttack.attacked)
                                    if i != k
                                        a = i == g ? f : pivotAttack.attacked[i]
                                        new_weights[a] += attacked_weights[a]
					newatt[a] += 1
                                    end
                                end

                                newlb, mst_edges, nothing = def(tails, heads, new_weights, followercons, [], newatt)
                                if newlb > pivotAttack.cost
                                    newatt = [pivotAttack.attacked[p] for p in eachindex(pivotAttack.attacked) if p != g && p != k]
                                    push!(newatt, f)
                                    sort!(newatt)
                                    newPivot = deepcopy(Attack(newatt, mst_edges, newlb, newlb))
                                    pivotAttack = deepcopy(newPivot)
                                    if newlb >= lb
                                        if !find_attack(newPivot, currAttacks)
                                            push!(currAttacks, deepcopy(newPivot))
                                        end
                                        newtree = Tree(deepcopy(mst_edges), newlb)
                                        if length(newTrees) < maxTrees || (newlb < newTrees[end].cost && !find_tree(newtree, newTrees))
                                            push!(newTrees, newtree)
                                            sort!(newTrees, lt = (u, v) -> u.cost < v.cost)
                                            edgeInTree[newtree] = falses(nedges)
                                            for h in newtree.edges
                                                edgeInTree[newtree][h] = true
                                            end
                                        end
                                    end
                                    if newlb > bestAttack.cost
                                        bestAttack = deepcopy(newPivot)
                                    end
                                    improved = true
                                    nonimprovits = 0
                                    break
                                else
                                    nt = Tree(deepcopy(mst_edges), newlb)
                                    if length(randTreesToEval) < maxTrees || (newlb < randTreesToEval[end].cost && !find_tree(nt, randTreesToEval))
                                        push!(randTreesToEval, nt)
                                        sort!(randTreesToEval, lt = (u, v) -> u.cost < v.cost)
                                        edgeInTree[nt] = falses(nedges)
                                        for h in nt.edges
                                            edgeInTree[nt][h] = true
                                        end
                                    end
                                    #println("number of new trees = $(length(newTrees))")
                                    if length(newTrees) < maxTrees || (newlb < newTrees[end].cost && !find_tree(nt, newTrees))
                                        push!(newTrees, nt)
                                        sort!(newTrees, lt = (u, v) -> u.cost < v.cost)
                                        edgeInTree[nt] = falses(nedges)
                                        for h in nt.edges
                                            edgeInTree[nt][h] = true
                                        end
                                    end
                                end
                            end
                            if improved continue
                            end
                        end
                    end
                end
            end
        end
    end
    sort!(currAttacks, lt = (u, v) -> u.cost > v.cost)
    for att in currAttacks
        lcons = sum(leadercons[e] for e in att.attacked)
        if lcons > budget
            println("infeasible attack in heuristic")
            exit(0)
        end
    end
    ub, currAttacks, nothing, newTrees
end

function check_attacks(currAttacks, leadercons, budget)
    for att in currAttacks
        edges = att.attacked
        cons = sum(leadercons[e] for e in edges)
        if cons > budget
            println("attack exceeds budget by $(cons - budget) units")
            return false
        end
    end
    return true
end

function performShake(attack, cands, tails, heads, weights, leadercons, followercons, budget, attacked_weights, def, is_symmetric, numshakes)
    nedges = length(tails)
    attacked = copy(attack.attacked)
    nonattacked = [e for e in 1 : nedges if !in(e, attacked)]
    numshakes = minimum([numshakes, length(attacked), length(nonattacked)])
    toswap = collect(1 : length(attacked))
    nshake = 0
    cons = sum(leadercons[e] for e in attacked)
    while nshake < numshakes
        feas = false
        i = e = 0
        nit = 0
        while !feas && nit < 100 * numshakes
            nit += 1
            i = rand(gParams.rng, toswap)
            e = rand(gParams.rng, nonattacked)
            newcons = cons - leadercons[attacked[i]] + leadercons[e]
            feas = newcons <= budget
        end
        if !feas break
        end
        nshake += 1
        cons += leadercons[e] - leadercons[attacked[i]]
        push!(nonattacked, attacked[i])
        filter!(x -> x != e, nonattacked)
        attacked[i] = e
    end
    new_weights = copy(weights)
    newatt = zeros(nedges)
    for a in attacked
        new_weights[a] += attacked_weights[a]
	newatt[a] += 1
    end

    newlb, mst_edges, nothing = def(tails, heads, new_weights, followercons, [], newatt)
    newatt = Attack(attacked, mst_edges, newlb, newlb)
    return newatt
end

function net_restricted_interdict_powbinsearch(tails,
                                    heads,
                                    weights,
                                    leadercons,
                                    followercons,
                                    attacked_weights,
                                    firstnew,
                                    trees,
                                    budget,
                                    lb,
                                    ub,
                                    attacks,
                                    def::Function,
                                    is_symmetric, tilim)
    newAttacks = []
    newTrees = []
    init_time = Dates.now()
    newlb, newub, newAttacksPow, opt, newTreesPow, stat = net_restricted_interdict_powsearch(tails,
                                        heads,
                                        weights,
                                        leadercons,
                                        followercons,
                                        attacked_weights,
                                        firstnew,
                                        trees,
                                        budget,
                                        lb,
                                        ub,
                                        attacks,
                                        def::Function,
                                        is_symmetric, tilim)
    if stat == :timelimit
        return newub, newAttacks, opt, newTrees, stat
    end
    for t in newTreesPow
        if !find_tree(t, newTrees)
            push!(newTrees, t)
        end
    end
    for att in newAttacksPow
        if !find_attack(att, newAttacks)
            push!(newAttacks, att)
        end
    end
    gap = (newub - newlb) / max(abs(newub), abs(newlb)) * 100
    if !isempty(newAttacks) && gap < 2
        return newub, newAttacks, opt, newTrees, stat
    end
    end_time_ps = Dates.now()
    elapsed_time_ps = Dates.value(end_time_ps - init_time) / 1000
    tilim -= elapsed_time_ps
    ub, newAttacksBS, opt, newTreesBS, stat = net_restricted_interdict_binsearch(tails,
                                        heads,
                                        weights,
                                        leadercons,
                                        followercons,
                                        attacked_weights,
                                        firstnew,
                                        trees,
                                        budget,
                                        newlb,
                                        newub,
                                        attacks,
                                        def::Function,
                                        is_symmetric, tilim)
    if stat == :timelimit
        return ub, newAttacksBS, opt, newTreesBS, stat
    end
    for t in newTreesBS
        if !find_tree(t, newTrees)
            push!(newTrees, t)
        end
    end
    for att in newAttacksBS
        if !find_attack(att, newAttacks)
            push!(newAttacks, att)
        end
    end
    return ub, newAttacks, opt, newTrees, stat
end

function net_restricted_interdict_powsearch(tails,
                                    heads,
                                    weights,
                                    leadercons,
                                    followercons,
                                    attacked_weights,
                                    firstnew,
                                    trees,
                                    budget,
                                    lb,
                                    ub,
                                    attacks,
                                    def::Function,
                                    is_symmetric, tilim)
    newAttacks = []
    newTrees = []
    newOpt = []
    currAttacks = copy(attacks)
    currTrees = copy(trees)
	initTime = Dates.now()
    base = 3
    step = 1
    initub = ub
    newlb = newub = ub
    while newlb > lb
        ftime = Dates.now()
        elapsed_time = round(Int64, Dates.value(ftime - initTime) / 100) / 10
        newtilim = tilim - elapsed_time
        newlb = max(lb, initub - step)
    #    println("ps in range [$newlb, $newub]")
        step *= base
        newub, newattacks, newopt, newtrees, stat = net_restricted_interdict_mip(tails,
                                                                            heads,
                                                                            weights,
                                                                            leadercons,
                                                                            followercons,
                                                                            attacked_weights,
                                                                            firstnew,
                                                                            currTrees,
                                                                            budget,
                                                                            newlb,
                                                                            newub,
                                                                            currAttacks,
                                                                            def::Function,
                                                                            is_symmetric, newtilim)
        ftime = Dates.now()
        elapsed_time = round(Int64, Dates.value(ftime - initTime) / 100) / 10
        for att in newattacks
            if !find_attack(att, newAttacks)
                push!(newAttacks, att)
            end
            if !find_attack(att, currAttacks)
                push!(currAttacks, att)
            end
        end
        for t in newtrees
            if !find_tree(t, newTrees)
                push!(newTrees, t)
            end
            if !find_tree(t, currTrees)
                push!(currTrees, t)
            end
        end
    	if stat == :timelimit
    		return newlb, newub, newAttacks, newopt, newTrees, stat
    	elseif isempty(newattacks)
            newub = newlb - 1
        else
            return newlb, newub, newAttacks, newopt, newTrees, :feasible
        end
    end
    return lb, lb, newAttacks, newOpt, newTrees, :feasible
end


function net_restricted_interdict_binsearch(tails,
                                    heads,
                                    weights,
                                    leadercons,
                                    followercons,
                                    attacked_weights,
                                    firstnew,
                                    trees,
                                    budget,
                                    lb,
                                    ub,
                                    attacks,
                                    def::Function,
                                    is_symmetric, tilim)

    newAttacks = []
    newTrees = []
    newOpt = []
    currAttacks = copy(attacks)
    currTrees = copy(trees)
	initTime = Dates.now()
    for att in attacks
        kl = sum(leadercons[e] for e in att.attacked)
        if kl > budget
            println("infeasible in binsearch before mip")
            exit(0)
        end
    end

    while lb <= ub
        gap = (ub - lb) / max(abs(ub), abs(lb)) * 100
        if !isempty(newAttacks) && gap < 2
            return ub, newAttacks, newOpt, newTrees, :feasible
        end
        ftime = Dates.now()
        elapsed_time = round(Int64, Dates.value(ftime - initTime) / 100) / 10
        newtilim = tilim - elapsed_time
        mid = ceil(Int64, (3 * ub + 2 * lb) / 5)
    #    println("bs in range [$lb, $ub]")
        newub, newattacks, newOpt, newtrees, stat = net_restricted_interdict_mip(tails,
                                                                            heads,
                                                                            weights,
                                                                            leadercons,
                                                                            followercons,
                                                                            attacked_weights,
                                                                            firstnew,
                                                                            currTrees,
                                                                            budget,
                                                                            mid,
                                                                            ub,
                                                                            currAttacks,
                                                                            def::Function,
                                                                            is_symmetric, newtilim)
        ftime = Dates.now()
        elapsed_time = round(Int64, Dates.value(ftime - initTime) / 100) / 10
        for att in newattacks
            kl = sum(leadercons[e] for e in att.attacked)
            if kl > budget
                println("infeasible at return of mip")
                exit(0)
            end
            if !find_attack(att, newAttacks)
                push!(newAttacks, att)
            end
            if !find_attack(att, currAttacks)
                push!(currAttacks, att)
            end
        end
        for att in currAttacks
            kl = isempty(att.attacked) ? 0 : sum(leadercons[e] for e in att.attacked)
            if kl > budget
                println("infeasible in binsearch before mip, $kl > $budget")
                exit(0)
            end
        end
        for t in newtrees
            if !find_tree(t, newTrees)
                push!(newTrees, t)
            end
            if !find_tree(t, currTrees)
                push!(currTrees, t)
            end
        end
    	if stat == :timelimit
    		return newub, newAttacks, newOpt, newTrees, stat
    	elseif isempty(newattacks)
        	ub = mid - 1
        else
            lb = mid
	    end
    end
	stat = isempty(newAttacks) ? :optimal : :feasible
    return ub, newAttacks, newOpt, newTrees, :stat
end

function net_restricted_interdict_mip(tails,
                                    heads,
                                    weights,
                                    leadercons,
                                    followercons,
                                    attacked_weights,
                                    firstnew,
                                    trees,
                                    budget,
                                    lb,
                                    ub,
                                    attacks,
                                    def::Function,
                                    is_symmetric, tilim)
	init_time = Dates.now()
    nnodes = max(maximum(tails), maximum(heads))
    nedges = length(tails)
    # m = Model(optimizer_with_attributes(CPLEX.Optimizer,
    #                                 "CPXPARAM_Threads" => 1,
    #                                 "CPXPARAM_MIP_Tolerances_AbsMIPGap" => 0.0,
    #                                 "CPXPARAM_MIP_Tolerances_MIPGap" => 0.0,
    #                                 "CPXPARAM_MIP_Tolerances_LowerCutoff" => lb,
    #                                 "CPXPARAM_Emphasis_MIP" => 1,
    #                                 "CPXPARAM_MIP_Strategy_VariableSelect" => 2,
    #                                 "CPXPARAM_ScreenOutput" => 0))

    m = Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer,
                                    "Threads" => 1,
                                    "MIPGapAbs" => 0.0,
                                    "MIPGap" => 0.0,
                                    "Cutoff" => lb - 1e-7,
                                    # "CPXPARAM_MIP_Tolerances_LowerCutoff" => lb,
                                    # "CPXPARAM_Emphasis_MIP" => 1,
                                    # "CPXPARAM_MIP_Strategy_VariableSelect" => 2,
                                    "OutputFlag" => 0))
    # m = Model(solver = CplexSolver(CPXPARAM_Threads = 1,
    #                                 CPXPARAM_MIP_Tolerances_AbsMIPGap = 0.0,
    #                                 CPXPARAM_MIP_Tolerances_MIPGap = 0.0,
    #                                 CPXPARAM_MIP_Tolerances_LowerCutoff = lb,
    #                                 CPXPARAM_Emphasis_MIP = 1,
    #                                 CPXPARAM_MIP_Strategy_VariableSelect = 2,
    #                                 CPXPARAM_ScreenOutput = 0))

    # m = Model(solver = GurobiSolver(Threads = 1,
    #                                 MIPGapAbs = 0,
    #                                 MIPGap = 0,
    #                             #    Presolve = 0,
    #                                 OutputFlag = 0))#,
    #                             #    MIPFocus = 3))

    # if ub < typemax(Int64)
    #     @variable(m, lambda <= ub + 1)
    # else
    #     @variable(m, lambda <= 1e+7)
    # end
    @variable(m, lambda <= ub + 0.1 * abs(ub))
    @variable(m, x[1 : nedges], Bin)
    @objective(m, Max, lambda)
    @constraint(m, JuMP.dot(leadercons, x) <= budget)

    rlb = lb
    attackCuts = AttackCut[]
    let
        bestatt = []
        for att in attacks
            remobjs = objectsFitInKnapsack(att.attacked, leadercons, budget)
            kl = isempty(att.attacked) ? 0 : sum(leadercons[e] for e in att.attacked)
            if kl > budget
                println("infeasible attack at the beginning: $kl > $budget")
                exit(0)
            end
            if isempty(remobjs)#length(att.attacked) >= budget
                new_weights = copy(weights)
		newatt = zeros(nedges)
                for e in att.attacked
                    new_weights[e] += attacked_weights[e]
		    newatt[e] += 1
                end

                newlb, mst_edges, nothing = def(tails, heads, new_weights, followercons, [], newatt)
                if newlb >= ub
                    sols = [att]
                    println("solution will not change")
                    return ub, sols, att, [], :feasible
                elseif newlb < lb
                    xcoefs = [in(e, att.attacked) ? 1 : 0 for e in 1 : nedges]
                    push!(attackCuts, AttackCut(0, xcoefs, length(xcoefs) - 1))
                elseif newlb > lb
                    rlb = newlb
                    bestatt = att
                end
            end
        end
        if bestatt != []
            for e in bestatt.attacked
                setvalue(x[e], 1)
            end
        end
    end

    let
        newatts = [att for att in attacks if isempty(objectsFitInKnapsack(att.attacked, leadercons, budget))]
        svi = separate_super_valid_inequalities(newatts,
                                                tails,
                                                heads,
                                                weights,
                                                attacked_weights,
                                                leadercons,
                                                followercons,
                                                rlb,
                                                def,
                                                is_symmetric)
        for (mstedges, b) in svi
            xcoefs = zeros(nedges)
            for e in mstedges
                xcoefs[e] = 1
            end
            push!(attackCuts, AttackCut(0, -xcoefs, -b))
        end
    end

    newTrees = []
    sols = []
    optsol = Attack([], [], typemin(Int64), typemax(Int64))
    function lazycb(cb)
        # xvals = getvalue(x)
        # lambdaval = getvalue(lambda)
        xvals = [callback_value(cb, x[i]) for i in eachindex(x)]
        lambdaval = callback_value(cb, lambda)
        if abs(lambdaval - round(lambdaval)) < 1e-7
            lambdaval = round(lambdaval)
        end

        for e in 1 : nedges
            if abs(xvals[e] - round(xvals[e])) < 1e-7
                xvals[e] = round(xvals[e])
            end
        end
        bestbound = ub
        # bestbound = min(ub, floor(MathProgBase.cbgetbestbound(cb) + 1e-7))
        cbatt = nothing
        let
            attacked = [e for e in 1 : nedges if xvals[e] > 1e-7]
            new_weights = copy(weights)
	    newatt = zeros(nedges)
            for e in attacked
                new_weights[e] += attacked_weights[e]
		newatt[e] += 1
            end
            newlb, mst_edges, nothing = def(tails, heads, new_weights, followercons, [], newatt)
            kl = isempty(attacked) ? 0 : sum(leadercons[e] for e in attacked)
            if kl <= budget
                push!(attacks, Attack(attacked, mst_edges, newlb, lambdaval))
                cbatt = attacks[end]
            end
        end
        let
            bencuts, newtrees = separate_benders_cuts(tails,
                                                    heads,
                                                    weights,
                                                    leadercons,
                                                    followercons,
                                                    attacked_weights,
                                                    def,
                                                    is_symmetric,
                                                    xvals,
                                                    lambdaval,
                                                    bestbound)
            # println("number of cuts found = $(length(bencuts))")
            for c in bencuts
                lcoef = c.lambda_coef
                xcoefs = c.edge_coefs
                rhs = c.rhs
                # println("lcoef = $lcoef")
                # println("xcoefs = $xcoefs")
                # println("rhs = $rhs")
                # println("viol = $(lcoef * lambdaval + JuMP.dot(xvals, xcoefs) - rhs)")
                # @lazyconstraint(cb, lcoef * lambda + JuMP.dot(x, xcoefs) <= rhs)
                newlazy = @build_constraint(lcoef * lambda + JuMP.dot(x, xcoefs) <= rhs)
                # print("adding cut...")
                MOI.submit(m, MOI.LazyConstraint(cb), newlazy)
                # println("done!")
            end
            for t in newtrees
                if !find_tree(t, newTrees)
                    push!(newTrees, t)
                end
            end
            if !isempty(bencuts) return
            end
        end
        if sum(abs.(xvals - round.(xvals))) < 1e-5
            # lcons = sum(leadercons[e] for e in cbatt.attacked)
            # if lcons > budget
            #     println("infeasible attack in lazy")
            #     exit(0)
            # end
            if !isnothing(cbatt)
                push!(sols, cbatt)
                if optsol.cost < cbatt.cost
                    optsol = cbatt
                end
            end
        end
    end

    function cutcb(cb)
        # if MathProgBase.cbgetexplorednodes(cb) > 1000 return
        # end
        # xvals = getvalue(x)
        # lambdaval = getvalue(lambda)
        xvals = [callback_value(cb, s) for s in x]
        lambdaval = callback_value(cb, lambda)
        bestbound = ub
        # bestbound = min(ub, floor(Int64, MathProgBase.cbgetbestbound(cb) + 1e-7))
        if abs(lambdaval - round(lambdaval)) < 1e-7
            lambdaval = round(lambdaval)
        end
        for e in 1 : nedges
            if abs(xvals[e] - round(xvals[e])) < 1e-7
                xvals[e] = round(xvals[e])
            end
        end
        foundcuts = false
        let
            viols = []
            for (i, acut) in enumerate(attackCuts)
                lcoef = acut.lambda_coef
                xcoef = acut.edge_coefs
                rhs = acut.rhs
                viol = (lcoef * lambdaval + JuMP.dot(xcoef, xvals) - rhs) / abs(rhs)
                if viol > 1e-1
                    push!(viols, (i, viol))
                end
            end
            sort!(viols, lt = (u, v) -> u[2] > v[2])
            if length(viols) > 10
                resize!(viols, 10)
            end
            for (i, v) in viols
                acut = attackCuts[i]
                newcut = @build_constraint(acut.lambda_coef * lambda + JuMP.dot(x, acut.edge_coefs) <= acut.rhs)
                # @usercut(cb, acut.lambda_coef * lambda + JuMP.dot(x, acut.edge_coefs) <= acut.rhs)
                MOI.submit(m, MOI.UserCut(cb), newcut)
                foundcuts = true
            end
        end
        if foundcuts return
        end
        let
            attedges = [e for e in 1 : nedges if xvals[e] > 1e-7]
            let
                sort!(lt = (u, v) -> xvals[u] < xvals[v], attedges)
                sum = 0
                while sum + xvals[attedges[1]] < 1 - 1e-1
                    sum += xvals[popfirst!(attedges)]
                end
            end

            attvec = [Attack(attedges, [], 0, 0)]
            svi = separate_super_valid_inequalities(attvec,
                                                    tails,
                                                    heads,
                                                    weights,
                                                    attacked_weights,
                                                    leadercons,
                                                    followercons,
                                                    rlb,
                                                    def,
                                                    is_symmetric)
            for (mstedges, b) in svi
                viol = b - sum(xvals[e] for e in mstedges)
                if viol > 1e-1
                    newcut = @build_constraint(sum(x[e] for e in mstedges) >= b)
                    # @usercut(cb, sum(x[e] for e in mstedges) >= b)
                    MOI.submit(m, MOI.UserCut(cb), newcut)
                    foundcuts = true
                end
            end
        end
        if foundcuts return
        end
        let
            bencuts, newtrees = separate_benders_cuts(tails,
                                                    heads,
                                                    weights,
                                                    leadercons,
                                                    followercons,
                                                    attacked_weights,
                                                    def,
                                                    is_symmetric,
                                                    xvals,
                                                    lambdaval,
                                                    bestbound)
            for c in bencuts
                lcoef = c.lambda_coef
                xcoefs = c.edge_coefs
                rhs = c.rhs
                newcut = @build_constraint(lcoef * lambda + JuMP.dot(x, xcoefs) <= rhs)
                # @usercut(cb, lcoef * lambda + JuMP.dot(x, xcoefs) <= rhs)
                MOI.submit(m, MOI.UserCut(cb), newcut)
                foundcuts = true
            end
        end
        if foundcuts return
        end
    end

	stat = :feasible
    function infocb(cb)
        xvals = round.(Int64, getvalue(x))
        lambdaval = round(Int64, getvalue(lambda))
        bestbound = MathProgBase.cbgetbestbound(cb)
        bestobj = MathProgBase.cbgetobj(cb)
        attack = [e for e in 1 : nedges if xvals[e] > 1e-7]
        new_weights = copy(weights)
        for e in 1 : nedges
            new_weights[e] += xvals[e] * attacked_weights[e]
        end

        newlb, mst_edges, nothing = def(tails, heads, new_weights, followercons, [], xvals)
        if newlb >= lb
            mst_unattacked = sum(weights[e] for e in mst_edges)
            newatt = Attack(attack, mst_edges, newlb, lambdaval)
            if optsol.cost < newlb
                optsol = newatt
            end
            if !find_attack(newatt, sols)
                push!(sols, newatt)
            end
            newtree = Tree(deepcopy(mst_edges), mst_unattacked)
            if !find_tree(newtree, newTrees)
                push!(newTrees, newtree)
            end
            if newlb > rlb
                rlb = newlb
            end
        end
    end

    function boundcheckercb(cb)
    	current_time = Dates.now()
    	elapsed_time = Dates.value(current_time - init_time) / 1000
    	if elapsed_time > tilim
    		stat = :timelimit
    		return JuMP.StopTheSolver
    	end
        if !isempty(sols)
            return JuMP.StopTheSolver
            bestbound = min(ub, ceil(Int64, MathProgBase.cbgetbestbound(cb) - 1e-7))
            relgap = abs(rlb - bestbound) / abs(rlb) * 100
            absgap = abs(rlb - bestbound)
            if length(sols) >= 2 || relgap < 5
                return JuMP.StopTheSolver
            end
        end
    end

    MOI.set(m, MOI.LazyConstraintCallback(), lazycb)
    MOI.set(m, MOI.UserCutCallback(), cutcb)
    # addcutcallback(m, cutcb)
    # addlazycallback(m, lazycb; fractional = false)
    # addinfocallback(m, infocb; when = :MIPSol)
    # addinfocallback(m, boundcheckercb; when = :MIPNode)
    optimize!(m)
    # status = solve(m; suppress_warnings = true)
	# objbound = ub
    objbound = ub
    term_status = termination_status(m)
	if in(term_status, [MOI.INFEASIBLE, MOI.OBJECTIVE_LIMIT])
        # println("has values = $(has_values(m))")
		objbound = lb - 1
		stat = :optimal
	else
        # println(term_status)
        optx = round.(Int64, value.(x))
        optlambda = round(Int64, value(lambda))
        attack = [e for e in 1 : nedges if optx[e] > 1e-7]
        new_weights = copy(weights)#[weights[e] + xvals[e] * attacked_weights[e] for e in 1 : nedges]
        for e in 1 : nedges
            new_weights[e] += optx[e] * attacked_weights[e]
        end
        newlb, mst_edges, nothing = def(tails, heads, new_weights, followercons, [], optx)
        mst_unattacked = sum(weights[e] for e in mst_edges)
        optatt = Attack(attack, mst_edges, newlb, optlambda)
        kl = sum(leadercons[e] for e in optatt.attacked)
        if kl > budget
            println("error in optimal attack in MIP")
            exit(0)
        end
        push!(sols, optatt)
        if optatt.cost > optsol.cost
            optsol = optatt
        end
		objbound = min(ub, objective_bound(m))
        if abs(objbound - round(Int64, objbound)) < 1e-7
            objbound = round(Int64, objbound)
        end
		if !isempty(sols)
			stat = :feasible
		end
	end
    # for att in sols
    #     lcons = sum(leadercons[e] for e in att.attacked)
    #     if lcons > budget
    #         println("infeasible attack in lazy")
    #         exit(0)
    #     end
    # end
    # println("lb = $lb and ub = $ub")
    # println("bbound = $objbound")
    # println("sols = $sols")
    # println("optsol = $optsol")
    # println("newtrees = $newTrees")
    # println("stat = $stat")
	objbound, sols, optsol, newTrees, stat
end
function separate_benders_cuts(tails,
                                heads,
                                weights,
                                leadercons,
                                followercons,
                                attacked_weights,
                                def,
                                is_symmetric,
                                xvals,
                                lambda,
                                ub)
    nnodes = max(maximum(tails), maximum(heads))
    nedges = length(tails)
    cuts = []
    trees = []
    new_weights = 1.0 * copy(weights)
    for e in 1 : nedges
        new_weights[e] += xvals[e] * attacked_weights[e]
    end
    newlb, mst_edges, nothing = def(tails, heads, new_weights, followercons, [], xvals)
    viol = lambda - newlb
    if viol > 1e-1
        mst_unattacked = sum(weights[e] for e in mst_edges)
        cutcoefs = zeros(nedges)
        for e in mst_edges
            ce = attacked_weights[e]
            if ce > ub - mst_unattacked
                ce = max(0, ub - mst_unattacked)
            end
            cutcoefs[e] = -ce
        end
        push!(cuts, AttackCut(1, cutcoefs, mst_unattacked))
        newtree = Tree(deepcopy(mst_edges), mst_unattacked)
        if !find_tree(newtree, trees)
            push!(trees, newtree)
        end
    end
    cuts, trees
end
function separate_super_valid_inequalities(attacks,
                                            tails,
                                            heads,
                                            weights,
                                            attacked_weights,
                                            leadercons,
                                            followercons,
                                            lb,
                                            def,
                                            is_symmetric)
    cuts = []
    nedges = length(tails)
    for att in attacks
        attacked = att.attacked
        new_weights = copy(weights)
	newatt = zeros(nedges)
        for e in attacked
            new_weights[e] += attacked_weights[e]
	    newatt[e] += 1
        end
        mstlb, mst_edges, nothing = def(tails, heads, new_weights, followercons, [], newatt)
        nonatt_lb = sum(weights[e] for e in mst_edges)
        sort!(lt = (u, v) -> attacked_weights[u] > attacked_weights[v], mst_edges)
        b = 0
        while nonatt_lb < lb
            b += 1
            nonatt_lb += attacked_weights[mst_edges[b]]
        end
        if b > 0
            push!(cuts, (mst_edges, b))
        end
    end
    cuts
end
