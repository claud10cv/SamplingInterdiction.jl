using SparseArrays
using Dates

function mst_interdiction_sampling(tails, heads, weights, attacked_weights, leadercons, budget; max_time = 604800)
    nedges = length(tails)
#    demands = ones(Int64, nedges)
    followercons = zeros(Int64, nedges)
    return abstract_interdiction_sampling("mst",
                                        tails,
                                        heads,
                                        weights,
                                        leadercons,
                                        followercons,
                                        attacked_weights,
                                        budget,
                                        max_time,
                                        mst;
                                        is_symmetric = true)
end

function cflp_interdiction_sampling(xnat, caps, dems, fcosts, attfcosts, asscosts, budget; max_time = 604800)
    nfacs = length(caps)
    ncusts = length(dems)
    tails = [i for i in 1 : nfacs]
    heads = [i for i in 1 : nfacs]
    weights = [fcosts[i] for i in 1 : nfacs]
    attweights = [attfcosts[i] for i in 1 : nfacs]
    leadercons = ones(Int64, nfacs)
    followercons = [caps[i] for i in 1 : nfacs]
    for i in 1 : nfacs, j in 1 : ncusts
        push!(tails, i)
        push!(heads, j + nfacs)
        push!(weights, asscosts[i, j])
        push!(attweights, 0)
        push!(leadercons, budget + 1)
        push!(followercons, dems[j])
    end
    return abstract_interdiction_sampling("cflp",
                                        tails,
                                        heads,
                                        weights,
                                        leadercons,
                                        followercons,
                                        attweights,
                                        budget,
                                        max_time,
                                        cflp_abs(xnat);
                                        is_symmetric = true)
end

function knapsack_interdiction_sampling(vals, att_vals, leadercons, followercons, cap, budget; max_time = 604800)
    nnodes = length(vals)
    tails = ones(Int64, nnodes)
    heads = collect(1 : nnodes)
    return abstract_interdiction_sampling("knapsack",
                                        tails,
                                        heads,
                                        -vals,
                                        leadercons,
                                        followercons,
                                        att_vals,
                                        budget,
                                        max_time,
                                        knapsack_demcap(cap);
                                        is_symmetric = true)
end
function sp_interdiction_sampling(tails, heads, weights, attacked_weights, leadercons, budget, s, t; max_time = 604800)
    nedges = length(tails)
    followercons = zeros(Int64, nedges)
    return abstract_interdiction_sampling("spp",
                                        tails,
                                        heads,
                                        weights,
                                        leadercons,
                                        followercons,
                                        attacked_weights,
                                        budget,
                                        max_time,
                                        spp_st(s, t);
                                        is_symmetric = false)
end

# function matching_interdiction_sampling(tails, heads, weights, attacked_weights, leadercons, budget; max_time = 604800)
#     nedges = length(tails)
#     followercons = zeros(Int64, nedges)
#     return abstract_interdiction_sampling("matching",
#                                         tails,
#                                         heads,
#                                         weights,
#                                         leadercons,
#                                         followercons,
#                                         attacked_weights,
#                                         budget,
#                                         max_time,
#                                         matching;
#                                         is_symmetric = true)
# end

function bipartite_matching_interdiction_sampling(tails, heads, weights, attacked_weights, leadercons, budget; max_time = 604800)
    nedges = length(tails)
    followercons = zeros(Int64, nedges)
    return abstract_interdiction_sampling("bipartite matching",
                                        tails,
                                        heads,
                                        -weights,
                                        leadercons,
                                        followercons,
                                        attacked_weights,
                                        budget,
                                        max_time,
                                        bip_matching;
                                        is_symmetric = true)
end

function abstract_interdiction_sampling(ptype, tails, heads, weights, leadercons, followercons, attacked_weights, budget, max_time, def::Function; is_symmetric = true)
#    srand(rng, 20170901)
    Random.seed!(gParams.rng, 20170901)
    gParams.heuristic_restricted_interdiction = true
    nnodes = max(maximum(tails), maximum(heads))
    nedges = length(tails)
    println("number of nodes = $nnodes")
    println("number of edges = $nedges")

    edge_lower_bounds = typemin(Int64) * ones(Int64, nedges)
    bestAttack = []
    let
        # if is_symmetric
        #     I = vcat(tails, heads)
        #     J = vcat(heads, tails)
        #     W = vact(weights, weights)
        #     Z = vcat(collect(1 : nedges), collect(1 : nedges))
        # else
        #     I = tails
        #     J = heads
        #     W = weights
        #     Z = collect(1 : nedges)
        # end
        lb, path, nothing = def(tails, heads, weights, followercons, [], zeros(nedges))
    #    println(path)
        bestAttack = Attack([], path, lb, lb)
    end

    println(bestAttack)
    optval, optsol = abstract_iterative_interdiction_sampling(bestAttack,
                                                            tails,
                                                            heads,
                                                            weights,
                                                            leadercons,
                                                            followercons,
                                                            attacked_weights,
                                                            budget,
                                                            max_time,
                                                            def,
                                                            edge_lower_bounds,
                                                            is_symmetric)

end

function abstract_iterative_interdiction_sampling(initAttack,
                                                    tails,
                                                    heads,
                                                    weights,
                                                    leadercons,
                                                    followercons,
                                                    attacked_weights,
                                                    budget,
                                                    max_time,
                                                    def,
                                                    edge_lower_bounds,
                                                    is_symmetric)

    init_time = Dates.now()
    nnodes = max(maximum(tails), maximum(heads))
    nedges = length(tails)
    lb, ub = initAttack.cost, 100 * abs(initAttack.cost)

#    draw_grid(tails, heads, weights, attacked_weights)
    let
        newweights = weights + attacked_weights
        ub, nothing, nothing = def(tails, heads, newweights, followercons, [], ones(nedges))
    end
    partial_edges = []
    added_edge = falses(nedges)
    edge_index = zeros(Int64, nedges)
    trees = []
    attacks = []
    algAttacks = []
    algPaths = []
    report = []
    let
        newatt_mst = []
        knapLoad = 0
        # println(initAttack.mst)
        for e in initAttack.mst
            if !added_edge[e]
                push!(partial_edges, e)
                knapLoad += leadercons[e]
                added_edge[e] = true
                edge_index[e] = length(partial_edges)
            end
            push!(newatt_mst, edge_index[e])
        end
        sort!(newatt_mst)
        push!(algPaths, copy(partial_edges))
        push!(trees, Tree(newatt_mst, sum(weights[initAttack.mst])))
        remobjs = objectsFitInKnapsack(partial_edges, leadercons, budget)
        while !isempty(remobjs)#knapsack_load(partial_edges, demands) < budget
            newweights = copy(weights)
            newmst = []
	    newatt = zeros(nedges)
            for e in partial_edges
                newweights[e] += attacked_weights[e]
		newatt[e] += 1
            end
            nothing, mstedges, nothing = def(tails, heads, newweights, followercons, [], newatt)
            nadded = 0
            for e in mstedges
                if !added_edge[e]
                    push!(partial_edges, e)
                    knapLoad += leadercons[e]
                    added_edge[e] = true
                    edge_index[e] = length(partial_edges)
                    nadded += 1
                end
                push!(newmst, edge_index[e])
            end
            if nadded > 0
                push!(trees, Tree(newmst, sum(weights[mstedges])))
            else
                cands = shuffle(gParams.rng, 1 : nedges)
                for e in cands
                    if knapsack_load(partial_edges, leadercons) >= budget
                        break
                    elseif !added_edge[e]
                        push!(partial_edges, e)
                        added_edge[e] = true
                        edge_index[e] = length(partial_edges)
                    end
                end
                break
            end
            remobjs = objectsFitInKnapsack(partial_edges, leadercons, budget)
        end
    end
    bestAttack = initAttack
    numits = 0
    firstnew = 1
    while ub > lb
        numits += 1
        restricted_tails = [tails[e] for e in partial_edges]
        restricted_heads = [heads[e] for e in partial_edges]
        restricted_weights = [weights[e] for e in partial_edges]
        restricted_leadercons = [leadercons[e] for e in partial_edges]
        restricted_followercons = [followercons[e] for e in partial_edges]
        restricted_attacked_weights = [attacked_weights[e] for e in partial_edges]

        for e in partial_edges
            if edge_lower_bounds[e] > ub
                println("edge not needed")
            end
        end
	    ftime = Dates.now()
        elapsed_time = round(Int64, Dates.value(ftime - init_time) / 100) / 10
        tilim = max_time - elapsed_time
        newub, sols, opt, newtrees, status = net_restricted_interdict(restricted_tails,
                                                    restricted_heads,
                                                    restricted_weights,
                                                    restricted_leadercons,
                                                    restricted_followercons,
                                                    restricted_attacked_weights,
                                                    firstnew,
                                                    trees,
                                                    budget,
                                                    lb,
                                                    ub,
                                                    attacks,
                                                    def,
                                                    is_symmetric,
                                                    tilim)

        ub = min(ub, newub)
        tokeep = [e for e in 1 : nedges if edge_lower_bounds[e] <= ub]
        sols_full = []
        added_to_alg = false
        for att in sols
            new_weights = copy(weights)
            lcons = sum(restricted_leadercons[e] for e in att.attacked)
            if lcons > budget
                println("infeasible attack")
                exit(0)
            end
            lcons2 = 0
	    attack = zeros(Int64, nedges)
            for e in att.attacked
                f = partial_edges[e]
		attack[f] = 1
                lcons2 += leadercons[f]
                new_weights[f] += attacked_weights[f]
            end
            if lcons != lcons2
                println("leader consumptions differ")
                exit(0)
            end
            newmstcost, newmst, nothing = def(tails, heads, new_weights, followercons, [], attack)
            let
                full_attacked = [partial_edges[e] for e in att.attacked]
                full_mst = newmst
                if !added_to_alg
                    push!(algAttacks, copy(full_attacked))
                    push!(algPaths, copy(full_mst))
                    added_to_alg = true
                end
                full_cost = newmstcost
                full_ub = att.ub
                att_full = Attack(full_attacked, full_mst, full_cost, full_ub)
                push!(sols_full, att_full)
            end
            lb = max(lb, newmstcost)
        end

        firstnew = length(partial_edges) + 1
        nedges_added = 0
        let
            if !isempty(sols_full)
                params = [(1, 0, 1)]
                solsSet = Set()
                for (alpha, beta, gamma) in params
                    added = 0
                    sort!(sols_full, lt = (u, v) -> u.ub * alpha + u.cost * beta > v.ub * alpha + v.cost * beta)
                    largestdiff = alpha * sols_full[1].ub + beta * sols_full[1].cost
                    for s in sols_full
                        diff = alpha * s.ub + beta * s.cost
                        gap = abs(largestdiff - diff) / max(1e-10, abs(largestdiff)) * 100
                    #    println("gap = $gap")
                        if added <= 0 || gap < gamma
                            push!(solsSet, s)
                            added += 1
                        end
                    end
                end
                solsVec = collect(solsSet)
                for s in solsVec
                    this_edges_added = 0
                    for e in s.mst
                        if !added_edge[e]
                            nedges_added += 1
                            this_edges_added += 1
                            added_edge[e] = true
                            push!(partial_edges, e)
                            edge_index[e] = length(partial_edges)
                        end
                    end
                    nt = edge_index[s.mst]
                    sort!(nt)
                    if !in(0, nt)
                        newtree = Tree(deepcopy(nt), s.cost)
                        if !find_tree(newtree, newtrees)
                            push!(newtrees, newtree)
                        end
                    end
                end
            end
        end

        for t in newtrees
            if !find_tree(t, trees)
                push!(trees, deepcopy(t))
            end
        end

        for s in sols_full
            if s.cost > bestAttack.cost
                bestAttack = s
            end
            let
                attacked = edge_index[s.attacked]
                lcons = sum(leadercons[e] for e in s.attacked)
                lcons2 = sum(restricted_leadercons[e] for e in attacked)
                if lcons != lcons2
                    println("full and partial leader cons differ")
                    exit(0)
                end
                if lcons > budget
                    println("leader cons infeasible in full")
                    exit(0)
                end
                tree = edge_index[s.mst]
                att = Attack(deepcopy(attacked), deepcopy(tree), s.cost, s.ub)
                if !find_attack(att, attacks)
                    push!(attacks, att)
                end
            end
        end
        ftime = Dates.now()
        elapsed_time = round(Int64, Dates.value(ftime - init_time) / 100) / 10
        printub = ub < typemax(Int64) ? ub : "inf"
        gap = ub < typemax(Int64) ? round(Int64, (ub - lb) / abs(lb) * 1000) / 10 : "x"
        println("I: $numits\tL: $lb\tU: $printub\tG: $gap\tE: $(length(partial_edges))\tE+: $nedges_added\tT: $elapsed_time")
        push!(report, (lb, ub, length(partial_edges), elapsed_time))
    	if status == :timelimit
    		break
    	end
    end

#    draw_algorithm(tails, heads, algPaths, algAttacks)
    report, bestAttack
end

function compute_gap(lb, ub)
    round(Int64, (ub - lb) / abs(ub) * 10000) / 100
end
