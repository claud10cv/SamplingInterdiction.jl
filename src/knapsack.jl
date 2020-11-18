using JuMP
using Gurobi#, CPLEX
using PisingerKnapsack

function knapsack_demcap(cap)::Function
    fixed_knapsack(tails, heads, weights, demands, adjlist, attack) = knapsack(tails, heads, weights, demands, adjlist, attack, cap)
    return fixed_knapsack
end
function knapsack(tails, heads, weights, demands, adjlist, attack, cap)
    inds = [i for i in 1 : length(weights) if weights[i] < -1e-7]
    noval = [i for i in 1 : length(weights) if abs(weights[i]) <= 1e-7 && attack[i] > 1e-7]
    if isempty(inds) return 0, [], []
    end
    #println("inds = $inds")
    optval, chosen = knapsack(-weights[inds], demands[inds], cap)
    chosen = inds[chosen]
    Q = sum(demands[chosen])
    sort!(noval, lt = (x, y) -> attack[x] > attack[y])
    for v in noval
	if Q + demands[v] <= cap
#		println("adding element")
		push!(chosen, v)
		Q += demands[v]
		optval -= weights[v]
	end
    end
    sort!(chosen)
    -optval, chosen, []
end
function knapsack(vals, demands, cap)
	return knapsack_minknap(vals, demands, cap)
    ram_dp = (length(vals) + 1) * (cap + 1) * sizeof(Float64) * 1e-6
    if ram_dp >= 1024
#	println("running mip")
        return knapsack_mip(vals, demands, cap)
    else
#	println("running dp")
        return knapsack_dp(vals, demands, cap)
    end
end

function knapsack_minknap(vals, demands, cap)
	nitems = length(vals)
#	println(vals)
#	println(demands)
#	println(cap)
	sumdem = sum(demands)
	if sumdem <= cap
		return sum(vals), [i for i in 1 : nitems]
	else
		obj, sol = minknap(vals, demands, cap)
		chosen = [i for i in 1 : nitems if sol[i] == 1]
		obj, chosen 
	end
end

function knapsack_mip(vals, demands, cap)
    println("solving MIP knapsack")
    m = Model(optimizer_with_attributes(Gurobi.Optimizer,
                                    "Threads" => 1,
                                    "MIPGapAbs" => 0.0,
                                    "MIPGap" => 0.0,
                                    "OutputFlag" => 0))

    # m = Model(witsolver = CplexSolver(CPXPARAM_Threads = 1,
    #                                 CPXPARAM_MIP_Tolerances_AbsMIPGap = 0,
    #                                 CPXPARAM_MIP_Tolerances_MIPGap = 0,
    #                                 CPXPARAM_ScreenOutput = 0))

    # m = Model(solver = GurobiSolver(Threads = 1,
    #                                 MIPGapAbs = 0,
    #                                 MIPGap = 0,
    #                                 OutputFlag = 0))

    nnodes = length(vals)
    @variable(m, x[1 : nnodes], Bin)
    @constraint(m, dot(demands, x) <= cap)
    @objective(m, Max, dot(x, vals))
    optimize!(m)
    xvals = round.(Int64, value.(x))
    optval = dot(xvals, vals)
    if abs(optval - round(Int64, optval)) < 1e-5
        optval = round(Int64, optval)
    end
    chosen = [i for i in 1 : nnodes if xvals[i] > 0]
    optval, chosen
end

function knapsack_dp(vals, demands, cap)
    nnodes = length(vals)
    sup = [i for i in 1 : nnodes if vals[i] > 1e-5]
    nonsup = [i for i in 1 : nnodes if vals[i] <= 1e-5]
    nsup = length(sup)
    v = zeros(nsup + 1, cap + 1)
#    println("size: $(sizeof(v))")
    for u in 0 : nsup
        for q in 0 : cap
            if u == 0 || q == 0
                v[u + 1, q + 1] = 0
            elseif demands[sup[u]] > q
                v[u + 1, q + 1] = v[u, q + 1]
            else
                v[u + 1, q + 1] = max(v[u, q + 1], vals[sup[u]] + v[u, q + 1 - demands[sup[u]]])
            end
        end
    end
    chosen = []
    optval = v[nsup + 1, cap + 1]
    if abs(optval - round(Int64, optval)) < 1e-5
        optval = round(Int64, optval)
    end
    let
        u = nsup
        q = cap
        while u > 0 && q > 0
            if v[u + 1, q + 1] != v[u, q + 1]
                push!(chosen, u)
                q = q - demands[sup[u]]
                u -= 1
            else
                u -= 1
            end
        end
    end
    supchosen = sup[chosen]
    availcap = cap
    if !isempty(supchosen)
        availcap -= sum(demands[u] for u in supchosen)
    end
    shuffle!(nonsup)
    for u in nonsup
        if demands[u] <= availcap
            push!(supchosen, u)
            availcap -= demands[u]
        end
    end
    optval, supchosen
end
