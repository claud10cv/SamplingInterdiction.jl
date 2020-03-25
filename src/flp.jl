using JuMP, MathOptInterface
using Gurobi#, CPLEX

const MOI = MathOptInterface

function cflp(nfacs, caps, fcost, ncusts, dems, asscost, probtype)
    m = Model(optimizer_with_attributes(Gurobi.Optimizer,
                            "Threads" => 1,
                            "MIPGapAbs" => 0,
                            "MIPGap" => 0,
                            "OutputFlag" => 0))

    # m = Model(solver = CplexSolver(CPXPARAM_Threads = 1,
    #                                 CPXPARAM_MIP_Tolerances_AbsMIPGap = 0,
    #                                 CPXPARAM_MIP_Tolerances_MIPGap = 0,
    #                                 CPXPARAM_ScreenOutput = 0))
    # MOI.set(m, MOI.RawParameter("CPXPARAM_Threads"), 1)
    # MOI.set(m, MOI.RawParameter("CPXPARAM_MIP_Tolerances_AbsMIPGap"), 0)
    # MOI.set(m, MOI.RawParameter("CPXPARAM_MIP_Tolerances_MIPGap"), 0)
    # MOI.set(m, MOI.RawParameter("CPXPARAM_ScreenOutput"), 0)

    @variable(m, y[1 : nfacs], Bin)
    if probtype == :sscflp
        @variable(m, x[1 : nfacs, 1 : ncusts], Bin)
    else
        @variable(m, 0 <= x[1 : nfacs, 1 : ncusts] <= 1)
    end
    @objective(m, Min, dot(y, fcost) + sum(asscost[i, j] * x[i, j] for i in 1 : nfacs, j in 1 : ncusts if asscost[i, j] > -1e-5))
    @constraint(m, deg[j = 1 : ncusts], sum(x[i, j] for i in 1 : nfacs if asscost[i, j] > -1e-5) >= 1)
    if probtype != :uflp
        @constraint(m, fcap[i = 1 : nfacs], sum(dems[j] * x[i, j] for j in 1 : ncusts if asscost[i, j] > -1e-5) - caps[i] * y[i] <= 0)
    end

    function strass(cb)
        xvals = [callback_value(cb, xx) for xx in x]
        yvals = [callback_value(cb, yy) for yy in y]
        mostviol = []
        for i in 1 : nfacs, j in 1 : ncusts
            viol = xvals[i, j] - yvals[i]
            if viol > 1e-1
                push!(mostviol, (i, j, viol))
            end
        end
        sort!(mostviol, lt = (u, v) -> u[3] > v[3])
        if length(mostviol) > 10
            resize!(mostviol, 10)
        end
        for (i, j, v) in mostviol
            con = @build_constraint(x[i, j] - y[i] <= 0)
            MOI.submit(m, MOI.LazyConstraint(cb), con)
#            @lazyconstraint(cb, x[i, j] - y[i] <= 0)
        end
    end
    MOI.set(m, MOI.LazyConstraintCallback(), strass)
#    addlazycallback(m, strass; fractional = true)
    optimize!(m)
    stats = termination_status(m)
    yvals = round.(Int64, value.(y))
    xvals = value.(x)
    for e in eachindex(xvals)
        if abs(xvals[e] - round(Int64, xvals[e])) < 1e-5
            xvals[e] = round(Int64, xvals[e])
        end
    end
    objval = dot(yvals, fcost) + JuMP.dot(xvals, asscost)
    if abs(objval - round(Int64, objval)) < 1e-5
        objval = round(Int64, objval)
    end
    usedfacs = [i for i in 1 : nfacs if yvals[i] > 1e-5]
    usedass = [(i, j) for i in 1 : nfacs, j in 1 : ncusts if xvals[i, j] > 1e-5 && asscost[i, j] > -1e-5]
#    println("using $(length(usedfacs)) depots")
    objval, usedfacs, usedass
end

function cflp_(tails, heads, weights, demands, adjlist, probtype)
    nedges = length(tails)
    incmx = sparse(tails, heads, collect(1 : nedges))
    facs = [tails[e] for e in 1 : nedges if tails[e] == heads[e]]
    nfacs = length(facs)
    maxnodes = max(maximum(tails), maximum(heads))
    facsid = zeros(Int64, maxnodes)
    for (e, t) in enumerate(facs)
        facsid[t] = e
    end
    custs = [heads[e] for e in 1 : nedges if tails[e] != heads[e]]
    sort!(custs)
    unique!(custs)
    ncusts = length(custs)
    custsid = zeros(Int64, maxnodes)
    for (e, h) in enumerate(custs)
        custsid[h] = e
    end
    caps = zeros(nfacs)
    fcost = zeros(nfacs)
    dems = zeros(ncusts)
    asscosts = -1 * ones(nfacs, ncusts)
    for e in 1 : nedges
        t, h = tails[e], heads[e]
        if t == h
            f = facsid[t]
            if f > 0
                caps[f] = demands[e]
                fcost[f] = weights[e]
            end
        else
            f = facsid[t]
            j = custsid[h]
            if f > 0 && j > 0
                asscosts[f, j] = weights[e]
                dems[j] = demands[e]
            end
        end
    end
    objval, usedfacs, usedass = cflp(nfacs, caps, fcost, ncusts, dems, asscosts, probtype)
    edges = []
    for f in usedfacs
        push!(edges, incmx[facs[f], facs[f]])
    end
    for (t, h) in usedass
        push!(edges, incmx[facs[t], custs[h]])
    end
    objval, edges, []
end

function cflp_abs(xnat)::Function
    fixed_cflp(tails, heads, weights, demands, adjlist, attack) = cflp_(tails, heads, weights, demands, adjlist, xnat)
    return fixed_cflp
end
