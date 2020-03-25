function cflp_orlib_file_parser(filename; cap = -1)
    nfacs = ncusts = 0
    demands = []
    caps = []
    fcosts = []
    asscosts = []
    M = 1e+2
    open(filename) do f
        let
            line = readline(f)
            tok = split(line, " "; keepempty = false)
            nfacs = parse(Int64, tok[1])
            ncusts = parse(Int64, tok[2])
            demands = zeros(Int64, ncusts)
            caps = zeros(Int64, nfacs)
            fcosts = zeros(Int64, nfacs)
            asscosts = zeros(Int64, nfacs, ncusts)
        end
        let
            for i in 1 : nfacs
                line = readline(f)
                tok = split(line, " "; keepempty = false)
                if cap > 0
                    caps[i] = cap
                else
                    caps[i] = parse(Int64, tok[1])
                end
                fcosts[i] = round(Int64, M * parse(Float64, tok[2]))
            end
        end
        let
            for j in 1 : ncusts
                demands[j] = parse(Int64, readline(f))
                ntoks = 0
                while ntoks < nfacs
                    line = readline(f)
                    tok = split(line, " "; keepempty = false)
                    for i in eachindex(tok)
                        asscosts[ntoks + i, j] = round(Int64, M * parse(Float64, tok[i]))
                    end
                    ntoks += length(tok)
                end
            end
        end
    end
    nfacs, ncusts, demands, caps, fcosts, asscosts
end

function cflp_contardosefair_file_reader(filename)
    nfacs = ncusts = 0
    fcosts = []
    caps = []
    attfcosts = []
    demands = []
    asscosts = []
    open(filename) do f
        let
            line = readline(f)
            tok = split(line)
            nfacs = parse(Int64, tok[3])
            ncusts = parse(Int64, tok[4])
            fcosts = zeros(Int64, nfacs)
            attfcosts = zeros(Int64, nfacs)
            caps = zeros(Int64, nfacs)
            demands = zeros(Int64, ncusts)
            asscosts = -1 * ones(Int64, nfacs, ncusts)
        end
        for i in 1 : nfacs
            line = readline(f)
            tok = split(line)
            caps[i] = parse(Int64, tok[2])
            fcosts[i] = parse(Int64, tok[3])
            attfcosts[i] = parse(Int64, tok[4])
        end
        for n in 1 : ncusts
            line = readline(f)
            tok = split(line)
            demands[n] = parse(Int64, tok[2])
        end
        for i in 1 : nfacs, j in 1 : ncusts
            line = readline(f)
            tok = split(line)
            t, h = parse(Int64, tok[2]), parse(Int64, tok[3])
            asscosts[t, h] = parse(Int64, tok[4])
        end
    end
    nfacs, ncusts, demands, caps, fcosts, attfcosts, asscosts
end

function cflp_boccia_file_parser(filename)
    nfacs = ncusts = 0
    demands = []
    caps = []
    fcosts = []
    asscosts = []
    M = 1e+2
    open(filename) do f
        let
            line = readline(f)
            tok = split(line, "\t")
            ncusts = parse(Int64, tok[1])
            nfacs = parse(Int64, tok[2])
            demands = zeros(Int64, ncusts)
            caps = zeros(Int64, nfacs)
            fcosts = zeros(Int64, nfacs)
            asscosts = zeros(Int64, nfacs, ncusts)
        end
        let
            ntoks = 0
            while ntoks < ncusts
                line = readline(f)
                tok = split(line, "\t"; keepempty = false)
                for i in eachindex(tok)
                    demands[ntoks + i] = parse(Int64, tok[i])
                end
                ntoks += length(tok)
            end
        end
        let
            ntoks = 0
            while ntoks < nfacs
                line = readline(f)
                tok = split(line, "\t"; keepempty = false)
                for i in eachindex(tok)
                    caps[ntoks + i] = parse(Int64, tok[i])
                end
                ntoks += length(tok)
            end
        end
        let
            ntoks = 0
            while ntoks < nfacs
                line = readline(f)
                tok = split(line, "\t"; keepempty = false)
                for i in eachindex(tok)
                    fcosts[ntoks + i] = round(Int64, M * parse(Float64, tok[i]))
                end
                ntoks += length(tok)
            end
        end
        let
            i = j = 1
            while i <= nfacs
                line = readline(f)
                tok = split(line, "\t"; keepempty = false)
                for k in eachindex(tok)
                    asscosts[i, j] = round(Int64, M * parse(Float64, tok[k]))
                    j += 1
                    if j > ncusts
                        j = 1
                        i += 1
                    end
                end
            end
        end
    end
    nfacs, ncusts, demands, caps, fcosts, asscosts
end

function spp_file_reader(filename)
    filepath = filename#"./instances/shortestpath/data/$filename"
    tails = []
    heads = []
    weights = []
    weightsp = []
    nnodes = nedges = 0
    s = t = 0
    open(filepath) do f
        for line in eachline(f)
            token = split(line, " ")
            lt = token[1]
            if lt == "p"
                nnodes = parse(Int64, token[3])
                nedges = parse(Int64, token[4])
            elseif lt == "n"
                if token[3] == "s"
                    s = parse(Int64, token[2]) + 1
                elseif token[3] == "t"
                    t = parse(Int64, token[2]) + 1
                end
            elseif lt == "a"
                tail = parse(Int64, token[2]) + 1
                head = parse(Int64, token[3]) + 1
                weight = parse(Int64, token[4])
                attack = parse(Int64, token[5])
                push!(tails, tail)
                push!(heads, head)
                push!(weights, weight)
                push!(weightsp, attack)
            end
        end
    end
    tails, heads, weights, weightsp, s, t
end

function spp_road_file_reader(filename)
    filepath = filename
    nnodes = nedges = 0
    s = t = 0
    f = open(filename)
    let
        line = readline(f)
        token = split(line, " ")
        nnodes = parse(Int64, token[3])
        nedges = parse(Int64, token[4])
    end
    let
        line = readline(f)
        token = split(line, " ")
        s = parse(Int64, token[1])
        line = readline(f)
        token = split(line, " ")
        t = parse(Int64, token[1])
    end

    addedTSPEdge = falses(nnodes)

    tails = zeros(Int64, 0)
    heads = zeros(Int64, 0)
    weights = zeros(Int64, 0)
    weightsp = zeros(Int64, 0)
    for e in 1 : nedges
        line = readline(f)
        token = split(line, " ")
	tail = parse(Int64, token[1])
        head = parse(Int64, token[2])
        if (head - tail == 1 || (tail == nnodes && head == 1))
            if addedTSPEdge[tail]
                continue
            else
                addedTSPEdge[tail] = true
            end
        end
        push!(tails, tail)
        push!(heads, head)
        push!(weights, parse(Int64, token[3]))
        push!(weightsp, parse(Int64, token[4]))
    end
    close(f)
    tails, heads, weights, weightsp, s, t
end

function mst_file_reader(filename)
    filepath = filename#"./instances/shortestpath/data/$filename"
    tails = []
    heads = []
    weights = []
    weightsp = []
    nnodes = nedges = 0
    open(filepath) do f
        for line in eachline(f)
            token = split(line, " ")
            lt = token[1]
            if lt == "p"
                nnodes = parse(Int64, token[3])
                nedges = parse(Int64, token[4])
            elseif lt == "e"
                tail = parse(Int64, token[2]) + 1
                head = parse(Int64, token[3]) + 1
                weight = parse(Int64, token[4])
                attack = parse(Int64, token[5])
                push!(tails, tail)
                push!(heads, head)
                push!(weights, weight)
                push!(weightsp, attack)
            end
        end
    end
    tails, heads, weights, weightsp
end

function knapsack_bkpinsv2_file_reader(filename)
    nnodes = 0
    followerdems = Int64[]
    vals = Int64[]
    cap = 0
    budget = 0
    open(filename) do f
        ln = 0
        for line in eachline(f)
            ln += 1
            token = split(line, " "; keep = false)
            if ln == 1 continue
            elseif ln == 2
                nnodes = parse(Int64, token[1])
            elseif ln == 3
                budget = parse(Int64, token[1])
            elseif ln == 4
                cap = parse(Int64, token[1])
            elseif ln == 5
                for t in token
                    if t == "" continue
                    else
                        push!(vals, parse(Int64, t))
                    end
                end
            elseif ln == 6
                for t in token
                    if t == "" continue
                    else
                        push!(followerdems, parse(Int64, t))
                    end
                end
            end
        end
    end
    leaderdems = ones(Int64, nnodes)
    attacked = copy(vals)
    vals, attacked, leaderdems, followerdems, cap, budget
end

function knapsack_mibs_file_reader(mpsfilename, auxfilename)
    nnodes = 0
    demands = Int64[]
    vals = Int64[]
    cap = 0
    budget = 0
    open(auxfilename) do f
        for line in eachline(f)
            token = split(line, " "; keep = false)
            lt = token[1]
            if lt == "N"
                nnodes = parse(Int64, token[2])
            elseif lt == "LO"
                weight = -parse(Int64, token[2])
                push!(vals, weight)
            end
        end
    end
    open(mpsfilename) do f
        for line in eachline(f)
            token = split(line, " "; keep = false)
            println(token)
            lt = token[1]
            if lt[1] == 'y' && token[2]  == "KF"
                push!(demands, parse(Int64, token[3]))
            elseif lt == "rhs" && token[2] == "KL"
                budget = parse(Int64, token[3])
            elseif lt == "rhs" && token[2] == "KF"
                cap = parse(Int64, token[3])
            end
        end
    end
    attacked = copy(vals)
    vals, attacked, demands, cap, budget
end

function knapsack_coral_file_reader(mpsfilename, auxfilename)
    nnodes = 0
    leadercons = Int64[]
    followercons = Int64[]
    vals = Int64[]
    cap = 0
    budget = 0
    open(auxfilename) do f
        for line in eachline(f)
            token = split(line, [' ', '.']; keep = false)
            lt = token[1]
            if lt == "IC"
                v = parse(Int64, token[2])
                push!(followercons, v)
            elseif lt == "IB"
                cap = parse(Int64, token[2])
            elseif lt == "LO"
                push!(vals, -parse(Int64, token[2]))
            end
        end
    end
    open(mpsfilename) do f
        for line in eachline(f)
            token = split(line, [' ', '.']; keep = false)
            println(token)
            lt = token[1]
            if length(lt) >= 3 && SubString(lt, 1, 3) == "C00"
                push!(leadercons, parse(Int64, token[5]))
            elseif lt == "RHS" && length(token) >= 3
                budget = parse(Int64, token[3])
            end
        end
    end
    attacked = copy(vals)
    vals, attacked, leadercons, followercons, cap, budget
end

function knapsack_file_reader(filename)
    filepath = filename#"./instances/shortestpath/data/$filename"
    vals = Int64[]
    attacked = Int64[]
    followercons = Int64[]
    leadercons = Int64[]
    nnodes = 0
    cap = budget = 0
    open(filepath) do f
        for line in eachline(f)
            token = split(line, " ")
            lt = token[1]
            if lt == "p"
                nnodes = parse(Int64, token[3])
            elseif lt == "c"
                if token[2] == "f"
                    cap = parse(Int64, token[3])
                elseif token[2] == "l"
                    budget = parse(Int64, token[3])
                end
            elseif lt == "n"
                v = parse(Int64, token[2])
                a = parse(Int64, token[3])
                fc = parse(Int64, token[4])
                lc = parse(Int64, token[5])
                push!(vals, v)
                push!(attacked, a)
                push!(followercons, fc)
                push!(leadercons, lc)
            end
        end
    end
    vals, attacked, followercons, leadercons, cap, budget
end

function perfect_matching_file_reader(filename)
    filepath = filename#"./instances/shortestpath/data/$filename"
    tails = []
    heads = []
    weights = []
    attacked = []
    nnodes = 0
    nedges = 0
    open(filepath) do f
        for line in eachline(f)
            token = split(line, " ")
            lt = token[1]
            if lt == "p"
                nnodes = parse(Int64, token[3])
                nedges = parse(Int64, token[4])
            elseif lt == "e"
                t = parse(Int64, token[2]) + 1
                h = parse(Int64, token[3]) + 1
                c = parse(Int64, token[4])
                d = parse(Int64, token[5])
                push!(tails, t)
                push!(heads, h)
                push!(weights, c)
                push!(attacked, d)
            end
        end
    end
    tails, heads, weights, attacked
end
