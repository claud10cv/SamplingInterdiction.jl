using Random

mutable struct Tree
    edges::Array{Int64, 1}
    cost::Number
end

mutable struct Attack
    attacked::Array{Int64, 1}
    mst::Array{Int64, 1}
    cost::Number
    ub::Number
end

struct AttackCut
    lambda_coef::Number
    edge_coefs::Array{Number, 1}
    rhs::Number
end

struct Label
    edge_status::Array{Int64, 1} # 0 = free; 1 = attacked; 2 = cannot attack
    optedges::Array{Int64, 1}
    parents::Array{Int64, 1}
    depth::Int64
    lowerbound::Int64
    freenodes::Set
    numfree::Int64
    attacked::Array{Int64, 1}
    numattacked::Int64
    ubattack::Int64
end

struct FortificationLabel
    edge_status::Array{Int64, 1}
    optedges::Array{Int64, 1}
    depth::Int64
    upperbound::Int64
    freenodes::Set
    numfree::Int64
    defended::Array{Int64, 1}
    numdefended::Int64
end

mutable struct Params
    rng
    heuristic_restricted_interdiction
end

gParams = Params(MersenneTwister(1), true)
