# SamplingInterdiction.jl
Sampling for interdiction problems

## Usage

### 0-1 Knapsack interdiction
```julia
v = rand(1 : 100, 100)
vp = [rand(1 : x) for x in v]
lw = rand(1 : 100, 100)
fw = rand(1 : 100, 100)
lK = rand(100 : 200)
fK = rand(200 : 500)
report, opt = SamplingInterdiction.knapsack_interdiction_sampling(v, vp, lw, fw, fK, lK)
```

### Shortest path interdiction
```julia
tails = [1, 1, 2, 3, 3, 4, 4, 4, 5, 5]
heads = [2, 3, 3, 4, 2, 2, 3, 4, 3, 4]
cost = rand(1 : 5, 10)
dcost = rand(1 : 3, 10)
lw = rand(1 : 3, 10)
budget = 2
s, t = 1, 5
report, opt = SamplingInterdiction.sp_interdiction_sampling(tails, heads, cost, dcost, lw, budget, s, t)
```

### Facility location interdiction
```julia
nat = :uflp # use :sscflp for single-source facility location with interdiction
caps = rand(100 : 300, 10)
dems = rand(1 : 10, 100)
fcosts = rand(100 : 200, 10)
dfcosts = rand(20 : 50, 10)
asscosts = rand(5 : 10, 10, 100)
budget = 2
report, opt = SamplingInterdiction.cflp_interdiction_sampling(nat, caps, dems, fcosts, dfcosts, asscosts, budget)
```

