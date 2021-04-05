using Pkg

using CSV
using DataFrames
using ForwardDiff
using LinearAlgebra
using Optim

using Random, Distributions
using Statistics
using StatsPlots

using Plots

##

# x = LinRange(1,10,100)
# y = x.^2
# println("1")
# println(y.^2)


# plot(x,y.^6)

cache = Dict()
cache[1] = 1
cache[2] = 2
function fib(n)
    # println("fib on ", n, " called")
    if n in keys(cache)
        return cache[n]
    end
    cache[n] = fib(n-1) + fib(n-2)
    return cache[n]
end

# println(fib(1))
# println(fib(2))
# println(fib(3))
# println(fib(4))
println(fib(800))