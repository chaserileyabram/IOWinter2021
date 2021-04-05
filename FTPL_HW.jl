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


function fib(n)
    println("on ", n, "called")
    if n == 1
        return 1
    elseif n == 2
        return 2
    end

    return fib(n-1) + fib(n-2)
end


