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
# Pkg.add("Parameters")
using Parameters

##

# Parameters (need to check against a paper)
rho = 0.05
gamma = 1
epsilon = 1.1
theta = 5
varphi = 1
nu = 0.5


# Policy
@with_kw struct Policyparam
    @deftype Float64
    # MP
    i_star = 0.05 # Should probably make match natural real
    phi_pi = 1.1
    phi_b = 0.0
    phi_x = 0.0
    # FP
    tau_star = 0
    psi_pi = 0.0
    psi_b = 0.0
    psi_x = 0.0
end


# Utility
function u(c, g) 
    if g == 1
        return log(c)
    else
        return (c^(1-g) - 1)/(1-g)
    end
end

# Marginal Utility
function up(c,g)
    if g == 1
        return 1/c
    else
        return c^(-g)
    end
end

function dpi(p)












