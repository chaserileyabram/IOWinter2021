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
rho_0 = 0.05
rho_1 = 0.03
gamma = 1
epsilon = 1.1
theta = 50
varphi = 1
nu = 0.5
kappa = (epsilon-1)/theta*(gamma+nu)

# phi_pi spiral bound
phi_pi_sp = 1 + rho^2/(4*kappa)


# Policy
@with_kw struct Policyparam @deftype Float64
    # MP
    i_star = rho_0 # Should probably make match natural real
    phi_pi = 0.9
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

# pi dot
function dpi(pi,x)
    return rho*pi - kappa*x
end

# x dot
function dx(pi,x,b,pp)
    @unpack_Policyparam pp

    return i_star + phi_pi*pi + phi_x*x + phi_b*b - pi - rho
end

# b dot
function db(pi,x,b,pp)
    @unpack_Policyparam pp

    return (i_star + phi_pi*pi + phi_x*x + phi_b*b - pi)*b - (tau_star + psi_pi*pi + psi_x*x + psi_b*b)
end

T = 200
# Get time paths
pis = ones(T+1)
xs = ones(T+1)
bs = ones(T+1)

pi_0 = 0
x_0 = 0
b_0 = 1

pis[1] = pi_0
xs[1] = x_0
bs[1] = b_0

dt = 0.1

policy_p = Policyparam()

# Run IRFs
for t = 2:T+1
    pis[t] = pis[t-1] + dt*dpi(pis[t-1],xs[t-1])
    xs[t] = xs[t-1] + dt*dx(pis[t-1],xs[t-1],bs[t-1],policy_p)
    bs[t] = bs[t-1] + dt*db(pis[t-1],xs[t-1],bs[t-1],policy_p)
end

# println("pis: ", pis)
# println("---")
# println("xs: ", xs)
# println("---")
# println("bs: ", bs)

ts = LinRange(0,T,T+1)
tvc_term = exp.(-rho *ts) .* bs

pi_plot = plot(ts, pis, title = "pi")
x_plot = plot(ts, xs, title = "x")
b_plot = plot(ts, bs, title = "b", label = "b")
tvc_plot = plot(ts, tvc_term, title = "TVC", label = "TVC term")

display(pi_plot)
display(x_plot)
display(b_plot)
display(tvc_plot)

println("phi_pi: ", policy_p.phi_pi)
println("phi_pi_sp: ", phi_pi_sp)
if policy_p.phi_pi > phi_pi_sp
    println("Spiral time!")
end









