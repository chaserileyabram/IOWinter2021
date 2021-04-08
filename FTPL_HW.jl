using Pkg
# Pkg.add("DifferentialEquations")

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

using DifferentialEquations

##

Random.seed!(1234)

# Parameters (need to check against a paper)
rho_0 = 0.05
rho_1 = 0.04
gamma = 1.0
epsilon = 10.0
theta = 100.0
varphi = 2.2 # currently not used
nu = 1.0
kappa = (epsilon-1)/theta*(gamma+nu)

# phi_pi spiral bound
phi_pi_sp = 1 + rho_1^2/(4*kappa)


# Policy
@with_kw struct Policyparam @deftype Float64
    # MP
    i_star = rho_0 # Should probably make match natural real
    phi_pi = 0.0 # 1.25
    phi_b = 0.0
    phi_x = 0.0
    # FP
    tau_star = 1.0
    psi_pi = 0.0
    psi_b = 2.0
    psi_x = 0.0
end


# Utility
# function u(c, g) 
#     if g == 1
#         return log(c)
#     else
#         return (c^(1-g) - 1)/(1-g)
#     end
# end

# # Marginal Utility
# function up(c,g)
#     if g == 1
#         return 1/c
#     else
#         return c^(-g)
#     end
# end

# # pi dot
# function dpi(pi,x)
#     return rho_1*pi - kappa*x
# end

# # x dot
# function dx(pi,x,b,pp)
#     @unpack_Policyparam pp

#     return (i_star + phi_pi*pi + phi_x*x + phi_b*b - pi - rho_1)/gamma
# end

# # b dot
# function db(pi,x,b,pp)
#     @unpack_Policyparam pp

#     return (i_star + phi_pi*pi + phi_x*x + phi_b*b - pi)*b - (tau_star + psi_pi*pi + psi_x*x + psi_b*b)
# end

# T = 2000
# # Get time paths
# pis = ones(T+1)
# xs = ones(T+1)
# bs = ones(T+1)

# pi_0 = 0
# x_0 = 0
# b_0 = 0

# pis[1] = pi_0
# xs[1] = x_0
# bs[1] = b_0

# dt = 0.1

# policy_p = Policyparam()

# # Run IRFs
# for t = 2:T+1
#     # if t > 500
#     #     policy_p = Policyparam(phi_pi = 0.9)
#     # end
#     pis[t] = pis[t-1] + dt*dpi(pis[t-1],xs[t-1])
#     xs[t] = xs[t-1] + dt*dx(pis[t-1],xs[t-1],bs[t-1],policy_p)
#     bs[t] = bs[t-1] + dt*db(pis[t-1],xs[t-1],bs[t-1],policy_p)
# end

# println("pis: ", pis)
# println("---")
# println("xs: ", xs)
# println("---")
# println("bs: ", bs)

# ts = LinRange(0,T,T+1)
# tvc_term = exp.(-rho_1 *ts) .* bs

# pi_plot = plot(ts, pis, title = "pi")
# x_plot = plot(ts, xs, title = "x")
# b_plot = plot(ts, bs, title = "b", label = "b")
# tvc_plot = plot(ts, tvc_term, title = "TVC", label = "TVC term")

# display(pi_plot)
# display(x_plot)
# display(b_plot)
# display(tvc_plot)

# println("phi_pi: ", policy_p.phi_pi)
# println("phi_pi_sp: ", phi_pi_sp)
# if policy_p.phi_pi > phi_pi_sp
#     println("Spiral time!")
# end

mu = 0.5
pp = Policyparam()

function nkftpl!(du,u,p,t)

    rho(s) = rho_0 + exp(-mu*s)*(rho_1 - rho_0)

    # pi
    du[1] = rho(t)*u[1] - kappa*u[2]

    # x
    du[2] = 1/gamma*(pp.i_star + pp.phi_pi*u[1] - u[1] - rho(t))

    # b
    du[3] = (pp.i_star + pp.phi_pi*u[1] - u[1])*u[3] - pp.tau_star - pp.psi_b*u[3]
end

b_0 = 0.0
b_T = 0.0
x_0 = 0.0

# Need to play with these
function bc!(residual, u, p, t)
    # inital debt
    # residual[4] = u[1][3] - b_0

    # non-explosive output gap
    residual[1] = u[end][2]

    # final debt
    residual[2] = u[end][3] - b_T

    # initial output gap
    # residual[3] = u[1][2] - x_0

    # non-explosive inflation
    residual[3] = u[end][1]
end

T = 20
tspan = (0.0, T)
tspan_range = LinRange(0,T,T+1)

# u0 = [0.0; 0.0; 1.0]
# prob = ODEProblem(nkftpl!, u0, tspan)
# sol = solve(prob)

# p_pi = plot(sol, vars = (1), title = "pi")
# display(p_pi)
# p_x = plot(sol, vars = (2), title = "x")
# display(p_x)
# p_b = plot(sol, vars = (3), title = "b")
# display(p_b)

bvp = BVProblem(nkftpl!, bc!, [0,0,0], tspan)
sol_bc = solve(bvp, GeneralMIRK4(), dt=0.5)

p_pi_bc = plot(sol_bc, vars = (1), title = "pi_bc")
display(p_pi_bc)
p_x_bc = plot(sol_bc, vars = (2), title = "x_bc")
display(p_x_bc)
p_b_bc = plot(sol_bc, vars = (3), title = "b_bc")
display(p_b_bc)
println(sol_bc)

p_pi_x_bc = plot(sol_bc, vars = (1,2), 
title = "pi_x_bc",
xlabel = "pi", ylabel = "x",
arrow = true)
display(p_pi_x_bc)

# p_all_bc  = plot(sol_bc, vars = (1,2,3), title = "all_bc")
# display(p_all_bc)

