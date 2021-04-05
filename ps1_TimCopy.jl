# Chase Abram
using Pkg
# Pkg.add("DataFrames")
using CSV
using DataFrames
using ForwardDiff
using LinearAlgebra
using Optim

using Random

# Load data
df = DataFrame()
df = CSV.read("psetOne.csv", DataFrame)

df = CSV.File("psetOne.csv") |> DataFrame

# names(df)

# Get market t = 17
df17 = df[df."Market" .== 17, :]

# Store variables for market 17
X = convert(Matrix, df17[:, filter(x -> x in ["Constant", "EngineSize", "SportsBike", "Brand2", "Brand3"], names(df17))])
Br = convert(Matrix, df17[:, filter(x -> x in ["Brand2", "Brand3"], names(df17))])
Br = Bool.(Br)
p_given = df17."Price"
df17


# Shares
function s_logit(p, alpha, beta, X, xi)
    
    # Compute numerators for each s_j
    nums = exp.(X*beta - alpha.*p + xi)
    
    return nums./(1 + sum(nums))
end

# Jacobian of shares
function s_logit_jac(p, alpha, beta, X, xi)
    # Get shares
    s = s_logit(p, alpha, beta, X, xi)
    
    # Jacobian is the outer product times alpha in this case...
    jac = alpha .* s * s'
    
    # ... except for an extra term on the diag
    for j in 1: length(p)
        jac[j,j] -= alpha*s[j]
    end
    
    return jac
end

# Hessian of shares (n by n by n tensor)
function s_logit_hess(p, alpha, beta, X, xi)
    
    # Initialize
    hess = zeros(length(p), length(p), length(p))
    
    # Shares
    s = s_logit(p, alpha, beta, X, xi)

    # Jacobian of shares
    sj = s_logit_jac(p, alpha, beta, X, xi)
    
    # Fill in
    for j in 1:length(p)
        for k in 1:length(p)
            for l in 1:length(p)
                
                if j == k
                    # on diag of jacobian
                    hess[j,k,l] = alpha*sj[l,j]*(2*s[j] - 1)
                else
                    # off diag of jacobian
                    hess[j,k,l] = alpha*(sj[l,j]*s[k] + s[j]*sj[l,k])
                end
            end
        end
    end
    
    return hess
end

# Ownership structure
function Omega_star(br)
    
    # Initialize
    Om = zeros(size(br,1),size(br,1))
    
    for i in 1:size(Om,1)
        for j in 1:size(Om,2)
            if br[i,1] && br[j,1]
                # Brand2
                Om[i,j] = 1
            elseif br[i,2] && br[j,2]
                # Brand3
                Om[i,j] = 1
            elseif (!br[i,1]&& !br[j,1]) && (!br[i,2]&& !br[j,2])
                # Brand1
                Om[i,j] = 1
            end
        end
    end
    
    return Om
end

# Omega in Nevo (2001)
function Omega(p, alpha, beta, X, xi, br)
    return Omega_star(br) .* -s_logit_jac(p, alpha, beta, X, xi)
end

# Fixed point expression (want zero at FP)
function fp_logit(p, alpha, beta, X, xi, br)
    return Omega(p, alpha, beta, X, xi, br) \ s_logit(p, alpha, beta, X, xi) 
    #s_logit(p, alpha, beta, X, xi) - Omega(p, alpha, beta, X, xi, br)*p
end

# (Jacobian of Omega)*p (this is used in the gradient of the FP expression)
function jac_om_p(p, alpha, beta, X, xi, br)
    # Initialize
    jacomp = zeros(length(p), length(p))
    
    # Jacobian of shares
    sj = s_logit_jac(p, alpha, beta, X, xi)
    
    # Hessian of shares
    sh = s_logit_hess(p, alpha, beta, X, xi)
    
    # Ownership
    oms = Omega_star(br)
    
    # Fill in
    for i in 1:length(p)
        for k in 1:length(p)
            jacomp[i,k] = -(sh[i,:,k]' * (oms[i,:] .* p))[1]
        end
    end
    return jacomp
end

# Jacobian of fixed point expression
function J_fp_logit(p, alpha, beta, X, xi, br)
    
    # Jacobian of shares
    sj = s_logit_jac(p, alpha, beta, X, xi)
    
    # (Jacobian of Omega)*p
    jacomp = jac_om_p(p, alpha, beta, X, xi, br)
    
    # Omega
    om = Omega(p, alpha, beta, X, xi, br)
#     println("sj: ", sj)
#     println("jacomp: ", jacomp)
#     println("om: ", om)
    return sj - jacomp - om
end
using NLsolve
# Solve for fixed point
function solve_fp_logit(p0, alpha, beta, X, xi, br, tol = 1e-14, maxiter = 1000)
    # Initialize price and update storage
    p = zeros(length(p0),1)
    pnew = zeros(length(p0),1)
    
    p .= p0
    
    # Initialize Jacobian
    J = zeros(length(p0), length(p0))
    
    # Iterator
    it = 0
    
    # Difference in updating
    diff = Inf
    println("p0: ", p0)

    # # This is a comparison to check I'm not being a moron, both converge
    # fixedpoint(x->fp_logit(x, alpha, beta, X, xi, br),
    #            ones(7)*0.0, show_trace = true, method = :newton )
    
    # Continue until converged or maximum iterations
    while diff > tol && it < maxiter
#         println("p: ", p)
        
        # Jacobian of FP expression
        #J .= J_fp_logit(p, alpha, beta, X, xi, br)

        J = ForwardDiff.jacobian( x->fp_logit(x, alpha, beta, X, xi, br) - x, p )
#         println("J: ", J)
#         println("J inv: ", inv(J))
#         println("J under fp: ", J \ fp_logit(p, alpha, beta, X, xi, br))
#         println("size(p): ", size(p))
#         println("size(J under fp): ", size(J \ fp_logit(x, alpha, beta, X, xi, br)))
        
        # Update
        pnew = p - J \ (fp_logit(p, alpha, beta, X, xi, br) - p)
        it += 1
        
#         if it > 20
#             println("s: ", s_logit(p, alpha, beta, X, xi))
#             println("J: ", J)
#             println("det(J): ", det(J))
#             println("fp_logit: ", fp_logit(p, alpha, beta, X, xi, br))
#         end
        
        println("pnew: ", pnew)
        println("Function Value: ", maximum(abs.(fp_logit(pnew, alpha, beta, X, xi, br)-p)) )
        
        # Get difference
        diff = maximum(abs.(pnew .- p))
        println("diff: ", diff)
        
        p = pnew
    end
    
    println("Took ", it, " iterations")
    println("final fp_logit: ", fp_logit(p, alpha, beta, X, xi, br))
    return p
end

#p_init = ones(size(X,1),1) .* 1.0
p_init = p_given
# p .= p_given
# p[7] = 1
alpha = 3
beta = [1 1 2 -1 1]'


# seed maintains results across runs
rng = MersenneTwister(1234);
xi = randn(rng, size(df17,1)) .* 0.0;
# println("xi: ", xi)

# println("p: ", p)
# println("delta: ", X*beta - alpha .* p + xi)

s_logit(p_init, alpha, beta, X, xi)
println("shares: ", s_logit(p_init, alpha, beta, X, xi))
# println("sum shares: ", sum(s_logit(p, alpha, beta, X, xi)))

s_logit_jac(p_init, alpha, beta, X, xi)
# println("shares_jac: ", s_logit_jac(p, alpha, beta, X, xi))

# s_logit_hess(p, alpha, beta, X, xi)

# Omega_star(Br)

# Omega(p, alpha, beta, X, xi, Br)

# fp_logit(p, alpha, beta, X, xi, Br)' * fp_logit(p, alpha, beta, X, xi, Br)
# println("fp_logit: ", fp_logit(p, alpha, beta, X, xi, Br))
# fp_logit(p, alpha, beta, X, xi, Br)

# J_fp_logit(p, alpha, beta, X, xi, Br)

out = solve_fp_logit(p_init, alpha, beta, X, xi, Br)
# println("out: ", out)
# println("eval fp_logit: ", fp_logit(out, alpha, beta, X, xi, Br))

# op = optimize(x -> (fp_logit(x, alpha, beta, X, xi, Br)' * fp_logit(x, alpha, beta, X, xi, Br))[1,1], p, LBFGS())

# pmin = op.minimizer
# println("p min: ", pmin)
# println("min: ", op.minimum)
# println("fp_logit: ", fp_logit(pmin, alpha, beta, X, xi, Br))

# NM: opam: [0.3720742469890362; 0.3724406061806239; 0.3447723065789802; 
# 0.3424326539927414; 2.4893064284464703; 2.5105275013368717; 4.366510826679992]

# BFGS: opam: [0.37346117252339095; 0.37346116923202516; 5.8690207005673365; 
# 0.34160603728869476; 2.5287477985474984; 2.5287470025674392; 2.528747016448347]

# From Sam: 
# 0.3938563345938473
# 0.3938563345938473
# 0.3496634377039666
# 0.3496634413963176
# 1.6632109730370337
# 1.6632109730370337
# 1.6632541435947288





# function f1(x, y, z) 
#     return (x[1] - 5 - y)^2 + (x[2] - 3 - z)^2
# end

# function g1(storage, x, y, z)
#     storage[1] = 2*(x[1] - 5 - y)
#     storage[2] = 2*(x[2] - 3 - z)
# end

# function g2(x,y,z)
#     out = zeros(2)
#     out[1] = 2*(x[1] - 5 - y)
#     out[2] = 2*(x[2] - 3 - z)
    
#     return out
# end

# function tester()
#     a = 1000
#     b = -1000
#     f(x) = f1(x,a,b)
# #     g!(storage, x) = g1(storage, x, a, b)
#     g!(x) = g2(x,a,b)
    
#     out = basic_newton_zero(f, g!, [-200.0, 1000.0])
    
# #     opt = optimize(f, g!, [100.0, -20.0], LBFGS())
# #     return opt, opt.minimizer, opt.minimum
#     return out
# end

# tester()



p_init = ones(size(X,1),1) .* 1.0
# p_init[1] = 50
# p_init[2] = 30
# p_init[3] = 90
# p_init[4] = -12
# p_init[5] = -40
# p_init[6] = 0
# p_init[7] = 1
s(x) = s_logit(x, alpha, beta, X, xi)
# println(s_logit_jac(p, alpha, beta, X, xi))
ForwardDiff.jacobian(s, p_init)

# ForwardDiff.gradient(fp, p)


# fp(x) = fp_logit_2(x, alpha, beta, X, xi, Br)
# println(fp_logit_2_jac(p, alpha, beta, X, xi, Br))
# ForwardDiff.gradient(fp, p)

fp(x) = fp_logit(x, alpha, beta, X, xi, Br)
J_ex = J_fp_logit(p_init, alpha, beta, X, xi, Br)
println(J_ex)
# println(maximum(abs.(ForwardDiff.jacobian(fp, p) .- J_ex)))
ForwardDiff.jacobian(fp, p_init)
# println(inv(J_ex))

# Compute probabilities of agent i choosing product j
function pr(delta::AbstractVector{T}, X, sigma, zeta) where T
    
    # Initialize utilities
    ubar = zeros(T, size(zeta,2), size(delta,1))
    
    # Compute utilities
    for i in 1:size(ubar, 1)
        for j in 1:size(ubar, 2)
            ubar[i,j] = delta[j] + X[j,:]' * sigma * zeta[:,i]
        end
    end
    
    # Initialize probabilities
    p = zeros(T, size(ubar,1), size(ubar,2))
    
    # Compute probabilities
    for i in 1:size(p,1)
        for j in 1:size(p,2)
            p[i,j] = exp(ubar[i,j])/(1 + sum(exp.(ubar[i,:])))
        end
    end
    
    return p
end



# Shares
function sHat(delta::AbstractVector{T}, X, sigma, zeta) where T
    
    # Compute probs
    p = pr(delta, X, sigma, zeta)
    
    # Initialize shares
    s = zeros(T, size(delta))
    
    # Compute share
    for j in 1:length(s)
        s[j] = 1/size(zeta,2) * sum(p[:,j])
    end
    
    return s
end


# Test cases
nI = 20
nJ = 3
nN = 5

delta_1 = zeros(3)
X_1 = zeros(3,nN)
sigma_1 = 0.0 .* I(size(X,2))
zeta_1 = zeros(nN, nI)

delta_2 = zeros(3)
delta_2[1] = 40
delta_2[2] = 20
delta_2[3] = 20
X_2 = zeros(nJ,nN)
sigma_2 = 0.0 .* I(size(X,2)).* 0.1
zeta_2 = zeros(nN, nI)

sigma_3 = 0.1 .* I(nN)
delta_3 = zeros(nJ)
X_3a = zeros(nJ,nN)

rng = MersenneTwister(1234)
zeta_3 = randn(rng, nN, nI)

rng = MersenneTwister(1234)
X_3b = randn(rng, size(X_3a))
X_3c = X_3b .* 10
X_3d = abs.(X_3c)


s1 = sHat(delta_1, X_1, sigma_1, zeta_1)
s2 = sHat(delta_2, X_2, sigma_2, zeta_2)
s3a = sHat(delta_3, X_3a, sigma_3, zeta_3)
s3b = sHat(delta_3, X_3b, sigma_3, zeta_3)
s3c = sHat(delta_3, X_3c, sigma_3, zeta_3)
s3d = sHat(delta_3, X_3d, sigma_3, zeta_3)

println("s1: ", s1)
println("s2: ", s2)
println("s3a: ", s3a)
println("s3b: ", s3b)
println("s3c: ", s3c)
println("s3d: ", s3d)

# Jacobian of shares
function sHat_jac(delta, X, sigma, zeta)
    # Get probabilities
    p = pr(delta, X, sigma, zeta)
    
    # Get shares
    s = sHat(delta, X, sigma, zeta)
    
    # Initialize Jacobian
    sj = zeros(size(delta,1), size(delta,1))
    
    # Fill in
    for j in 1:size(sj,1)
        for k in 1:size(sj,2)
            if j == k
                sj[j,k] = 1/size(zeta,2) * sum(p[:,j] .* (1 .- p[:,j]))
            else
                sj[j,k] = 1/size(zeta,2) * sum(-p[:, j] .* p[:, k])
            end
        end
    end
    
    return sj
end

# Jacobian of log of shares
function sHat_log_jac(delta, X, sigma, zeta)
    # shares
    s = sHat(delta, X, sigma, zeta)
    
    # jacobian of shares
    sj = sHat_jac(delta, X, sigma, zeta)
    
    # initialize jacobian of log of share
    slj = zeros(size(delta,1), size(delta,1))
    
    # Fill in
    for j in 1:size(slj, 1)
        for k in 1:size(slj,2)
            slj[j,k] = sj[j,k]/s[k]
        end
    end
    
    return slj
end

s1_jac = sHat_jac(delta_1, X_1, sigma_1, zeta_1)
println("s1_jac: ", s1_jac)
s1_log_jac = sHat_log_jac(delta_1, X_1, sigma_1, zeta_1)
println("s1_log_jac: ", s1_log_jac)

function sHat_inv(s, X, sigma, zeta, tol = 1e-14, maxiter = 1e5)
    
    # Initialize
    delta = zeros(size(X,1))
    shat = zeros(size(delta))
    
    diff = maximum(log.(s) - log.(sHat(delta, X, sigma, zeta)))
    it = 0
    
    while diff > tol && it < maxiter

        if it % 1000 == 0
            println(diff)
        end
        
#         println("it: ", it)
        
        shat = sHat(delta, X, sigma, zeta)
#         println("shat: ", shat)
        
#         if diff > 1e-3
#             println("Contraction")
            inc = log.(s) - log.(shat)
#         else
#             println("Newton")
            
#             slj = sHat_log_jac(delta, X, sigma, zeta)
# #             println("slj: ", slj)
#             inc = -slj \ log.(shat)
#         end
        
        delta += inc
        
        diff = maximum(abs.(inc))
#         println("diff: ", diff)
        it += 1
    end
    println("Loop exited")
    
    println("it: ", it)
    println("diff: ", diff)
    
    
    return delta
end

println("s1_inv: ", sHat_inv(s1, X_1, sigma_1, zeta_1))
println("s2_inv: ", sHat_inv(s2, X_2, sigma_2, zeta_2))
println("s3a_inv: ", sHat_inv(s3a, X_3a, sigma_3, zeta_3))
# println("s3b_inv: ", sHat_inv(s3b, X_3b, sigma_3, zeta_3))
# println("s3c_inv: ", sHat_inv(s3c, X_3c, sigma_3, zeta_3))
# println("s3d_inv: ", sHat_inv(s3d, X_3d, sigma_3, zeta_3))

delta_1 = ones(3)
delta_2 = 2
delta_3 = 40
X_1 = 10 .* randn(3,nN)
sigma_1 = 0.1 .* I(size(X,2))
zeta_1 = randn(nN, 10000)

shl(x) = log.(sHat(x, X_1, sigma_1, zeta_1))
println(sHat_log_jac(delta_1, X_1, sigma_1, zeta_1))
ForwardDiff.jacobian(shl, delta_1)


