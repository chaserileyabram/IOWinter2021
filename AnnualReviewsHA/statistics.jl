# Chase Abram
# For computing statistics

include("ModelDiscrete.jl")
include("stationarydist.jl")
using Statistics
using Pkg
# Pkg.add("QuadGK")
using QuadGK

# Mean over all types of income
function mean_income(m)
    println("mean_income: ", sum(m.statdist .* (m.yyFgrid + m.yyPgrid + m.yyTgrid)))
    return sum(m.statdist .* (m.yyFgrid + m.yyPgrid + m.yyTgrid))
end

# Variance in log of income (includes all types)
function var_log_gross_income(m)
    println("var_log_gross_income: ", sum(m.statdist .* log.(m.yyFgrid + m.yyPgrid + m.yyTgrid).^2) - sum(m.statdist .* log.(m.yyFgrid + m.yyPgrid + m.yyTgrid))^2)
    return sum(m.statdist .* log.(m.yyFgrid + m.yyPgrid + m.yyTgrid).^2) - sum(m.statdist .* log.(m.yyFgrid + m.yyPgrid + m.yyTgrid))^2
end

# Mean wealth
function mean_wealth(m)
    println("mean_wealth: ", sum(m.statdist .* m.aagrid))
    return sum(m.statdist .* m.aagrid)
end

# Prob. wealth <= a
function wealth_lesseq(m,a)
    println("wealth less than or equal to ", a, ": ", sum(m.statdist[m.aagrid .<= a]))
    return sum(m.statdist[m.aagrid .<= a])
end

# Prob. sav <= a
function a_tom_lesseq(m,a)
    println("savings less than or equal to ", a, ": ", sum(m.statdist[m.a_tom .<= a]))
    return sum(m.statdist[m.a_tom .<= a])
end

# Prob. wealth <= k*income
function wealth_lesseq_fi(m,f)
    # println("wealth less than or equal to ", f, "*income: ", sum(m.statdist[m.aagrid .<= (f .* (m.yyFgrid + m.yyPgrid + m.yyTgrid))]))
    return sum(m.statdist[m.aagrid .<= (f .* (m.yyFgrid + m.yyPgrid + m.yyTgrid))])
end

# Wealth level at quantile q
function wealth_quantile(m,q)
    if q <= wealth_lesseq(m,m.agrid[1])
        wbc = wealth_lesseq(m,m.agrid[1])
        # println("Requested ", q,", but ", wbc, " at borrowing constraint")
        return wbc
    end

    m.tmp_itp = interpolate((cumsum(sum(m.statdist, dims=(2,3,4,5))[:,1,1,1,1]),) , m.agrid, Gridded(Linear()))
    # println("wealth at quantile ", q, ": ", m.tmp_itp(q))
    return m.tmp_itp(q)
end

# Wealth held by top q
function wealth_held(m, q)
    # println("Top ", q, " of wealthholders hold proportion: ", sum((m.statdist .* m.aagrid)[m.aagrid .> wealth_quantile(m,1-q)])/sum(m.statdist .* m.aagrid))
    return sum((m.statdist .* m.aagrid)[m.aagrid .> wealth_quantile(m,1-q)])/sum(m.statdist .* m.aagrid)
end

function gini(m)
    # println("gini: ", 1 - sum(cumsum(sum(m.statdist, dims=(2,3,4,5))[:,1,1,1,1])))
    # return 1 - 2*sum(cumsum(sum(m.statdist, dims=(2,3,4,5))[:,1,1,1,1]))

    lorenz(z) = 1 - wealth_held(m,z)
    return 1 - 2*quadgk(lorenz, 0, 1, rtol=1e-6)
end


# m0 = ModelDiscrete()
# setup_power_grids(m0)
# setup_income(m0)


mean_income(m0)

var_log_gross_income(m0)

mean_wealth(m0)

wealth_lesseq(m0,0.0)
wealth_lesseq(m0,1.0)

a_tom_lesseq(m0, 0.0)
a_tom_lesseq(m0, 1.0)

wealth_lesseq_fi(m0,0)
wealth_lesseq_fi(m0,1/6)
wealth_lesseq_fi(m0,1000)

wealth_quantile(m0, 0.0)
wealth_quantile(m0, 0.1)
wealth_quantile(m0, 0.25)
wealth_quantile(m0, 0.5)
wealth_quantile(m0, 0.9)
wealth_quantile(m0, 0.99)
wealth_quantile(m0, 0.999)

wealth_held(m0, 1)
wealth_held(m0, 0.5)
wealth_held(m0, 0.1)
wealth_held(m0, 0.01)

