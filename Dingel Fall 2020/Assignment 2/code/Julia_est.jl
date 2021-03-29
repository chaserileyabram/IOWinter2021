using StatFiles
using DataFrames
using FixedEffectModels
using RegressionTables

df = DataFrame(load("col_regfile09.dta"))

df = df[df.flow .> 0,:]
df.logflow = log.(df.flow)
df.logdist = log.(df.distw)

small_df = df[df.year .> 1999,:]
small_reg = reg(small_df, @formula(logflow ~ logdist + contig 
        + comlang_off + fe(iso_o)&fe(year) + fe(iso_d)&fe(year)), Vcov.robust())

@time r = reg(df, @formula(logflow ~ logdist + contig 
        + comlang_off + fe(iso_o)&fe(year) + fe(iso_d)&fe(year)), 
    Vcov.robust(),
)

regtable(r, renderSettings = latexOutput("Table_3_Julia.tex"))
r
